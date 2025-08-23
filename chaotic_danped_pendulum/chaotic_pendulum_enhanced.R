# chaotic_pendulum.R <- an enhanced V2 with 0-1 chaos tests, Poincare sections etc.
# Periodically forced, damped pendulum: RK4 integrator + phase portraits,
# Poincaré sections, bifurcation diagram vs A, and approximate LLE.
#
# Minimal deps: ggplot2, dplyr
suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("dplyr", quietly = TRUE);   if (!ok) stop("Install dplyr")
})
library(ggplot2)
library(dplyr)

# ---------------- Utilities ----------------

wrap_to_pi <- function(x) {
  # Map angle to (-pi, pi]
  y <- (x + pi) %% (2*pi)
  y[y == 0] <- 2*pi
  y - pi
}

unwrap_angle <- function(theta_wrapped) {
  # Unwrap a sequence of angles (radians) to a continuous signal
  # Assumes small step-to-step changes (true for small dt).
  unwrapped <- theta_wrapped
  for (i in 2:length(unwrapped)) {
    delta <- unwrapped[i] - unwrapped[i-1]
    if (delta >  pi) unwrapped[i:length(unwrapped)] <- unwrapped[i:length(unwrapped)] - 2*pi
    if (delta < -pi) unwrapped[i:length(unwrapped)] <- unwrapped[i:length(unwrapped)] + 2*pi
  }
  unwrapped
}

# ---------------- Model ----------------

pend_rhs <- function(t, x, pars) {
  # x = c(theta, omega), where omega = dtheta/dt
  theta <- x[1]; omega <- x[2]
  gamma  <- pars$gamma
  A      <- pars$A
  wdrive <- pars$omega
  
  dtheta <- omega
  domega <- -gamma*omega - sin(theta) + A*cos(wdrive*t)
  c(dtheta, domega)
}

rk4_step <- function(f, t, x, dt, pars) {
  k1 <- f(t,         x,            pars)
  k2 <- f(t + dt/2,  x + dt*k1/2,  pars)
  k3 <- f(t + dt/2,  x + dt*k2/2,  pars)
  k4 <- f(t + dt,    x + dt*k3,    pars)
  x + dt*(k1 + 2*k2 + 2*k3 + k4)/6
}

integrate_pendulum <- function(theta0=0.2, omega0=0.0,
                               pars=list(gamma=0.2, A=1.2, omega=2/3),
                               t_max=200, dt=0.01) {
  n  <- floor(t_max/dt) + 1
  t  <- seq(0, by=dt, length.out=n)
  x  <- matrix(0.0, nrow=n, ncol=2)
  x[1,] <- c(theta0, omega0)
  
  for (i in 2:n) {
    x[i,] <- rk4_step(pend_rhs, t[i-1], x[i-1,], dt, pars)
  }
  
  theta <- x[,1]
  omega <- x[,2]
  tibble::tibble(
    t = t,
    theta_unwrapped = theta,
    theta = wrap_to_pi(theta),
    omega = omega
  )
}

# Linear interpolation of state at arbitrary times
interp_state <- function(sim_df, t_query) {
  # Interpolate unwrapped theta and omega at t_query
  th <- approx(sim_df$t, sim_df$theta_unwrapped, xout=t_query, rule=2)$y
  om <- approx(sim_df$t, sim_df$omega,           xout=t_query, rule=2)$y
  tibble::tibble(
    t = t_query,
    theta_unwrapped = th,
    theta = wrap_to_pi(th),
    omega = om
  )
}

# ---------------- Diagnostics & Plots ----------------

plot_time_series <- function(sim_df, n_keep = NULL) {
  df <- sim_df
  if (!is.null(n_keep)) df <- df %>% dplyr::slice_tail(n=n_keep)
  p1 <- ggplot(df, aes(t, theta)) + geom_line() + labs(title="theta(t)", y="theta (wrapped)")
  p2 <- ggplot(df, aes(t, omega)) + geom_line() + labs(title="omega(t)", y="d theta / dt")
  list(p_theta = p1, p_omega = p2)
}

plot_phase <- function(sim_df, frac_discard=0.5) {
  n <- nrow(sim_df)
  start <- floor(frac_discard*n) + 1
  df <- sim_df[start:n,]
  ggplot(df, aes(theta, omega)) +
    geom_path(alpha=0.6) +
    labs(title="Phase portrait", x="theta (wrapped to [-pi,pi])", y="omega")
}

poincare_section <- function(sim_df, pars, n_skip_periods=200, n_take_periods=200) {
  Tdrive <- 2*pi/pars$omega
  t0 <- sim_df$t[1]
  t_skip <- t0 + n_skip_periods*Tdrive
  t_take <- t_skip + (0:(n_take_periods-1))*Tdrive
  interp_state(sim_df, t_take)
}

plot_poincare <- function(pmap_df) {
  ggplot(pmap_df, aes(theta, omega)) +
    geom_point(alpha=0.6, size=0.8) +
    labs(title="Poincaré section (strobing at drive period)", x="theta (wrapped)", y="omega")
}

# ---------------- Bifurcation vs A + LLE ----------------

# Approximate largest Lyapunov exponent via two-trajectory method (Benettin-style)
lyapunov_largest <- function(theta0=0.2, omega0=0.0,
                             pars=list(gamma=0.2, A=1.2, omega=2/3),
                             dt=0.01,
                             n_periods_total=800,
                             renorm_every_steps=50,
                             d0=1e-8) {
  Tdrive <- 2*pi/pars$omega
  n_steps <- floor((n_periods_total*Tdrive)/dt)
  
  x  <- c(theta0, omega0)
  # small random perturbation direction
  dir <- rnorm(2); dir <- dir / sqrt(sum(dir^2))
  y  <- x + d0*dir
  
  sum_log <- 0.0
  t <- 0.0
  k <- 0L
  
  for (i in 1:n_steps) {
    x <- rk4_step(pend_rhs, t, x, dt, pars)
    y <- rk4_step(pend_rhs, t, y, dt, pars)
    t <- t + dt
    
    if (i %% renorm_every_steps == 0) {
      delta <- y - x
      # keep theta difference continuous (states are unwrapped internally)
      d <- sqrt(sum(delta^2))
      if (d == 0) {
        # avoid log(0)
        d <- .Machine$double.eps
      }
      sum_log <- sum_log + log(d/d0)
      # renormalize
      y <- x + (delta / d) * d0
      k <- k + 1L
    }
  }
  # total time of k renormalizations
  tau <- k * renorm_every_steps * dt
  if (tau == 0) return(NA_real_)
  sum_log / tau
}

bifurcation_vs_A <- function(A_seq,
                             base_pars=list(gamma=0.2, A=1.2, omega=2/3),
                             theta0=0.2, omega0=0.0,
                             dt_factor=200,     # steps per drive period
                             n_skip_periods=300,
                             n_take_periods=120,
                             compute_lle=TRUE,
                             lle_periods=600) {
  out_pts <- list()
  lle_vals <- numeric(length(A_seq))
  names(lle_vals) <- as.character(A_seq)
  
  for (idx in seq_along(A_seq)) {
    A <- A_seq[idx]
    pars <- base_pars; pars$A <- A
    
    Tdrive <- 2*pi/pars$omega
    dt <- Tdrive / dt_factor
    t_max <- (n_skip_periods + n_take_periods + 5) * Tdrive # +5 buffer
    
    sim <- integrate_pendulum(theta0, omega0, pars, t_max=t_max, dt=dt)
    pmap <- poincare_section(sim, pars, n_skip_periods=n_skip_periods, n_take_periods=n_take_periods)
    out_pts[[idx]] <- tibble::tibble(A=A, theta=pmap$theta, omega=pmap$omega)
    
    if (compute_lle) {
      lle_vals[idx] <- lyapunov_largest(theta0, omega0, pars,
                                        dt = dt,
                                        n_periods_total = lle_periods)
    }
  }
  bif <- dplyr::bind_rows(out_pts)
  if (compute_lle) {
    lle_df <- tibble::tibble(A=A_seq, LLE=lle_vals)
  } else {
    lle_df <- NULL
  }
  list(bif=bif, lle=lle_df)
}

plot_bifurcation <- function(bif_df, yvar=c("theta","omega")) {
  yvar <- match.arg(yvar)
  ggplot(bif_df, aes(A, .data[[yvar]])) +
    geom_point(alpha=0.3, size=0.3) +
    labs(title=paste("Bifurcation diagram (y =", yvar, ")"),
         x="Forcing amplitude A", y=yvar)
}

plot_lle <- function(lle_df) {
  ggplot(lle_df, aes(A, LLE)) +
    geom_line() + geom_hline(yintercept=0, linetype="dashed") +
    labs(title="Largest Lyapunov exponent vs A", x="A", y="LLE")
}

# ---------------- Convenience demos ----------------

demo_phase <- function(pars=list(gamma=0.2, A=1.2, omega=2/3),
                       theta0=0.2, omega0=0.0,
                       n_periods=400, dt_factor=200) {
  Tdrive <- 2*pi/pars$omega
  dt <- Tdrive / dt_factor
  sim <- integrate_pendulum(theta0, omega0, pars, t_max=n_periods*Tdrive, dt=dt)
  
  plots <- plot_time_series(sim)
  print(plots$p_theta); print(plots$p_omega)
  Sys.sleep(5)
  print(plot_phase(sim, frac_discard=0.5))
  Sys.sleep(5)
  pmap <- poincare_section(sim, pars, n_skip_periods=200, n_take_periods=200)
  print(plot_poincare(pmap))
  Sys.sleep(5)
  invisible(list(sim=sim, pmap=pmap))
}

demo_bifurcation <- function(A_min=0.9, A_max=1.6, nA=120,
                             base_pars=list(gamma=0.2, A=1.2, omega=2/3),
                             yvar="theta",
                             compute_lle=TRUE) {
  A_seq <- seq(A_min, A_max, length.out=nA)
  res <- bifurcation_vs_A(
    A_seq,
    base_pars = base_pars,
    dt_factor = 200,
    n_skip_periods = 250,
    n_take_periods = 120,
    compute_lle = compute_lle,
    lle_periods = 500
  )
  print(plot_bifurcation(res$bif, yvar=yvar))
  Sys.sleep(10)
  if (!is.null(res$lle)) print(plot_lle(res$lle))
  invisible(res)
}

# ---------------- 0–1 test for chaos ---------------- (Gottwald - Melbourne test)
zero_one_test <- function(sim_df,
                          observable=c("theta","omega"),
                          frac_discard=0.5,
                          c_vals=NULL) {
  observable <- match.arg(observable)
  n <- nrow(sim_df)
  start <- floor(frac_discard*n) + 1
  v <- sim_df[[observable]][start:n]
  v <- as.numeric(v - mean(v))  # de-mean
  N <- length(v)
  if (is.null(c_vals)) {
    # Avoid resonances near 0, pi; pick a spread of c in (0, pi)
    set.seed(42)
    c_vals <- runif(25, min=0.1, max=pi-0.1)
  }
  
  Ks <- numeric(length(c_vals))
  names(Ks) <- sprintf("%.3f", c_vals)
  
  idx <- 1:N
  for (i in seq_along(c_vals)) {
    c <- c_vals[i]
    pc <- cumsum(v * cos(idx * c))
    qc <- cumsum(v * sin(idx * c))
    Mc <- pc*pc + qc*qc
    # Correlation of mean-square displacement with n: ~1 for chaos, ~0 for order
    Ks[i] <- suppressWarnings(cor(Mc, idx, method="pearson"))
  }
  K <- median(Ks, na.rm=TRUE)
  list(K=K, Ks=Ks, N=N)
}

demo_zero_one <- function(pars=list(gamma=0.2, A=1.2, omega=2/3),
                          theta0=0.2, omega0=0.0,
                          n_periods=800, dt_factor=200,
                          observable="theta") {
  Tdrive <- 2*pi/pars$omega
  dt <- Tdrive / dt_factor
  sim <- integrate_pendulum(theta0, omega0, pars, t_max=n_periods*Tdrive, dt=dt)
  z <- zero_one_test(sim, observable=observable, frac_discard=0.5)
  message(sprintf("0–1 test K ≈ %.3f (≈1 chaos, ≈0 regular) [N=%d]", z$K, z$N))
  invisible(z)
}

# Return Poincare section maps - super useful to look at ##

# ---------------- Return map on the Poincaré section ----------------
return_map <- function(pmap_df, var=c("theta","omega")) {
  var <- match.arg(var)
  v <- pmap_df[[var]]
  if (length(v) < 3) stop("Poincaré series too short for a return map.")
  df <- tibble::tibble(x = v[-length(v)], y = v[-1])
  ggplot(df, aes(x, y)) +
    geom_point(alpha=0.6, size=0.8) +
    labs(title=paste("Return map on Poincaré: ", var[nchar(var)>0]),
         x=paste0(var, "[n]"), y=paste0(var, "[n+1]"))
}

demo_return_map <- function(pars=list(gamma=0.2, A=1.2, omega=2/3),
                            theta0=0.2, omega0=0.0,
                            n_skip_periods=200, n_take_periods=400,
                            dt_factor=200, var="theta") {
  Tdrive <- 2*pi/pars$omega
  dt <- Tdrive / dt_factor
  sim <- integrate_pendulum(theta0, omega0, pars, t_max=(n_skip_periods+n_take_periods+10)*Tdrive, dt=dt)
  pmap <- poincare_section(sim, pars, n_skip_periods=n_skip_periods, n_take_periods=n_take_periods)
  print(plot_poincare(pmap))
  Sys.sleep(5)
  print(return_map(pmap, var=var))
  Sys.sleep(5)
  invisible(pmap)
}

# basins of attraction - important for characterising robustness of the attractor, fixed point 
# can be quantified via "area" / "generic volume" / volume fractions in N-D 

# ---------------- Basins of attraction over (theta0, omega0) grid ----------------
classify_attractor <- function(pmap_df, last_k=6, round_digits=2) {
  if (nrow(pmap_df) < last_k) return(NA_character_)
  tail_df <- pmap_df %>% dplyr::slice_tail(n=last_k)
  # Average last_k points to stabilize label
  th <- mean(tail_df$theta)
  om <- mean(tail_df$omega)
  paste0(round(th, round_digits), "_", round(om, round_digits))
}

basin_of_attraction_grid <- function(theta_range=c(-pi, pi),
                                     omega_range=c(-2, 2),
                                     n_theta=80, n_omega=80,
                                     pars=list(gamma=0.2, A=1.2, omega=2/3),
                                     dt_factor=160,
                                     n_skip_periods=200,
                                     n_take_periods=80,
                                     last_k=6,
                                     round_digits=2,
                                     quiet=FALSE) {
  thetas0 <- seq(theta_range[1], theta_range[2], length.out=n_theta)
  omegas0 <- seq(omega_range[1], omega_range[2], length.out=n_omega)
  
  Tdrive <- 2*pi/pars$omega
  dt <- Tdrive / dt_factor
  t_max <- (n_skip_periods + n_take_periods + 6) * Tdrive
  
  out <- vector("list", length= n_theta * n_omega)
  idx <- 1L
  
  for (i in seq_along(thetas0)) {
    for (j in seq_along(omegas0)) {
      th0 <- thetas0[i]; om0 <- omegas0[j]
      sim <- integrate_pendulum(theta0=th0, omega0=om0, pars=pars, t_max=t_max, dt=dt)
      pmap <- poincare_section(sim, pars, n_skip_periods=n_skip_periods, n_take_periods=n_take_periods)
      lab <- classify_attractor(pmap, last_k=last_k, round_digits=round_digits)
      out[[idx]] <- tibble::tibble(theta0=th0, omega0=om0, class_id=lab)
      idx <- idx + 1L
    }
    if (!quiet) message(sprintf("Row %d/%d done", i, n_theta))
  }
  dplyr::bind_rows(out)
}

plot_basin <- function(basin_df) {
  basin_df$class_id <- factor(basin_df$class_id)
  ggplot(basin_df, aes(x=theta0, y=omega0, fill=class_id)) +
    geom_tile() +
    guides(fill="none") +
    labs(title="Basins of attraction (coarse classes from last Poincaré points)",
         x=expression(theta[0]), y=expression(omega[0]))
}

demo_basin <- function(pars=list(gamma=0.2, A=1.35, omega=2/3),
                       theta_range=c(-pi, pi), omega_range=c(-2, 2),
                       n_theta=80, n_omega=80,
                       dt_factor=160,
                       n_skip_periods=220, n_take_periods=90,
                       last_k=6, round_digits=2) {
  basin <- basin_of_attraction_grid(theta_range, omega_range, n_theta, n_omega,
                                    pars, dt_factor, n_skip_periods, n_take_periods,
                                    last_k, round_digits, quiet=FALSE)
  print(plot_basin(basin))
  invisible(basin)
}



# ---------------- How to run (examples) ----------------
# From R:
demo_phase()                           # single-run: time series, phase portrait, Poincaré
Sys.sleep(5)
demo_bifurcation()   # sweep A, bifurcation + LLE overlay
Sys.sleep(5)
demo_zero_one(pars=list(gamma=0.2, A=1.4, omega=2/3), observable="theta")
Sys.sleep(5)
demo_return_map(pars=list(gamma=0.22, A=1.45, omega=2/3), var="theta")
Sys.sleep(5)
# Chaos-prone parameters show gnarly boundaries:
demo_basin(pars=list(gamma=0.2, A=1.4, omega=2/3),
           theta_range=c(-pi, pi), omega_range=c(-2, 2),
           n_theta=80, n_omega=80)

# personalise #
# Tweak params, e.g.:
#   demo_phase(pars=list(gamma=0.25, A=1.4, omega=2/3))
#   demo_bifurcation(A_min=0.8, A_max=1.8, nA=160,
#                    base_pars=list(gamma=0.22, A=1.2, omega=2/3), yvar="omega")
