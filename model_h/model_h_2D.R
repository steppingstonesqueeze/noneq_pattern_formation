# model_h_2d_fd.R  (robust V2)
# 2D Model H: CH + incompressible Navier–Stokes (periodic), real-space FD, semi-implicit.
# Elliptic solves (CH Helmholtz-biharmonic, velocity Helmholtz, pressure Poisson) via NA-safe matrix-free CG.
# Viz: phi raster + downsampled quiver of velocity; MP4 via av.
# Metrics: mass, CH free energy, kinetic energy, length scale L ~ sqrt(∫phi^2 / ∫|∇phi|^2).

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# ---------------- User params ----------------
nx <- 192; ny <- 192              # grid (even)
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny          # assume dx=dy
dt <- 0.01                        # small for advective CFL
nsteps <- 2000
save_every    <- 10               # frames cadence
metrics_every <- 10

# CH params
epsilon <- 0.02
Mmob    <- 1.0
m0      <- 0.0                    # target mean composition (conserved)
ic_amp  <- 0.10                   # random IC amplitude

# Fluid params
rho <- 1.0
nu  <- 0.5                        # kinematic viscosity (damps fast)

# CG tolerances (looser early, tighter late)
tol_initial <- 2e-6
tol_final   <- 1e-8

# Soft caps (safety nets; raise if you find them too tight)
cap_phi <- 5      # |phi|
cap_u   <- 50     # |ux|, |uy|
cap_p   <- 1e3    # |p|
# ---------------------------------------------

set.seed(7)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# ---------- periodic shifts ----------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

# ---------- derivatives ----------
D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap  <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

div2D <- function(ux, uy) D1x(ux) + D1y(uy)
gradp <- function(p) list(px = D1x(p), py = D1y(p))

# ---------- helpers ----------
sanitize <- function(U, cap = NA_real_) {
  U[!is.finite(U)] <- 0
  if (is.finite(cap)) {
    U[U >  cap] <-  cap
    U[U < -cap] <- -cap
  }
  U
}

# ---------- NA-safe matrix-free CG ----------
cg_solve <- function(apply_A, b, x0 = NULL, tol = 1e-8, maxit = 300, proj_zero_mean = FALSE) {
  safe_dot <- function(a, b) {
    s <- sum(as.numeric(a) * as.numeric(b), na.rm = TRUE)
    if (!is.finite(s)) 0 else s
  }
  make_finite <- function(M) { M[!is.finite(M)] <- 0; M }
  
  x <- if (is.null(x0)) matrix(0, nrow = nrow(b), ncol = ncol(b)) else x0
  if (proj_zero_mean) b <- b - mean(b)
  b <- make_finite(b)
  
  Ax <- make_finite(apply_A(x))
  r  <- make_finite(b - Ax)
  if (proj_zero_mean) r <- r - mean(r)
  p  <- r
  
  rr <- safe_dot(r, r)
  bnorm <- sqrt(max(safe_dot(b, b), 0)) + 1e-30
  
  for (it in 1:maxit) {
    Ap <- make_finite(apply_A(p))
    denom <- safe_dot(p, Ap) + 1e-30
    if (!is.finite(denom) || denom <= 0) {
      warning("CG: non-positive or non-finite denominator; returning current iterate")
      return(make_finite(x))
    }
    alpha <- rr / denom; if (!is.finite(alpha)) alpha <- 0
    
    x <- make_finite(x + alpha * p)
    r <- make_finite(r - alpha * Ap)
    if (proj_zero_mean) r <- r - mean(r)
    
    rr_new <- safe_dot(r, r)
    if (!is.finite(rr_new)) {
      warning("CG: non-finite residual; returning current iterate")
      return(make_finite(x))
    }
    if (sqrt(rr_new) <= tol * bnorm) return(make_finite(x))
    
    beta <- rr_new / (rr + 1e-30); if (!is.finite(beta)) beta <- 0
    p  <- make_finite(r + beta * p)
    rr <- rr_new
  }
  make_finite(x)
}

# ---------- init fields ----------
phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
ux  <- matrix(0, ny, nx)
uy  <- matrix(0, ny, nx)
p   <- matrix(0, ny, nx)

# ---------- viz helpers ----------
frames_dir <- "modelh_frames"
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L
quiv_stride <- 6  # downsample quiver grid

save_frame <- function(phi, ux, uy, ttag) {
  # scalar field
  df <- as.data.frame(as.table(phi)); colnames(df) <- c("y","x","val")
  rng <- range(df$val[is.finite(df$val)])
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(diff(rng)) < .Machine$double.eps) {
    df$val <- df$val + rnorm(nrow(df), sd = 1e-9)
  }
  p1 <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"), na.value = "grey80") +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("Model H: composition φ  (t=%.2f)", ttag))
  
  # quiver: speed background + arrows
  speed <- sqrt(ux^2 + uy^2)
  speed <- sanitize(speed)
  df2 <- as.data.frame(as.table(speed)); colnames(df2) <- c("y","x","spd")
  xs <- seq(1, ncol(ux), by = quiv_stride); ys <- seq(1, nrow(ux), by = quiv_stride)
  qx <- ux[ys, xs]; qy <- uy[ys, xs]
  qdf <- expand.grid(x = xs, y = ys); qdf$ux <- as.numeric(qx); qdf$uy <- as.numeric(qy)
  scl <- 0.6 * quiv_stride  # visual scale
  p2 <- ggplot(df2, aes(x, y, fill = spd)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#f7fbff","#c6dbef","#6baed6","#2171b5")) +
    geom_segment(data = qdf,
                 aes(x = x, y = y, xend = x + scl*ux, yend = y + scl*uy),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black",
                 arrow = arrow(length = grid::unit(0.08, "cm"))) +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("velocity |u| and arrows (t=%.2f)", ttag))
  
  fn1 <- file.path(frames_dir, sprintf("frame_phi_%05d.png", frame_id))
  fn2 <- file.path(frames_dir, sprintf("frame_u_%05d.png", frame_id))
  ggsave(fn1, p1, width = 5, height = 5, dpi = 120)
  ggsave(fn2, p2, width = 5, height = 5, dpi = 120)
}

encode_video <- function() {
  files1 <- list.files(frames_dir, pattern = "frame_phi_\\d+\\.png$", full.names = TRUE)
  files2 <- list.files(frames_dir, pattern = "frame_u_\\d+\\.png$",   full.names = TRUE)
  files1 <- files1[order(files1)]; files2 <- files2[order(files2)]
  if (length(files1)) av::av_encode_video(files1, output = "modelh_phi.mp4", framerate = 24)
  if (length(files2)) av::av_encode_video(files2, output = "modelh_u.mp4",   framerate = 24)
}

# ---------- metrics ----------
free_energy_CH <- function(phi) {
  g <- 0.25*(phi^2 - 1)^2
  e <- 0.5*epsilon^2 * (D1x(phi)^2 + D1y(phi)^2) + g
  sum(e) * dx * dy
}
length_scale <- function(phi) {
  num <- sum(phi^2) * dx * dy
  den <- sum(D1x(phi)^2 + D1y(phi)^2) * dx * dy + 1e-30
  sqrt(num / den)
}

ts <- data.frame(t = numeric(0), mass = numeric(0), F = numeric(0),
                 L = numeric(0), KE = numeric(0))

message(sprintf("IC: mass=%.4f, L=%.3f", mean(phi), length_scale(phi)))

# ---------- time loop ----------
for (step in 0:nsteps) {
  tol_now <- if (step < nsteps/2) tol_initial else tol_final
  
  # ---- CH step (stabilized IMEX) ----
  # (I - dt M ∇² + dt M ε² ∇⁴) φ^{n+1} = φ^n + dt [ -u·∇φ^n + M ∇²(φ^3)^n ]
  adv_term <- -(ux*D1x(phi) + uy*D1y(phi))                 # explicit advection
  rhs_phi  <- phi + dt * ( adv_term + Mmob * lap(phi^3) )  # no double -lap(phi)
  A_phi    <- function(U) U - dt*Mmob*lap(U) + dt*Mmob*epsilon^2*bilap(U)
  phi <- cg_solve(A_phi, b = rhs_phi, x0 = phi, tol = tol_now, maxit = 250)
  phi <- sanitize(phi, cap = cap_phi)
  # Mass recentering (numerical conservation)
  phi <- phi + (m0 - mean(phi))
  
  # ---- chemical potential for capillary force ----
  mu <- -epsilon^2 * lap(phi) + (phi^3 - phi)
  
  # ---- momentum: semi-implicit viscous + projection ----
  # explicit advection + capillary force
  advx <- ux*D1x(ux) + uy*D1y(ux)
  advy <- ux*D1x(uy) + uy*D1y(uy)
  fx <- -(phi)*D1x(mu) / rho
  fy <- -(phi)*D1y(mu) / rho
  
  rhs_ux <- ux + dt * ( -advx + fx )
  rhs_uy <- uy + dt * ( -advy + fy )
  
  # Helmholtz solve (I - dt nu ∇^2) u* = rhs
  A_u <- function(U) U - dt*nu*lap(U)
  ux_star <- cg_solve(A_u, b = rhs_ux, x0 = ux, tol = tol_now, maxit = 200)
  uy_star <- cg_solve(A_u, b = rhs_uy, x0 = uy, tol = tol_now, maxit = 200)
  ux_star <- sanitize(ux_star, cap = cap_u)
  uy_star <- sanitize(uy_star, cap = cap_u)
  
  # Projection: solve ∇^2 p = (rho/dt) ∇·u*  (zero-mean RHS)
  div_star <- div2D(ux_star, uy_star)
  div_star <- div_star - mean(div_star)
  p <- cg_solve(function(U) -lap(U), b = (rho/dt)*div_star, x0 = p, tol = tol_now, maxit = 300, proj_zero_mean = TRUE)
  p <- sanitize(p, cap = cap_p)
  
  # Update velocity
  gp <- gradp(p)
  ux <- ux_star - (dt/rho) * gp$px
  uy <- uy_star - (dt/rho) * gp$py
  ux <- sanitize(ux, cap = cap_u)
  uy <- sanitize(uy, cap = cap_u)
  
  # ---- output & metrics ----
  if (step %% metrics_every == 0) {
    mass <- mean(phi)
    F    <- free_energy_CH(phi)
    Lc   <- length_scale(phi)
    KE   <- 0.5 * sum(ux^2 + uy^2) * dx * dy
    ts <- rbind(ts, data.frame(t = step*dt, mass = mass, F = F, L = Lc, KE = KE))
    message(sprintf("t=%.2f  mass=%.4f  F=%.3e  L=%.3f  KE=%.3e",
                    step*dt, mass, F, Lc, KE))
  }
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(phi, ux, uy, ttag = step*dt)
  }
}

write.csv(ts, "modelh_timeseries.csv", row.names = FALSE)
encode_video()
message("Done: Model H (FD, robust V2)")
