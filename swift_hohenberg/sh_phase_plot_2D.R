# sh_phase_plot_fd.R — Swift–Hohenberg 2D phase map (FD semi-implicit + CG solver)
# Sweeps (r,g), simulates with a finite-difference semi-implicit scheme (no FFT in solver),
# computes spectral diagnostics (single FFT for analysis), classifies pattern types,
# and saves a phase map + CSV + thumbnails.
#
# Run: Rscript sh_phase_plot_fd.R
# Requires: ggplot2, av (av only if you later enable sample frame saves)

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
})
library(ggplot2)

# ----------------------- USER CONFIG -----------------------
# Parameter grid
RG_R <- c(0.05, 0.15, 0.25)     # linear growth parameter r
RG_G <- c(0.6, 1.0, 1.4)        # cubic coefficient g

RG_R <- seq(0.05, 0.55, by = 0.1)
RG_G <- seq(0.6, 1, by = 0.05)

# Simulation params (keep modest for sweep)
nx <- 192; ny <- 192            # even
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx / nx; dy <- Ly / ny    # assume dx=dy
dt <- 0.05
steps <- 1500                  # steps per (r,g)
eta <- 0.05                     # additive noise amplitude (eta*sqrt(dt)*N)
seed <- 42

# Optional: save intermediate PNGs per (r,g)? (slower)
save_samples <- FALSE
sample_every <- 100

# Ring detection / angular spectrum
ring_halfwidth_frac <- 0.08     # ±8% around k* for angle spectrum
n_theta_bins <- 90

# Classification thresholds (tweak as you like)
class_rules <- list(
  stripe_S2_thresh = 0.40,
  stripe_H6_max    = 0.25,
  hex_H6_min       = 0.30,
  iso_A_max        = 0.15,
  spots_skew_min   = 0.20       # |skew| threshold
)

# Output
csv_out <- "sh_phase_metrics.csv"
phase_png <- "sh_phase_map.png"
thumb_dir <- "thumbs"
if (!dir.exists(thumb_dir)) dir.create(thumb_dir, recursive = TRUE)

set.seed(seed)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# -------------------- PERIODIC FD OPERATORS ---------------------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# Implicit SPD operator for SH step: A(U) = U + dt * ∇^4 U
apply_A <- function(U) U + dt * bilap(U)

# Matrix-free Conjugate Gradient (2D)
cg_solve <- function(b, x0 = NULL, tol = 1e-8, maxit = 200) {
  x <- if (is.null(x0)) matrix(0, nrow = ny, ncol = nx) else x0
  r <- b - apply_A(x)
  p <- r
  rr <- sum(r*r)
  bnorm <- sqrt(sum(b*b)) + 1e-30
  for (it in 1:maxit) {
    Ap <- apply_A(p)
    alpha <- rr / (sum(p*Ap) + 1e-30)
    x <- x + alpha * p
    r <- r - alpha * Ap
    rr_new <- sum(r*r)
    if (sqrt(rr_new) <= tol * bnorm) break
    beta <- rr_new / (rr + 1e-30)
    p <- r + beta * p
    rr <- rr_new
  }
  x
}

# ------------------ SH TIME-STEPPER (FD, no FFT) ----------------
simulate_SH_fd <- function(rparam, gparam, steps, dt, u0 = NULL, save_png = FALSE, tag = "") {
  u <- if (is.null(u0)) matrix(0.01 * rnorm(nx*ny), nrow = ny, ncol = nx) else u0
  frames_dir <- NULL
  if (save_png) {
    frames_dir <- sprintf("sh_fd_frames_r%.3f_g%.3f", rparam, gparam)
    if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
  }
  
  if (save_png) {
    df <- as.data.frame(as.table(u)); colnames(df) <- c("y","x","val")
    p <- ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
      scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
      ggtitle(sprintf("t=%.2f  (r=%.2f, g=%.2f)", 0.0, rparam, gparam))
    ggsave(file.path(frames_dir, sprintf("frame_%05d.png", 0)), p, width = 4.5, height = 4.5, dpi = 120)
  }
  
  for (s in 1:steps) {
    # Explicit RHS: r u - 2 ∇^2 u - g u^3 + noise
    rhs <- u + dt * ( rparam * u - 2 * lap(u) - gparam * (u^3) ) +
      eta * sqrt(dt) * matrix(rnorm(nx*ny), nrow = ny, ncol = nx)
    
    # Implicit solve: (I + dt ∇^4) u^{n+1} = rhs
    u <- cg_solve(rhs, x0 = u, tol = 1e-8, maxit = 200)
    
    # Guards
    if (!all(is.finite(u))) u[!is.finite(u)] <- 0
    rng_now <- range(u[is.finite(u)])
    if (is.finite(rng_now[1]) && is.finite(rng_now[2]) && abs(rng_now[2] - rng_now[1]) < 1e-14) {
      u <- u + 1e-6 * matrix(rnorm(length(u)), nrow = ny, ncol = nx)
    }
    
    if (save_png && (s %% sample_every == 0)) {
      df <- as.data.frame(as.table(u)); colnames(df) <- c("y","x","val")
      p <- ggplot(df, aes(x = x, y = y, fill = val)) +
        geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
        scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
        ggtitle(sprintf("t=%.2f  (r=%.2f, g=%.2f)", s*dt, rparam, gparam))
      ggsave(file.path(frames_dir, sprintf("frame_%05d.png", s)), p, width = 4.5, height = 4.5, dpi = 120)
    }
  }
  u
}

# -------------------- k-GRID FOR DIAGNOSTICS --------------------
# We'll use FFT once per case to get S(k). This is for analysis only.
kx_vec <- (2*pi/Lx) * c(0:(nx/2-1), -nx/2:-1)
ky_vec <- (2*pi/Ly) * c(0:(ny/2-1), -ny/2:-1)
kx <- matrix(rep(kx_vec, each = ny), nrow = ny, ncol = nx)
ky <- matrix(rep(ky_vec, times = nx), nrow = ny, ncol = nx)
k2 <- kx^2 + ky^2

fft_rows <- function(u, inv = FALSE) t(apply(u, 1, function(x) fft(x, inverse = inv)))
fft_cols <- function(u, inv = FALSE) apply(u, 2, function(x) fft(x, inverse = inv))
fft2  <- function(u) fft_cols(fft_rows(u, inv = FALSE), inv = FALSE)

# ------------------ DIAGNOSTIC HELPERS --------------------
radial_avg <- function(S) {
  dkx <- 2*pi/Lx; dky <- 2*pi/Ly; dk <- min(dkx, dky)
  r <- sqrt(k2); bin <- floor(r/dk + 0.5)
  df <- data.frame(bin = as.vector(bin), S = as.vector(S))
  agg <- aggregate(S ~ bin, data = df, FUN = mean)
  k <- agg$bin * dk
  list(k = k, S = agg$S)
}

dominant_k <- function(S) {
  rs <- radial_avg(S)
  if (length(rs$k) < 3) return(list(kstar = NA_real_, idx = NA_integer_, k = rs$k, Sr = rs$S))
  pos <- which(rs$k > 0)
  if (length(pos) == 0) return(list(kstar = NA_real_, idx = NA_integer_, k = rs$k, Sr = rs$S))
  i <- pos[which.max(rs$S[pos])]
  list(kstar = rs$k[i], idx = i, k = rs$k, Sr = rs$S)
}

angular_spectrum <- function(Uhat, kstar, halfwidth_frac = 0.08, nbin = 90) {
  if (!is.finite(kstar) || kstar <= 0) return(list(theta = numeric(0), F = numeric(0)))
  r <- sqrt(k2)
  ring_lo <- (1 - halfwidth_frac) * kstar
  ring_hi <- (1 + halfwidth_frac) * kstar
  sel <- (r >= ring_lo) & (r <= ring_hi)
  if (!any(sel)) return(list(theta = numeric(0), F = numeric(0)))
  
  Uh <- Uhat[sel]; kx_s <- kx[sel]; ky_s <- ky[sel]
  pow <- Mod(Uh)^2
  theta <- atan2(ky_s, kx_s)  # (-pi, pi]
  edges <- seq(-pi, pi, length.out = nbin + 1)
  bins <- cut(theta, breaks = edges, include.lowest = TRUE, labels = FALSE)
  Ftheta <- numeric(nbin)
  for (i in seq_len(nbin)) Ftheta[i] <- sum(pow[bins == i])
  th_centers <- 0.5 * (edges[-1] + edges[-length(edges)])
  list(theta = th_centers, F = Ftheta)
}

fourier_mode_score <- function(theta, Ftheta, m) {
  if (length(theta) == 0) return(0)
  z <- sum(Ftheta * exp(-1i * m * theta))
  denom <- sum(abs(Ftheta)) + 1e-12
  as.numeric(Mod(z) / denom)
}

anisotropy_index <- function(Ftheta) {
  if (length(Ftheta) == 0) return(NA_real_)
  (max(Ftheta) - mean(Ftheta)) / (sum(Ftheta) + 1e-12)
}

field_skewness <- function(u) {
  v <- as.numeric(u); v <- v - mean(v)
  s2 <- mean(v^2); s3 <- mean(v^3)
  if (s2 <= 1e-16) return(0)
  s3 / (s2^(3/2))
}

classify_pattern <- function(S2, H6, A, skew, rules) {
  if (!is.finite(S2)) S2 <- 0
  if (!is.finite(H6)) H6 <- 0
  if (!is.finite(A))  A  <- 1
  if (!is.finite(skew)) skew <- 0
  
  if (S2 > rules$stripe_S2_thresh && H6 < rules$stripe_H6_max) return("stripes")
  if (H6 > rules$hex_H6_min) {
    if (skew >  rules$spots_skew_min)  return("hex_up")
    if (skew < -rules$spots_skew_min)  return("hex_down")
    return("hexagons")
  }
  if (A < rules$iso_A_max && abs(skew) > rules$spots_skew_min)
    return(ifelse(skew > 0, "spots_up", "spots_down"))
  "labyrinth"
}

# ------------------------- SWEEP -------------------------
results <- list()
thumb_records <- list()

for (rr in RG_R) {
  for (ggg in RG_G) {
    u_final <- simulate_SH_fd(rr, ggg, steps = steps, dt = dt, save_png = save_samples)
    
    # Diagnostics (single FFT)
    uh <- fft2(u_final); S <- Mod(uh)^2
    dk <- dominant_k(S)
    kstar <- dk$kstar
    ring <- angular_spectrum(uh, kstar, halfwidth_frac = ring_halfwidth_frac, nbin = n_theta_bins)
    S2 <- fourier_mode_score(ring$theta, ring$F, m = 2)
    H6 <- fourier_mode_score(ring$theta, ring$F, m = 6)
    A  <- anisotropy_index(ring$F)
    skew <- field_skewness(u_final)
    Lscale <- if (is.finite(kstar) && kstar > 0) 2*pi/kstar else NA_real_
    cls <- classify_pattern(S2, H6, A, skew, class_rules)
    
    results[[length(results)+1]] <- data.frame(
      r = rr, g = ggg, kstar = kstar, L = Lscale,
      S2 = S2, H6 = H6, A = A, skew = skew, class = cls
    )
    
    # Thumbnail
    df <- as.data.frame(as.table(u_final)); colnames(df) <- c("y","x","val")
    p <- ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
      scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"))
    thumb_fn <- file.path(thumb_dir, sprintf("thumb_r%.3f_g%.3f.png", rr, ggg))
    ggsave(thumb_fn, p, width = 2.2, height = 2.2, dpi = 120)
    thumb_records[[length(thumb_records)+1]] <- data.frame(r = rr, g = ggg, thumb = thumb_fn, class = cls)
    
    message(sprintf("[DONE] r=%.3f g=%.3f  => class=%s  (k*=%.3f, S2=%.2f, H6=%.2f, A=%.2f, skew=%.2f)",
                    rr, ggg, cls, kstar, S2, H6, A, skew))
  }
}

res_df <- do.call(rbind, results)
write.csv(res_df, csv_out, row.names = FALSE)
message("Saved metrics: ", csv_out)

# ---------------------- PHASE MAP PLOT --------------------
pal <- c(
  stripes   = "#1f78b4",
  hexagons  = "#33a02c",
  hex_up    = "#33a02c",
  hex_down  = "#b2df8a",
  spots_up  = "#e31a1c",
  spots_down= "#fb9a99",
  labyrinth = "#6a3d9a"
)

phase_plot <- ggplot(res_df, aes(x = factor(r), y = factor(g), fill = class)) +
  geom_tile(color = "grey30") +
  scale_fill_manual(values = pal) +
  labs(x = "r", y = "g", fill = "class",
       title = "Swift–Hohenberg Phase Map (2D, FD semi-implicit)") +
  theme_minimal(base_size = 12)

ggsave(phase_png, phase_plot, width = 7, height = 5, dpi = 140)
message("Saved phase map: ", phase_png)

print(res_df)
