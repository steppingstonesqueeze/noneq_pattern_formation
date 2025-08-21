# ks_1d_fd.R
# 1D Kuramoto–Sivashinsky (real-space finite-difference, semi-implicit, CG)
# PDE: u_t = -u_xx - u_xxxx - u u_x
# Scheme: implicit only in u_xxxx (SPD); explicit for -u_xx and nonlinear flux form.
# Step: (I + dt * D4) u^{n+1} = u^n + dt [ -D2(u^n) - 0.5 * D1( (u^n)^2 ) ]

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
N <- 512                  # grid (even)
L <- 2*pi                 # domain length
dx <- L / N
dt <- 0.001              # stable step
nsteps <- 15000
save_every <- 100
seed <- 3
out_prefix <- "ks1d_fd"
# -----------------------------------------------------

set.seed(seed)
stopifnot(N %% 2 == 0)

# periodic shift helpers
shiftL1 <- function(v) c(v[-1], v[1])
shiftR1 <- function(v) c(v[length(v)], v[-length(v)])

# centered differences (periodic)
D1 <- function(v) (shiftL1(v) - shiftR1(v)) / (2*dx)                   # first derivative
D2 <- function(v) (shiftL1(v) + shiftR1(v) - 2*v) / dx^2               # Laplacian in 1D
D4 <- function(v) D2(D2(v))                                            # biharmonic

# implicit operator A(v) = (I + dt * D4) v (SPD)
apply_A <- function(v) v + dt * D4(v)

# matrix-free CG for 1D SPD
cg_solve_vec <- function(b, x0 = NULL, tol = 1e-8, maxit = 500) {
  x <- if (is.null(x0)) rep(0, length(b)) else x0
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

# initial condition
u <- 0.1 * rnorm(N)

# output helpers
frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L
xgrid <- seq(0, L, length.out = N + 1)[- (N + 1)]

save_frame <- function(u, ttag) {
  df <- data.frame(x = xgrid, u = u)
  p <- ggplot(df, aes(x, u)) +
    geom_line(linewidth = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = sprintf("KS 1D (FD)  t=%.2f", ttag), x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 7, height = 3.5, dpi = 120)
}
encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

# time loop
for (step in 0:nsteps) {
  # Nonlinear (conservative) and linear explicit pieces at time n
  Nflux <- -0.5 * D1(u * u)         # = -u u_x in conservative form
  rhs <- u + dt * ( - D2(u) + Nflux )
  
  # implicit solve for D4
  u <- cg_solve_vec(rhs, x0 = u, tol = 1e-8, maxit = 400)
  
  # guards
  if (!all(is.finite(u))) u[!is.finite(u)] <- 0
  rng_now <- range(u[is.finite(u)])
  if (is.finite(rng_now[1]) && is.finite(rng_now[2]) && abs(rng_now[2] - rng_now[1]) < 1e-14) {
    u <- u + 1e-6 * rnorm(length(u))
  }
  
  # output
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    v <- mean( (u - mean(u))^2 ); m <- range(u[is.finite(u)])
    message(sprintf("t=%.2f, var=%.3e, min=%.4g, max=%.4g", step*dt, v, m[1], m[2]))
    save_frame(u, ttag = step*dt)
  }
}

encode_video()
message("Done: KS 1D (FD)")
