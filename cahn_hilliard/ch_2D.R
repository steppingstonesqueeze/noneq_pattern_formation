# ch_2d.R — 2D Cahn–Hilliard (finite-difference semi-implicit, CG solver)
# PDE: phi_t = M * [ -eps^2 ∇^4 phi  - ∇^2 phi  + ∇^2(phi^3) ]
# Scheme (semi-implicit, constant coefficients):
#   (I + dt*M*eps^2 ∇^4) phi^{n+1} = phi^n + dt*M*( ∇^2(phi^3)^n - ∇^2 phi^n )
# This is SPD -> solved by matrix-free Conjugate Gradient.
#
# Run: Rscript ch_2d.R

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
nx <- 192; ny <- 192         # grid (even). Start modest for speed.
Lx <- 2*pi; Ly <- 2*pi       # domain size
dx <- Lx / nx; dy <- Ly / ny
dt <- 0.001                   # stable & growing
nsteps <- 15000
save_every <- 100
epsilon <- 0.02              # interface width
Mmob <- 1.0
seed <- 1
out_prefix <- "ch_fd"
m0 <- 0.0                    # mean composition
ic_amp <- 0.10               # initial noise
# -----------------------------------------------------

set.seed(seed)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# -------- periodic shift helpers (wrap indices) --------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

# 5-point Laplacian (periodic)
lap <- function(U) {
  (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2  # assume dx=dy
}
bilap <- function(U) lap(lap(U))

# Matrix-free application of implicit operator: A(U) = U + dt*M*eps^2 * ∇^4 U
apply_A <- function(U) U + dt * Mmob * epsilon^2 * bilap(U)

# Conjugate Gradient solver (matrix-free SPD)
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

# ---------- init field ----------
phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
icrng <- range(as.numeric(phi))
message(sprintf("IC check: min=%.4g, max=%.4g", icrng[1], icrng[2]))

# ---------- output utils ----------
frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(mat, tag, title_txt = "Cahn–Hilliard (FD)") {
  # NA/constant-safe plot
  if (!any(is.finite(mat))) mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  rng <- range(mat[is.finite(mat)]); if (diff(rng) < .Machine$double.eps) mat <- mat + rnorm(length(mat), sd = 1e-9)
  df <- as.data.frame(as.table(mat)); colnames(df) <- c("y","x","val"); df$val <- as.numeric(df$val)
  p <- ggplot(df, aes(x, y, fill = val)) + geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"), na.value = "grey80") +
    theme_void() + theme(legend.position = "none") + ggtitle(sprintf("%s t=%s", title_txt, tag))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}
encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

# ---------- time loop ----------
for (step in 0:nsteps) {
  # RHS: phi^n + dt*M*( ∇^2(phi^3)^n - ∇^2 phi^n )
  rhs <- phi + dt * Mmob * ( lap(phi^3) - lap(phi) )
  
  # Implicit solve: (I + dt*M*eps^2 ∇^4) phi^{n+1} = rhs
  phi <- cg_solve(rhs, x0 = phi, tol = 1e-8, maxit = 200)
  
  # Guards
  if (!all(is.finite(phi))) phi[!is.finite(phi)] <- 0
  rng_now <- range(phi[is.finite(phi)])
  if (is.finite(rng_now[1]) && is.finite(rng_now[2]) && abs(rng_now[2] - rng_now[1]) < 1e-14) {
    phi <- phi + 1e-6 * matrix(rnorm(length(phi)), nrow = ny, ncol = nx)
  }
  
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    v <- mean((phi - mean(phi))^2); m <- range(phi[is.finite(phi)])
    message(sprintf("t=%.2f, var=%.3e, min=%.4g, max=%.4g", step*dt, v, m[1], m[2]))
    save_frame(phi, tag = step)
  }
}

encode_video()
message("Done: CH (FD)")
