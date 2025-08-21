# swift_hohenberg_2d_fd.R
# 2D Swift–Hohenberg (real-space finite-difference, semi-implicit, CG)
# PDE: u_t = r u - (1 + ∇^2)^2 u - g u^3 + eta * ξ
#      = (r - 1) u - 2 ∇^2 u - ∇^4 u - g u^3 + eta * ξ
# Scheme: implicit only in ∇^4 (SPD); explicit for r u, ∇^2 term, and -g u^3, plus noise.
# Step: (I + dt * ∇^4) u^{n+1} = u^n + dt[ r u^n - 2 ∇^2 u^n - g (u^n)^3 ] + eta * sqrt(dt) * N(0,1)

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
nx <- 192; ny <- 192              # grid (even)
Lx <- 2*pi; Ly <- 2*pi            # domain
dx <- Lx / nx; dy <- Ly / ny      # (assume dx=dy)
dt <- 0.01                        # stable step for explicit pieces
nsteps <- 15000
save_every <- 100
rparam <- 0.2                     # r
gparam <- 1.0                     # g
eta <- 0.05                       # noise amplitude (additive)
seed <- 4
out_prefix <- "sh2d_fd"
# -----------------------------------------------------

set.seed(seed)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# periodic shift helpers
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

# 5-point Laplacian and biharmonic (periodic)
lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# Implicit operator A(U) = (I + dt * ∇^4) U  (SPD)
apply_A <- function(U) U + dt * bilap(U)

# Matrix-free Conjugate Gradient for SPD systems
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

# --- initial condition: small noise ---
u <- matrix(0.01 * rnorm(nx*ny), nrow = ny, ncol = nx)
icrng <- range(as.numeric(u))
message(sprintf("IC check: min=%.4g, max=%.4g", icrng[1], icrng[2]))

# --- output helpers ---
frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(mat, tag, title_txt = "Swift–Hohenberg (FD)") {
  if (!any(is.finite(mat))) mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  rng <- range(mat[is.finite(mat)])
  if (diff(rng) < .Machine$double.eps) mat <- mat + rnorm(length(mat), sd = 1e-9)
  df <- as.data.frame(as.table(mat)); colnames(df) <- c("y","x","val"); df$val <- as.numeric(df$val)
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"), na.value = "grey80") +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("%s t=%s (r=%.2f, g=%.2f, eta=%.3f)", title_txt, tag, rparam, gparam, eta))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}
encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

# --- time loop ---
for (step in 0:nsteps) {
  # explicit RHS
  rhs <- u + dt * ( rparam * u - 2 * lap(u) - gparam * (u^3) ) +
    eta * sqrt(dt) * matrix(rnorm(nx*ny), nrow = ny, ncol = nx)
  
  # implicit solve
  u <- cg_solve(rhs, x0 = u, tol = 1e-8, maxit = 200)
  
  # guards
  if (!all(is.finite(u))) u[!is.finite(u)] <- 0
  rng_now <- range(u[is.finite(u)])
  if (is.finite(rng_now[1]) && is.finite(rng_now[2]) && abs(rng_now[2] - rng_now[1]) < 1e-14) {
    u <- u + 1e-6 * matrix(rnorm(length(u)), nrow = ny, ncol = nx)
  }
  
  # output
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    v <- mean((u - mean(u))^2); m <- range(u[is.finite(u)])
    message(sprintf("t=%.2f, var=%.3e, min=%.4g, max=%.4g", step*dt, v, m[1], m[2]))
    save_frame(u, tag = step)
  }
}

encode_video()
message("Done: SH 2D (FD)")
