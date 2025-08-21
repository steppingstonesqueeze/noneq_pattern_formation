# ohta_kawasaki_2d_fd.R
# 2D Ohta–Kawasaki (periodic), real-space FD, semi-implicit CH step.
# Nonlocal term via Poisson: -∇^2 ψ = φ - φ̄, solved by matrix-free CG (zero-mean enforced).
# Viz + MP4, and metrics including full OK free energy.

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# ------------- User params -------------
nx <- 192; ny <- 192
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny
dt <- 0.08
nsteps <- 2000
save_every    <- 10
metrics_every <- 10

epsilon <- 0.02
Mmob    <- 1.0
alpha   <- 0.06          # nonlocal strength
m0      <- 0.0           # mean composition
ic_amp  <- 0.10
# ---------------------------------------

set.seed(11)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# periodic shifts & operators
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# CG (2D)
cg_solve <- function(apply_A, b, x0 = NULL, tol = 1e-8, maxit = 300, proj_zero_mean = FALSE) {
  x <- if (is.null(x0)) matrix(0, nrow = nrow(b), ncol = ncol(b)) else x0
  if (proj_zero_mean) b <- b - mean(b)
  r <- b - apply_A(x); if (proj_zero_mean) r <- r - mean(r)
  p <- r; rr <- sum(r*r); bnorm <- sqrt(sum(b*b)) + 1e-30
  for (it in 1:maxit) {
    Ap <- apply_A(p)
    alpha <- rr / (sum(p*Ap) + 1e-30)
    x <- x + alpha * p
    r <- r - alpha * Ap
    if (proj_zero_mean) r <- r - mean(r)
    rr_new <- sum(r*r)
    if (sqrt(rr_new) <= tol * bnorm) break
    beta <- rr_new / (rr + 1e-30)
    p <- r + beta * p
    rr <- rr_new
  }
  if (proj_zero_mean) x <- x - mean(x)
  x
}

# init
phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)

frames_dir <- "ok_frames"
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(phi, tag) {
  df <- as.data.frame(as.table(phi)); colnames(df) <- c("y","x","val")
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("Ohta–Kawasaki φ  (t=%s, α=%.3f, m̄=%.2f)", tag, alpha, m0))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}
encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  if (length(files)) av::av_encode_video(files, output = "ok_run.mp4", framerate = 24)
}

# metrics
length_scale <- function(phi) {
  num <- sum(phi^2) * dx * dy
  den <- sum(D1x(phi)^2 + D1y(phi)^2) * dx * dy + 1e-30
  sqrt(num / den)
}
energy_OK <- function(phi, psi) {
  g <- 0.25*(phi^2 - 1)^2
  e_local <- 0.5*epsilon^2*(D1x(phi)^2 + D1y(phi)^2) + g
  e_nonlocal <- 0.5*alpha * (phi - mean(phi)) * psi
  sum(e_local) * dx * dy + sum(e_nonlocal) * dx * dy
}

ts <- data.frame(t = numeric(0), mass = numeric(0), E = numeric(0), L = numeric(0))

# time loop
for (step in 0:nsteps) {
  # solve Poisson for ψ: -∇^2 ψ = φ - φ̄ (zero-mean RHS)
  rhs_psi <- phi - mean(phi)
  psi <- cg_solve(function(U) -lap(U), b = rhs_psi, x0 = NULL, tol = 1e-8, maxit = 400, proj_zero_mean = TRUE)
  
  # CH-OK step: (I + dt M eps^2 ∇^4) φ^{n+1} = φ^n + dt M ∇^2 [ (φ^3 - φ) + alpha ψ ] - dt M ∇^2 φ^n
  mu_exp <- (phi^3 - phi) + alpha * psi
  rhs_phi <- phi + dt * Mmob * ( lap(mu_exp) - lap(phi) )
  phi <- cg_solve(function(U) U + dt*Mmob*epsilon^2*bilap(U),
                  b = rhs_phi, x0 = phi, tol = 1e-8, maxit = 250)
  
  # guards
  if (!all(is.finite(phi))) phi[!is.finite(phi)] <- 0
  
  # metrics & viz
  if (step %% metrics_every == 0) {
    Lc <- length_scale(phi)
    E  <- energy_OK(phi, psi)
    ts <- rbind(ts, data.frame(t = step*dt, mass = mean(phi), E = E, L = Lc))
    message(sprintf("t=%.2f  mass=%.4f  E=%.3e  L=%.3f", step*dt, mean(phi), E, Lc))
  }
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(phi, tag = step)
  }
}

write.csv(ts, "ok_timeseries.csv", row.names = FALSE)
encode_video()
message("Done: Ohta–Kawasaki (FD)")
