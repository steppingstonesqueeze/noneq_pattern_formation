# ks_1d.R — 1D KS with ETDRK4 (Kassam–Trefethen), frames + MP4 + spectrum metrics
# Run: Rscript ks_1d.R

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
N <- 512                 # grid (even)
L <- 2*pi                # domain length
dt <- 0.25               # time step
nsteps <- 4000           # total steps
save_every <- 5
seed <- 3
out_prefix <- "ks1d"
# -----------------------------------------------------

set.seed(seed)
stopifnot(N %% 2 == 0)

x <- seq(0, L, length.out = N + 1)[- (N + 1)]
k <- (2*pi/L) * c(0:(N/2), -N/2+1:-1)  # length N
ik <- 1i * k

# Linear operator Lk: -u_xx - u_xxxx => in Fourier: (+k^2 - k^4)
Lk <- (k^2 - k^4)

# ETDRK4 coefficients
Mquad <- 32L
E <- exp(dt * Lk)
E2 <- exp(dt * Lk / 2)
r <- exp(1i * pi * ((1:Mquad) - 0.5) / Mquad)
LR <- outer(Lk, r, function(a,b) dt * a + b)  # N x M
Q  <- dt * rowMeans((exp(LR/2) - 1) / LR)
f1 <- dt * rowMeans((-4 - LR + exp(LR) * (4 - 3*LR + LR^2)) / (LR^3))
f2 <- dt * rowMeans(( 2 + LR + exp(LR) * (-2 + LR)) / (LR^3))
f3 <- dt * rowMeans((-4 - 3*LR - LR^2 + exp(LR) * (4 - LR)) / (LR^3))

fft1 <- function(u) fft(u)
ifft1 <- function(U) Re(fft(U, inverse = TRUE) / N)

# Nonlinear term N(u) = -0.5 * d_x (u^2)
Nterm_hat <- function(u) {
  u2 <- u * u
  fft_u2 <- fft1(u2)
  -0.5 * ik * fft_u2
}

# IC: small random
u <- 0.1 * rnorm(N)

frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(u, ttag) {
  df <- data.frame(x = x, u = u)
  p <- ggplot(df, aes(x, u)) + geom_line(linewidth = 0.8) +
    theme_minimal(base_size = 12) + labs(title = sprintf("KS 1D  t=%.2f", ttag), x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 7, height = 3.5, dpi = 120)
}

encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

ts_t <- c(); ts_kstar <- c(); ts_lambda <- c()

# Time loop (ETDRK4)
u_hat <- fft1(u)
for (step in 0:nsteps) {
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(ifft1(u_hat), step*dt)
    # Spectrum metrics
    U <- u_hat
    P <- Mod(U)^2
    idx <- which(abs(k) > 0)
    if (length(idx) > 0) {
      kstar <- abs(k[idx])[which.max(P[idx])]
      lambda <- if (kstar > 0) 2*pi/kstar else NA_real_
      ts_t <- c(ts_t, step*dt); ts_kstar <- c(ts_kstar, kstar); ts_lambda <- c(ts_lambda, lambda)
    }
  }
  Nv  <- Nterm_hat(ifft1(u_hat))
  a   <- E2 * u_hat + Q * Nv
  Na  <- Nterm_hat(ifft1(a))
  b   <- E2 * u_hat + Q * Na
  Nb  <- Nterm_hat(ifft1(b))
  c   <- E2 * a + Q * (2*Nb - Nv)
  Nc  <- Nterm_hat(ifft1(c))
  u_hat <- E * u_hat + f1 * Nv + 2 * f2 * (Na + Nb) + f3 * Nc
}

encode_video()

if (length(ts_t) > 10) {
  png(sprintf("%s_metrics.png", out_prefix), width = 900, height = 400)
  par(mfrow = c(1,2), mar = c(4,4,2,1))
  plot(ts_t, ts_kstar, type = "l", xlab = "t", ylab = "k*(t)", main = "Dominant wavenumber")
  plot(ts_t, ts_lambda, type = "l", xlab = "t", ylab = "λ(t)=2π/k*", main = "Dominant wavelength")
  dev.off()
}
message("Done: KS 1D")
