# swift_hohenberg_2d.R — 2D Swift–Hohenberg with additive noise eta*sqrt(dt)*N(0,1)
# IMEX Euler (linear implicit), frames + MP4 + structure factor metrics
# Run: Rscript swift_hohenberg_2d.R

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
nx <- 256; ny <- 256
Lx <- 2*pi; Ly <- 2*pi
dt <- 0.2
nsteps <- 2500
save_every <- 5
rparam <- 0.2         # linear growth parameter r
gparam <- 1.0         # cubic nonlinearity coefficient g
eta <- 0.05           # noise amplitude (additive); term: eta*sqrt(dt)*N(0,1)
seed <- 4
dealias <- TRUE       # 2/3-rule low-pass on nonlinear term
out_prefix <- "sh2d"
# -----------------------------------------------------

set.seed(seed)
ny <- as.integer(ny); nx <- as.integer(nx)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

kx_vec <- (2*pi/Lx) * c(0:(nx/2-1), -nx/2:-1)
ky_vec <- (2*pi/Ly) * c(0:(ny/2-1), -ny/2:-1)
kx <- matrix(rep(kx_vec, each = ny), nrow = ny, ncol = nx)
ky <- matrix(rep(ky_vec, times = nx), nrow = ny, ncol = nx)
k2 <- kx^2 + ky^2
Lhat <- rparam - (1 - k2)^2  # Fourier symbol of linear operator

# 2/3 dealias mask
mask <- matrix(1, nrow = ny, ncol = nx)
if (dealias) {
  kx_cut <- max(abs(kx_vec)) * 2/3
  ky_cut <- max(abs(ky_vec)) * 2/3
  mask[abs(kx) > kx_cut | abs(ky) > ky_cut] <- 0
}

fft_rows <- function(u, inv = FALSE) t(apply(u, 1, function(x) fft(x, inverse = inv)))
fft_cols <- function(u, inv = FALSE) apply(u, 2, function(x) fft(x, inverse = inv))
fft2  <- function(u) fft_cols(fft_rows(u, inv = FALSE), inv = FALSE)
ifft2 <- function(U) Re(fft_cols(fft_rows(U, inv = TRUE), inv = TRUE) / (nx * ny))

radial_spectrum <- function(uhat_mod2) {
  dkx <- 2*pi/Lx; dky <- 2*pi/Ly; dk <- min(dkx, dky)
  r <- sqrt(k2); bin <- floor(r/dk + 0.5)
  df <- data.frame(bin = as.vector(bin), S = as.vector(uhat_mod2))
  agg <- aggregate(S ~ bin, data = df, FUN = mean)
  list(k = agg$bin * dk, S = agg$S)
}

# IC: small noise
u <- matrix(0.01 * rnorm(nx*ny), nrow = ny, ncol = nx)

frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(mat, tag) {
  df <- as.data.frame(as.table(mat)); colnames(df) <- c("y","x","val")
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("Swift–Hohenberg t=%s (r=%.2f, g=%.2f, eta=%.3f)", tag, rparam, gparam, eta))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}

encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

ts_t <- c(); ts_kstar <- c(); ts_L <- c()

for (step in 0:nsteps) {
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(u, step*dt)
    uhat <- fft2(u)
    S <- Mod(uhat)^2
    rs <- radial_spectrum(S)
    if (length(rs$k) > 1) {
      kpos <- which(rs$k > 0)
      if (length(kpos) > 0) {
        imax <- kpos[which.max(rs$S[kpos])]
        kstar <- rs$k[imax]
        L <- if (kstar > 0) 2*pi/kstar else NA_real_
        ts_t <- c(ts_t, step*dt); ts_kstar <- c(ts_kstar, kstar); ts_L <- c(ts_L, L)
      }
    }
  }
  
  # Nonlinear + noise in physical space
  Nphys <- -gparam * (u^3) + eta * sqrt(dt) * matrix(rnorm(nx*ny), nrow = ny, ncol = nx)
  
  # IMEX (linear implicit):
  # u_{n+1} - dt*L u_{n+1} = u_n + dt*N(u_n)
  rhs_hat <- fft2(u + dt * Nphys)
  denom <- (1 - dt * Lhat)
  unew_hat <- rhs_hat / denom
  u <- ifft2(unew_hat)
  
  # Optional de-aliasing (low-pass) on state after step
  if (dealias) {
    uhat <- fft2(u); uhat <- uhat * mask; u <- ifft2(uhat)
  }
}

encode_video()

if (length(ts_t) > 10) {
  png(sprintf("%s_metrics.png", out_prefix), width = 900, height = 400)
  par(mfrow = c(1,2), mar = c(4,4,2,1))
  plot(ts_t, ts_L, type = "l", xlab = "t", ylab = "L(t)=2π/k*", main = "Pattern scale vs time")
  plot(ts_t, ts_kstar, type = "l", xlab = "t", ylab = "k*(t)", main = "Dominant wavenumber")
  dev.off()
}
message("Done: SH 2D")
