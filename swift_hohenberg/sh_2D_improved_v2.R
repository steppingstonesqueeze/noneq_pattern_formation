# swift_hohenberg_2d_fd_v2.R
# 2D Swift–Hohenberg (real-space FD, semi-implicit, CG) + exponent metrics
# PDE: u_t = r u - (1 + ∇^2)^2 u - g u^3 + eta ξ
#     = (r - 1) u - 2 ∇^2 u - ∇^4 u - g u^3 + eta ξ
# Step: (I + dt ∇^4) u^{n+1} = u^n + dt[ r u^n - 2 ∇^2 u^n - g (u^n)^3 ] + eta √dt N(0,1)

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
nx <- 192; ny <- 192
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx / nx; dy <- Ly / ny   # assume dx=dy
dt <- 0.10
nsteps <- 2000
save_every    <- 10            # PNG cadence
metrics_every <- 10            # metrics cadence (uses one FFT)
rparam <- 0.25                 # try 0.20–0.35
gparam <- 1.0
eta    <- 0.05
seed <- 4
out_prefix <- "sh2d_fd_v2"
# -----------------------------------------------------

set.seed(seed); stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# periodic shifts
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

lap   <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# Implicit SPD operator A(U) = U + dt * ∇^4 U
apply_A <- function(U) U + dt * bilap(U)

# Matrix-free Conjugate Gradient
cg_solve <- function(b, x0 = NULL, tol = 1e-8, maxit = 200) {
  x <- if (is.null(x0)) matrix(0, nrow = ny, ncol = nx) else x0
  r <- b - apply_A(x); p <- r
  rr <- sum(r*r); bnorm <- sqrt(sum(b*b)) + 1e-30
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

# ----- init (small noise) -----
u <- matrix(0.01 * rnorm(nx*ny), nrow = ny, ncol = nx)
message(sprintf("IC check: min=%.4g, max=%.4g", min(u), max(u)))

# ----- output dirs -----
frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)

# ----- 2D FFT helpers for diagnostics only -----
fft_rows <- function(u, inv = FALSE) t(apply(u, 1, function(x) fft(x, inverse = inv)))
fft_cols <- function(u, inv = FALSE) apply(u, 2, function(x) fft(x, inverse = inv))
fft2  <- function(u) fft_cols(fft_rows(u, inv = FALSE), inv = FALSE)

kx_vec <- (2*pi/Lx) * c(0:(nx/2-1), -nx/2:-1)
ky_vec <- (2*pi/Ly) * c(0:(ny/2-1), -ny/2:-1)
kx <- matrix(rep(kx_vec, each = ny), nrow = ny, ncol = nx)
ky <- matrix(rep(ky_vec, times = nx), nrow = ny, ncol = nx)
k2 <- kx^2 + ky^2

radial_spectrum <- function(S) {
  dk <- min(2*pi/Lx, 2*pi/Ly)
  r  <- sqrt(k2); bin <- floor(r/dk + 0.5)
  df <- data.frame(bin = as.vector(bin), S = as.vector(S))
  df <- df[df$bin >= 0, , drop = FALSE]
  if (!nrow(df)) return(list(k=numeric(0), S=numeric(0)))
  agg <- aggregate(S ~ bin, data = df, FUN = mean)
  list(k = agg$bin * dk, S = agg$S)
}

ring_metrics <- function(U) {
  Uh <- fft2(U); S <- Mod(Uh)^2
  rs <- radial_spectrum(S)
  out <- list(A = sqrt(mean(U^2)), kstar = NA_real_, L = NA_real_, dk = NA_real_, xi = NA_real_)
  if (length(rs$k) > 2) {
    pos <- which(rs$k > 0)
    if (length(pos)) {
      imax <- pos[which.max(rs$S[pos])]
      kstar <- rs$k[imax]; out$kstar <- kstar; out$L <- if (kstar>0) 2*pi/kstar else NA_real_
      # ring width via second moment in ±20% window
      wfrac <- 0.20; klo <- (1-wfrac)*kstar; khi <- (1+wfrac)*kstar
      sel <- which(rs$k >= klo & rs$k <= khi)
      if (length(sel) >= 3) {
        kk <- rs$k[sel]; ww <- rs$S[sel]
        mu <- sum(kk*ww)/sum(ww); sig2 <- sum((kk-mu)^2*ww)/sum(ww)
        dk <- sqrt(max(sig2, 0)); out$dk <- dk; out$xi <- if (dk>0) 1/dk else NA_real_
      }
    }
  }
  out
}

# ----- save frame -----
frame_id <- 0L
save_frame <- function(mat, tag) {
  df <- as.data.frame(as.table(mat)); colnames(df) <- c("y","x","val"); df$val <- as.numeric(df$val)
  rng <- range(df$val[is.finite(df$val)])
  if (diff(rng) < .Machine$double.eps) df$val <- df$val + rnorm(nrow(df), sd = 1e-9)
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"), na.value = "grey80") +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("Swift–Hohenberg (FD) t=%s  r=%.2f g=%.2f η=%.2f", tag, rparam, gparam, eta))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}

encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  if (length(files)) av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

# ----- time series storage -----
ts_t <- c(); ts_A <- c(); ts_k <- c(); ts_L <- c(); ts_dk <- c(); ts_xi <- c()

# ----- time loop -----
for (step in 0:nsteps) {
  # explicit RHS
  rhs <- u + dt * ( rparam * u - 2 * lap(u) - gparam * (u^3) ) +
    eta * sqrt(dt) * matrix(rnorm(nx*ny), nrow = ny, ncol = nx)
  # implicit solve
  u <- cg_solve(rhs, x0 = u, tol = 1e-8, maxit = 200)
  
  if (!all(is.finite(u))) u[!is.finite(u)] <- 0
  
  if (step %% metrics_every == 0) {
    metr <- ring_metrics(u)
    ts_t  <- c(ts_t, step*dt); ts_A  <- c(ts_A, metr$A)
    ts_k  <- c(ts_k, metr$kstar); ts_L  <- c(ts_L, metr$L)
    ts_dk <- c(ts_dk, metr$dk);   ts_xi <- c(ts_xi, metr$xi)
    v <- mean((u - mean(u))^2); m <- range(u[is.finite(u)])
    message(sprintf("t=%.2f, A=%.3e, L=%.3f, xi=%.3f, var=%.3e, min=%.3g, max=%.3g",
                    step*dt, metr$A, metr$L, metr$xi, v, m[1], m[2]))
  }
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(u, tag = step)
  }
}

encode_video()

# ----- post: fit exponent n (L ~ t^n) on late-time window -----
ok <- is.finite(ts_t) & is.finite(ts_L) & ts_t > 0 & is.finite(log(ts_L))
tmin <- quantile(ts_t[ok], 0.5, na.rm = TRUE)   # use latter half for fit
win  <- ok & ts_t >= tmin
n_fit <- NA_real_; c_fit <- NA_real_
if (sum(win) >= 6) {
  fit <- lm(log(ts_L[win]) ~ log(ts_t[win]))
  n_fit <- coef(fit)[2]; c_fit <- coef(fit)[1]
  message(sprintf("Coarsening exponent n ≈ %.3f (L ~ t^n) over t ≥ %.2f", n_fit, tmin))
}

# ----- save CSV -----
df_ts <- data.frame(t = ts_t, A = ts_A, kstar = ts_k, L = ts_L, dk = ts_dk, xi = ts_xi)
write.csv(df_ts, sprintf("%s_timeseries.csv", out_prefix), row.names = FALSE)

# ----- plots -----
png(sprintf("%s_exponents.png", out_prefix), width = 1000, height = 450)
par(mfrow = c(1,2), mar = c(4,4,2,1))
# L(t) vs t
plot(ts_t, ts_L, log = "xy", pch = 16, cex = 0.7,
     xlab = "t", ylab = "L(t) = 2π/k*(t)", main = "Swift–Hohenberg coarsening")
if (is.finite(n_fit)) {
  tt <- ts_t[win]; yy <- exp(c_fit) * tt^n_fit
  lines(tt, yy, col = "red", lwd = 2)
  legend("topleft", bty="n", legend = sprintf("n ≈ %.3f", n_fit), text.col = "red")
}
# xi(t) vs t
plot(ts_t, ts_xi, log = "xy", pch = 16, cex = 0.7,
     xlab = "t", ylab = "ξ(t) ≈ 1/Δk", main = "Ring width length scale")
dev.off()

message("Saved: timeseries CSV + exponent plots")
