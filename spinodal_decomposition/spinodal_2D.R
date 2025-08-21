# spinodal_2d.R — Spinodal decomposition via CH (finite-difference semi-implicit, CG)
# Same scheme & solver as ch_2d.R, with quench IC near phi ~ 0.
#
# Run: Rscript spinodal_2d.R

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# -------------------- User params --------------------
nx <- 192; ny <- 192
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx / nx; dy <- Ly / ny
dt <- 0.001
nsteps <- 15000
save_every <- 100
epsilon <- 0.02
Mmob <- 1.0
seed <- 2
out_prefix <- "spinodal_fd"
m0 <- 0.0
noise_amp <- 0.05
# -----------------------------------------------------

set.seed(seed)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# periodic shifts
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

apply_A <- function(U) U + dt * Mmob * epsilon^2 * bilap(U)

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

# init (quench near 0 with small noise)
phi <- matrix(m0 + noise_amp * rnorm(nx*ny), nrow = ny, ncol = nx)
icrng <- range(as.numeric(phi)); message(sprintf("IC check: min=%.4g, max=%.4g", icrng[1], icrng[2]))

frames_dir <- sprintf("%s_frames", out_prefix)
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(mat, tag, title_txt = "Spinodal (FD)") {
  if (!any(is.finite(mat))) mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  rng <- range(mat[is.finite(mat)]); if (diff(rng) < .Machine$double.eps) mat <- mat + rnorm(length(mat), sd = 1e-9)
  df <- as.data.frame(as.table(mat)); colnames(df) <- c("y","x","val"); df$val <- as.numeric(df$val)
  p <- ggplot(df, aes(x, y, fill = val)) + geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#ffffbf","#f46d43","#a50026"), na.value = "grey80") +
    theme_void() + theme(legend.position = "none") + ggtitle(sprintf("%s t=%s", title_txt, tag))
  fn <- file.path(frames_dir, sprintf("frame_%05d.png", frame_id))
  ggsave(fn, p, width = 6, height = 6, dpi = 120)
}
encode_video <- function() {
  files <- list.files(frames_dir, pattern = "frame_\\d+\\.png$", full.names = TRUE)
  files <- files[order(files)]
  av::av_encode_video(files, output = sprintf("%s_run.mp4", out_prefix), framerate = 24)
}

for (step in 0:nsteps) {
  rhs <- phi + dt * Mmob * ( lap(phi^3) - lap(phi) )
  phi <- cg_solve(rhs, x0 = phi, tol = 1e-8, maxit = 200)
  
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
message("Done: Spinodal (FD)")
