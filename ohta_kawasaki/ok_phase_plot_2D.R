# ok_phase_plot_fd.R
# Ohta–Kawasaki phase sweep (no FFT anywhere).
# Periodic FD semi-implicit solver, structure-tensor anisotropy + skewness classifier.
# Outputs: ok_phase_metrics.csv, ok_phase_map.png, thumbs/

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
})
library(ggplot2)

# -------- sweep params --------
ALPHAS <- c(0.03, 0.06, 0.10)     # nonlocal strength
MEANS  <- c(0.00, 0.20, -0.20)    # average composition m̄
nx <- 160; ny <- 160
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny
dt <- 0.08
steps <- 1600
epsilon <- 0.02
Mmob <- 1.0
ic_amp <- 0.08
seed <- 17

# classifier thresholds
ANISO_STRIPE_T <- 0.40
SKEW_SPOT_T    <- 0.25

thumb_dir <- "thumbs_ok"; if (!dir.exists(thumb_dir)) dir.create(thumb_dir, recursive = TRUE)
csv_out <- "ok_phase_metrics.csv"; png_out <- "ok_phase_map.png"

set.seed(seed); stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# periodic ops
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]
D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# CG
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

simulate_OK <- function(alpha, m0, steps, dt) {
  phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
  for (s in 1:steps) {
    rhs_psi <- phi - mean(phi)                              # zero-mean RHS
    psi <- cg_solve(function(U) -lap(U), b = rhs_psi, x0 = NULL, tol = 1e-8, maxit = 250, proj_zero_mean = TRUE)
    mu_exp <- (phi^3 - phi) + alpha * psi
    rhs_phi <- phi + dt * Mmob * ( lap(mu_exp) - lap(phi) ) # explicit pieces
    phi <- cg_solve(function(U) U + dt*Mmob*epsilon^2*bilap(U),
                    b = rhs_phi, x0 = phi, tol = 1e-8, maxit = 250)
  }
  phi
}

# real-space metrics
length_scale <- function(phi) {
  num <- sum(phi^2) * dx * dy
  den <- sum(D1x(phi)^2 + D1y(phi)^2) * dx * dy + 1e-30
  sqrt(num / den)
}
skewness <- function(phi) {
  v <- as.numeric(phi); v <- v - mean(v); s2 <- mean(v^2); s3 <- mean(v^3)
  if (s2 <= 1e-16) 0 else s3 / (s2^(3/2))
}
anisotropy <- function(phi) {
  gx <- D1x(phi); gy <- D1y(phi)
  J11 <- mean(gx*gx); J22 <- mean(gy*gy); J12 <- mean(gx*gy)
  tr <- J11 + J22; disc <- sqrt(max((J11 - J22)^2 + 4*J12^2, 0))
  lam1 <- 0.5*(tr + disc); lam2 <- 0.5*(tr - disc)
  if (tr <= 1e-16) 0 else (lam1 - lam2) / (lam1 + lam2)
}
classify <- function(A, S) {
  if (!is.finite(A)) A <- 0; if (!is.finite(S)) S <- 0
  if (A > ANISO_STRIPE_T) return("stripes")
  if (S >  SKEW_SPOT_T)   return("spots_up")
  if (S < -SKEW_SPOT_T)   return("spots_down")
  "labyrinth"
}

# sweep
res <- list()
thumbs <- list()

for (a in ALPHAS) {
  for (m in MEANS) {
    phi <- simulate_OK(alpha = a, m0 = m, steps = steps, dt = dt)
    A <- anisotropy(phi); S <- skewness(phi); Lc <- length_scale(phi)
    cls <- classify(A, S)
    
    res[[length(res)+1]] <- data.frame(alpha = a, mbar = m, anisotropy = A, skew = S, L = Lc, class = cls)
    
    # thumbnail
    df <- as.data.frame(as.table(phi)); colnames(df) <- c("y","x","val")
    p <- ggplot(df, aes(x, y, fill = val)) +
      geom_raster() + coord_fixed() +
      scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
      theme_void() + theme(legend.position = "none")
    fn <- file.path(thumb_dir, sprintf("thumb_a%.3f_m%.2f.png", a, m))
    ggsave(fn, p, width = 2.2, height = 2.2, dpi = 120)
    thumbs[[length(thumbs)+1]] <- data.frame(alpha = a, mbar = m, thumb = fn)
    message(sprintf("[OK] alpha=%.3f  m̄=%.2f  -> %s  (A=%.2f, skew=%.2f, L=%.3f)", a, m, cls, A, S, Lc))
  }
}

res_df <- do.call(rbind, res)
write.csv(res_df, csv_out, row.names = FALSE)

# phase map
pal <- c(stripes = "#1f78b4", spots_up = "#e31a1c", spots_down = "#fb9a99", labyrinth = "#6a3d9a")
phase_plot <- ggplot(res_df, aes(x = factor(alpha), y = factor(mbar), fill = class)) +
  geom_tile(color = "grey30") +
  scale_fill_manual(values = pal) +
  labs(x = "α (nonlocal)", y = "m̄ (mean composition)", fill = "class",
       title = "Ohta–Kawasaki phase map (FD, no FFT)") +
  theme_minimal(base_size = 12)
ggsave(png_out, phase_plot, width = 7, height = 5, dpi = 140)

print(res_df)
message("Saved: ", csv_out, " and ", png_out)
