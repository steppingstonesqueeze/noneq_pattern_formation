# ok_phase_plot_fd.R  (robust V2)
# Ohta–Kawasaki phase sweep (no FFT anywhere).
# Periodic FD semi-implicit solver, stabilized IMEX, NA-safe CG, mass recentering.
# Classifier: structure-tensor anisotropy + skewness (real space).
# Outputs: ok_phase_metrics.csv, ok_phase_map.png, thumbs/

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
})
library(ggplot2)

# -------- sweep params --------
ALPHAS <- c(0.03, 0.06, 0.10)     # nonlocal strength α
MEANS  <- c(0.00, 0.20, -0.20)    # average composition m̄
nx <- 160; ny <- 160
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny          # assume dx = dy
dt <- 0.05                        # smaller, safer default
steps <- 1600
epsilon <- 0.02
Mmob <- 1.0
ic_amp <- 0.08
seed <- 17

# CG tolerances (looser early, tighter late)
tol_initial <- 2e-6
tol_final   <- 1e-8

# Soft caps (safety nets; raise if too tight)
cap_phi <- 5     # |phi| cap
cap_psi <- 50    # |psi| cap

# classifier thresholds
ANISO_STRIPE_T <- 0.40
SKEW_SPOT_T    <- 0.25

thumb_dir <- "thumbs_ok"; if (!dir.exists(thumb_dir)) dir.create(thumb_dir, recursive = TRUE)
csv_out <- "ok_phase_metrics.csv"; png_out <- "ok_phase_map.png"

set.seed(seed); stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# ---------- periodic ops ----------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]
D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# ---------- helpers ----------
sanitize <- function(U, cap = NA_real_) {
  U[!is.finite(U)] <- 0
  if (is.finite(cap)) {
    U[U >  cap] <-  cap
    U[U < -cap] <- -cap
  }
  U
}

# ---------- NA-safe CG (matrix-free) ----------
cg_solve <- function(apply_A, b, x0 = NULL, tol = 1e-8, maxit = 300, proj_zero_mean = FALSE) {
  safe_dot <- function(a, b) {
    s <- sum(as.numeric(a) * as.numeric(b), na.rm = TRUE)
    if (!is.finite(s)) 0 else s
  }
  make_finite <- function(M) { M[!is.finite(M)] <- 0; M }
  
  x <- if (is.null(x0)) matrix(0, nrow = nrow(b), ncol = ncol(b)) else x0
  if (proj_zero_mean) b <- b - mean(b)
  b <- make_finite(b)
  
  Ax <- make_finite(apply_A(x))
  r  <- make_finite(b - Ax)
  if (proj_zero_mean) r <- r - mean(r)
  p  <- r
  
  rr <- safe_dot(r, r)
  bnorm <- sqrt(max(safe_dot(b, b), 0)) + 1e-30
  
  for (it in 1:maxit) {
    Ap <- make_finite(apply_A(p))
    denom <- safe_dot(p, Ap) + 1e-30
    if (!is.finite(denom) || denom <= 0) {
      warning("CG: non-positive or non-finite denominator; returning current iterate")
      return(make_finite(x))
    }
    alpha <- rr / denom; if (!is.finite(alpha)) alpha <- 0
    
    x <- make_finite(x + alpha * p)
    r <- make_finite(r - alpha * Ap)
    if (proj_zero_mean) r <- r - mean(r)
    
    rr_new <- safe_dot(r, r)
    if (!is.finite(rr_new)) {
      warning("CG: non-finite residual; returning current iterate")
      return(make_finite(x))
    }
    if (sqrt(rr_new) <= tol * bnorm) return(make_finite(x))
    
    beta <- rr_new / (rr + 1e-30); if (!is.finite(beta)) beta <- 0
    p  <- make_finite(r + beta * p)
    rr <- rr_new
  }
  make_finite(x)
}

# ---------- one OK simulation (stabilized IMEX) ----------
simulate_OK <- function(alpha, m0, steps, dt) {
  phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
  psi <- matrix(0, nrow = ny, ncol = nx)  # warm-start Poisson
  for (s in 1:steps) {
    tol_now <- if (s < steps/2) tol_initial else tol_final
    
    # Poisson: -∇² ψ = φ - φ̄ (zero-mean RHS)
    rhs_psi <- phi - mean(phi)
    psi <- cg_solve(function(U) -lap(U), b = rhs_psi, x0 = psi, tol = tol_now, maxit = 350, proj_zero_mean = TRUE)
    psi <- sanitize(psi, cap = cap_psi)
    
    # Stabilized IMEX CH–OK:
    # (I - dt M ∇² + dt M ε² ∇⁴) φ^{n+1} = φ^n + dt M [ ∇²(φ^3) + α ∇² ψ ]
    rhs_phi <- phi + dt * Mmob * ( lap(phi^3) + alpha * lap(psi) )
    A_phi   <- function(U) U - dt*Mmob*lap(U) + dt*Mmob*epsilon^2*bilap(U)
    phi <- cg_solve(A_phi, b = rhs_phi, x0 = phi, tol = tol_now, maxit = 300)
    phi <- sanitize(phi, cap = cap_phi)
    
    # Mass recentering (conserved mean)
    phi <- phi + (m0 - mean(phi))
  }
  phi
}

# ---------- real-space metrics ----------
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

# ---------- sweep ----------
res <- list()
thumbs <- list()

for (a in ALPHAS) {
  for (m in MEANS) {
    phi <- simulate_OK(alpha = a, m0 = m, steps = steps, dt = dt)
    phi <- sanitize(phi, cap = cap_phi)
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
    message(sprintf("[OK] α=%.3f  m̄=%.2f  -> %s  (A=%.2f, skew=%.2f, L=%.3f)", a, m, cls, A, S, Lc))
  }
}

res_df <- do.call(rbind, res)
write.csv(res_df, csv_out, row.names = FALSE)

# ---------- phase map ----------
pal <- c(stripes = "#1f78b4", spots_up = "#e31a1c", spots_down = "#fb9a99", labyrinth = "#6a3d9a")
phase_plot <- ggplot(res_df, aes(x = factor(alpha), y = factor(mbar), fill = class)) +
  geom_tile(color = "grey30") +
  scale_fill_manual(values = pal) +
  labs(x = "α (nonlocal)", y = "m̄ (mean composition)", fill = "class",
       title = "Ohta–Kawasaki phase map (FD, stabilized IMEX, no FFT)") +
  theme_minimal(base_size = 12)
ggsave(png_out, phase_plot, width = 7, height = 5, dpi = 140)

print(res_df)
message("Saved: ", csv_out, " and ", png_out)
