# ohta_kawasaki_2d_fd.R  (robust V3)
# 2D Ohta–Kawasaki (periodic), real-space FD, semi-implicit CH step.
# Nonlocal: -∇^2 ψ = φ - φ̄ (zero-mean RHS), solved by matrix-free CG.
# Stabilized IMEX: (I - dt M ∇² + dt M ε² ∇⁴) φ^{n+1} = φ^n + dt M [ ∇²(φ^3) + α ∇² ψ ].
# Extras: NA/overflow guards, mass recentering to m0, MP4 + metrics.

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# ---------------- User params ----------------
nx <- 192; ny <- 192
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny        # assume dx=dy
dt <- 0.01                    # smaller, safer
nsteps <- 2000
save_every    <- 10
metrics_every <- 10

epsilon <- 0.02
Mmob    <- 1.0
alpha   <- 0.06                 # nonlocal strength
m0      <- 0.0                  # target mean (conserved)
ic_amp  <- 0.10

# CG tolerances (looser early, tighter late)
tol_initial <- 2e-6
tol_final   <- 1e-8

# Soft cap for overflow protection (post-update clamp)
cap_phi <- 5     # |phi| cap (raise to 8–10 if needed)
cap_psi <- 50    # |psi| cap (psi can be larger; adjust if needed)
# ---------------------------------------------

set.seed(11)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# ---------- periodic shifts & operators ----------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap  <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

# ---------- utils ----------
sanitize <- function(U, cap) {
  U[!is.finite(U)] <- 0
  if (!missing(cap) && is.finite(cap)) {
    U[U >  cap] <-  cap
    U[U < -cap] <- -cap
  }
  U
}

# ---------- NA-safe matrix-free CG ----------
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
    
    alpha <- rr / denom
    if (!is.finite(alpha)) alpha <- 0
    
    x <- make_finite(x + alpha * p)
    r <- make_finite(r - alpha * Ap)
    if (proj_zero_mean) r <- r - mean(r)
    
    rr_new <- safe_dot(r, r)
    if (!is.finite(rr_new)) {
      warning("CG: non-finite residual; returning current iterate")
      return(make_finite(x))
    }
    if (sqrt(rr_new) <= tol * bnorm) return(make_finite(x))
    
    beta <- rr_new / (rr + 1e-30)
    if (!is.finite(beta)) beta <- 0
    
    p  <- make_finite(r + beta * p)
    rr <- rr_new
  }
  make_finite(x)
}

# ---------- init ----------
phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
psi <- matrix(0, ny, nx)  # warm-start Poisson

# ---------- viz ----------
frames_dir <- "ok_frames"
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L

save_frame <- function(phi, tag) {
  df <- as.data.frame(as.table(phi)); colnames(df) <- c("y","x","val")
  rng <- range(df$val[is.finite(df$val)])
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(diff(rng)) < .Machine$double.eps) {
    df$val <- df$val + rnorm(nrow(df), sd = 1e-9)
  }
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"), na.value = "grey80") +
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

# ---------- metrics ----------
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
message(sprintf("IC: mass=%.4f, L=%.3f", mean(phi), length_scale(phi)))

# ---------- time loop ----------
for (step in 0:nsteps) {
  tol_now <- if (step < nsteps/2) tol_initial else tol_final
  
  # Poisson: -∇² ψ = φ - φ̄ (zero-mean RHS)
  rhs_psi <- phi - mean(phi)
  psi <- cg_solve(function(U) -lap(U), b = rhs_psi, x0 = psi, tol = tol_now, maxit = 400, proj_zero_mean = TRUE)
  psi <- sanitize(psi, cap = cap_psi)
  
  # IMEX CH-OK:
  # (I - dt M ∇² + dt M ε² ∇⁴) φ^{n+1} = φ^n + dt M [ ∇²(φ^3) + α ∇² ψ ]
  rhs_phi <- phi + dt * Mmob * ( lap(phi^3) + alpha * lap(psi) )
  A_phi   <- function(U) U - dt*Mmob*lap(U) + dt*Mmob*epsilon^2*bilap(U)
  phi <- cg_solve(A_phi, b = rhs_phi, x0 = phi, tol = tol_now, maxit = 300)
  phi <- sanitize(phi, cap = cap_phi)
  
  # Hard enforce mean (mass conservation to m0)
  phi <- phi + (m0 - mean(phi))
  
  # metrics & viz
  if (step %% metrics_every == 0) {
    Lc <- length_scale(phi)
    E  <- energy_OK(phi, psi)
    ms <- mean(phi)
    ts <- rbind(ts, data.frame(t = step*dt, mass = ms, E = E, L = Lc))
    message(sprintf("t=%.2f  mass=%.4f  E=%s  L=%.3f",
                    step*dt, ms, format(E, digits = 4), Lc))
  }
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(phi, tag = step)
  }
}

write.csv(ts, "ok_timeseries.csv", row.names = FALSE)
encode_video()
message("Done: Ohta–Kawasaki (FD, V3)")
