# model_h_2d_fd.R
# 2D Model H: CH + incompressible Navier–Stokes (periodic), real-space FD, semi-implicit.
# Elliptic solves (CH biharmonic, Helmholtz for viscous step, Poisson for pressure) via matrix-free CG.
# Viz: phi raster + downsampled quiver of velocity; MP4 via av.
# Metrics: mass, CH free energy, kinetic energy, length scale L ~ sqrt(∫phi^2 / ∫|∇phi|^2).

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
  ok <- requireNamespace("av", quietly = TRUE); if (!ok) stop("Install av")
})
library(ggplot2)

# ---------------- User params ----------------
nx <- 160; ny <- 160              # grid (even)
Lx <- 2*pi; Ly <- 2*pi
dx <- Lx/nx; dy <- Ly/ny          # assume dx=dy
dt <- 0.05                        # small for advective CFL
nsteps <- 2000
save_every    <- 10               # frames cadence
metrics_every <- 10

# CH params
epsilon <- 0.02
Mmob    <- 1.0
m0      <- 0.0                    # mean composition
ic_amp  <- 0.10                   # random IC amplitude

# Fluid params
rho <- 1.0
nu  <- 0.5                        # kinematic viscosity (damps fast)
# ---------------------------------------------

set.seed(7)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# ---------- periodic shifts ----------
shiftL <- function(U) U[, c(2:ncol(U), 1)]
shiftR <- function(U) U[, c(ncol(U), 1:(ncol(U)-1))]
shiftU <- function(U) U[c(2:nrow(U), 1), ]
shiftD <- function(U) U[c(nrow(U), 1:(nrow(U)-1)), ]

# ---------- derivatives ----------
D1x <- function(U) (shiftL(U) - shiftR(U)) / (2*dx)
D1y <- function(U) (shiftU(U) - shiftD(U)) / (2*dy)
lap <- function(U) (shiftL(U) + shiftR(U) + shiftU(U) + shiftD(U) - 4*U) / dx^2
bilap <- function(U) lap(lap(U))

div2D <- function(ux, uy) D1x(ux) + D1y(uy)
gradp <- function(p) list(px = D1x(p), py = D1y(p))

# ---------- CG solvers ----------
cg_solve <- function(apply_A, b, x0 = NULL, tol = 1e-8, maxit = 300, proj_zero_mean = FALSE) {
  x <- if (is.null(x0)) matrix(0, nrow = nrow(b), ncol = ncol(b)) else x0
  if (proj_zero_mean) { b <- b - mean(b) }  # enforce solvability for periodic Poisson
  Ax <- apply_A(x); r <- b - Ax; if (proj_zero_mean) r <- r - mean(r)
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

# Operators:
A_CH   <- function(U) U + dt * Mmob * epsilon^2 * bilap(U)        # (I + dt M eps^2 ∇^4)
A_Helm <- function(U) U - dt * (nu) * lap(U)                      # (I - dt nu ∇^2)
A_Pois <- function(U) -lap(U)                                     # -∇^2 (SPD on zero-mean)

# ---------- init fields ----------
phi <- matrix(m0 + ic_amp*(2*matrix(runif(nx*ny), ny, nx) - 1), nrow = ny, ncol = nx)
ux  <- matrix(0, ny, nx)
uy  <- matrix(0, ny, nx)
p   <- matrix(0, ny, nx)

# ---------- viz helpers ----------
frames_dir <- "modelh_frames"
if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
frame_id <- 0L
quiv_stride <- 6  # downsample quiver grid

save_frame <- function(phi, ux, uy, ttag) {
  # scalar field
  df <- as.data.frame(as.table(phi)); colnames(df) <- c("y","x","val")
  p1 <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("Model H: composition φ  (t=%.2f)", ttag))
  
  # quiver: speed background + arrows
  speed <- sqrt(ux^2 + uy^2)
  df2 <- as.data.frame(as.table(speed)); colnames(df2) <- c("y","x","spd")
  # arrow grid
  xs <- seq(1, ncol(ux), by = quiv_stride); ys <- seq(1, nrow(ux), by = quiv_stride)
  qx <- ux[ys, xs]; qy <- uy[ys, xs]
  qdf <- expand.grid(x = xs, y = ys); qdf$ux <- as.numeric(qx); qdf$uy <- as.numeric(qy)
  scl <- 0.6 * quiv_stride  # visual scale
  p2 <- ggplot(df2, aes(x, y, fill = spd)) +
    geom_raster() + coord_fixed() +
    scale_fill_gradientn(colors = c("#f7fbff","#c6dbef","#6baed6","#2171b5")) +
    geom_segment(data = qdf,
                 aes(x = x, y = y, xend = x + scl*ux, yend = y + scl*uy),
                 inherit.aes = FALSE, linewidth = 0.4, color = "black",
                 arrow = arrow(length = unit(0.08, "cm"))) +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("velocity |u| and arrows (t=%.2f)", ttag))
  
  fn1 <- file.path(frames_dir, sprintf("frame_phi_%05d.png", frame_id))
  fn2 <- file.path(frames_dir, sprintf("frame_u_%05d.png", frame_id))
  ggsave(fn1, p1, width = 5, height = 5, dpi = 120)
  ggsave(fn2, p2, width = 5, height = 5, dpi = 120)
}

encode_video <- function() {
  files1 <- list.files(frames_dir, pattern = "frame_phi_\\d+\\.png$", full.names = TRUE)
  files2 <- list.files(frames_dir, pattern = "frame_u_\\d+\\.png$",   full.names = TRUE)
  files1 <- files1[order(files1)]; files2 <- files2[order(files2)]
  if (length(files1)) av::av_encode_video(files1, output = "modelh_phi.mp4", framerate = 24)
  if (length(files2)) av::av_encode_video(files2, output = "modelh_u.mp4",   framerate = 24)
}

# ---------- metrics ----------
free_energy_CH <- function(phi) {
  g <- 0.25*(phi^2 - 1)^2
  e <- 0.5*epsilon^2 * (D1x(phi)^2 + D1y(phi)^2) + g
  sum(e) * dx * dy
}
length_scale <- function(phi) {
  num <- sum(phi^2) * dx * dy
  den <- sum(D1x(phi)^2 + D1y(phi)^2) * dx * dy + 1e-30
  sqrt(num / den)
}

ts <- data.frame(t = numeric(0), mass = numeric(0), F = numeric(0),
                 L = numeric(0), KE = numeric(0))

# ---------- time loop ----------
for (step in 0:nsteps) {
  
  # ---- CH step with advection u·∇phi (explicit), ∇^4 implicit ----
  mu_exp <- (phi^3 - phi)           # f'(phi)
  rhs_phi <- phi + dt * ( - (ux*D1x(phi) + uy*D1y(phi)) + Mmob * ( lap(mu_exp) - lap(phi) ) )
  phi <- cg_solve(function(U) U + dt*Mmob*epsilon^2*bilap(U),
                  b = rhs_phi, x0 = phi, tol = 1e-8, maxit = 200)
  
  # chemical potential for force
  mu <- -epsilon^2 * lap(phi) + (phi^3 - phi)
  
  # ---- momentum: semi-implicit viscous + projection ----
  # explicit advection + capillary force
  advx <- ux*D1x(ux) + uy*D1y(ux)
  advy <- ux*D1x(uy) + uy*D1y(uy)
  fx <- -(phi)*D1x(mu) / rho
  fy <- -(phi)*D1y(mu) / rho
  
  rhs_ux <- ux + dt * ( -advx + fx )
  rhs_uy <- uy + dt * ( -advy + fy )
  
  # Helmholtz solve (I - dt nu ∇^2) u* = rhs
  ux_star <- cg_solve(function(U) U - dt*nu*lap(U), b = rhs_ux, x0 = ux, tol = 1e-8, maxit = 200)
  uy_star <- cg_solve(function(U) U - dt*nu*lap(U), b = rhs_uy, x0 = uy, tol = 1e-8, maxit = 200)
  
  # Projection: solve ∇^2 p = (rho/dt) ∇·u*
  div_star <- div2D(ux_star, uy_star)
  # Enforce zero-mean RHS (periodic solvability)
  div_star <- div_star - mean(div_star)
  p <- cg_solve(function(U) -lap(U), b = (rho/dt)*div_star, x0 = p, tol = 1e-8, maxit = 300, proj_zero_mean = TRUE)
  
  # Update velocity
  gp <- gradp(p)
  ux <- ux_star - (dt/rho) * gp$px
  uy <- uy_star - (dt/rho) * gp$py
  
  # ---- output & metrics ----
  if (step %% metrics_every == 0) {
    mass <- mean(phi)
    F    <- free_energy_CH(phi)
    Lc   <- length_scale(phi)
    KE   <- 0.5 * sum(ux^2 + uy^2) * dx * dy
    ts <- rbind(ts, data.frame(t = step*dt, mass = mass, F = F, L = Lc, KE = KE))
    message(sprintf("t=%.2f  mass=%.4f  F=%.3e  L=%.3f  KE=%.3e",
                    step*dt, mass, F, Lc, KE))
  }
  if (step %% save_every == 0) {
    frame_id <- frame_id + 1L
    save_frame(phi, ux, uy, ttag = step*dt)
  }
}

write.csv(ts, "modelh_timeseries.csv", row.names = FALSE)
encode_video()
message("Done: Model H (FD)")
