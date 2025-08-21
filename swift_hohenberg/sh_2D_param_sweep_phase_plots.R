# sh_phase_plot.R — Phase map for Swift–Hohenberg (2D)
# Sweeps (r,g), simulates to late time, computes diagnostics, classifies pattern,
# and saves a phase map + CSV + thumbnails.
#
# Run: Rscript sh_phase_plot.R
# Requires: ggplot2

suppressPackageStartupMessages({
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install ggplot2")
})

library(ggplot2)

# ----------------------- USER CONFIG -----------------------
# Parameter grid
RG_R <- c(0.05, 0.15, 0.25)     # linear growth parameter r
RG_G <- c(0.5, 1.0, 1.5)        # cubic nonlinearity g

# Simulation params
nx <- 160; ny <- 160            # keep modest for sweeping; must be even
Lx <- 2*pi; Ly <- 2*pi
dt <- 0.2
steps <- 1600                   # total time steps per (r,g)
save_samples <- FALSE           # set TRUE to also drop evolution frames per (r,g)
sample_every <- 100
eta <- 0.05                     # additive noise amplitude (eta*sqrt(dt)*N(0,1))
seed <- 42
dealias <- TRUE                 # 2/3-rule low-pass on state each step

# Ring detection / angular spectrum
ring_halfwidth_frac <- 0.08     # width around k* for angle spectrum (fraction of k*)
n_theta_bins <- 90              # number of angular bins to analyze
class_rules <- list(
  stripe_S2_thresh = 0.40,
  stripe_H6_max    = 0.25,
  hex_H6_min       = 0.30,
  iso_A_max        = 0.15,
  spots_skew_min   = 0.20       # absolute skewness threshold
)

# Output
csv_out <- "sh_phase_metrics.csv"
phase_png <- "sh_phase_map.png"
thumb_dir <- "thumbs"
if (!dir.exists(thumb_dir)) dir.create(thumb_dir, recursive = TRUE)

set.seed(seed)
stopifnot(nx %% 2 == 0, ny %% 2 == 0)

# -------------------- FOURIER HELPERS ---------------------
kx_vec <- (2*pi/Lx) * c(0:(nx/2-1), -nx/2:-1)
ky_vec <- (2*pi/Ly) * c(0:(ny/2-1), -ny/2:-1)
kx <- matrix(rep(kx_vec, each = ny), nrow = ny, ncol = nx)
ky <- matrix(rep(ky_vec, times = nx), nrow = ny, ncol = nx)
k2 <- kx^2 + ky^2

fft_rows <- function(u, inv = FALSE) t(apply(u, 1, function(x) fft(x, inverse = inv)))
fft_cols <- function(u, inv = FALSE) apply(u, 2, function(x) fft(x, inverse = inv))
fft2  <- function(u) fft_cols(fft_rows(u, inv = FALSE), inv = FALSE)
ifft2 <- function(U) Re(fft_cols(fft_rows(U, inv = TRUE), inv = TRUE) / (nx * ny))

# Dealias mask (2/3 rule)
mask <- matrix(1, nrow = ny, ncol = nx)
if (dealias) {
  kx_cut <- max(abs(kx_vec)) * 2/3
  ky_cut <- max(abs(ky_vec)) * 2/3
  mask[abs(kx) > kx_cut | abs(ky) > ky_cut] <- 0
}

# ------------------ DIAGNOSTIC HELPERS --------------------
radial_avg <- function(S) {
  # returns k values and radially averaged S(k)
  dkx <- 2*pi/Lx; dky <- 2*pi/Ly; dk <- min(dkx, dky)
  r <- sqrt(k2); bin <- floor(r/dk + 0.5)
  df <- data.frame(bin = as.vector(bin), S = as.vector(S))
  agg <- aggregate(S ~ bin, data = df, FUN = mean)
  k <- agg$bin * dk
  list(k = k, S = agg$S)
}

dominant_k <- function(S) {
  rs <- radial_avg(S)
  if (length(rs$k) < 3) return(list(kstar = NA_real_, idx = NA_integer_, k = rs$k, Sr = rs$S))
  pos <- which(rs$k > 0)
  if (length(pos) == 0) return(list(kstar = NA_real_, idx = NA_integer_, k = rs$k, Sr = rs$S))
  i <- pos[which.max(rs$S[pos])]
  list(kstar = rs$k[i], idx = i, k = rs$k, Sr = rs$S)
}

angular_spectrum <- function(Uhat, kstar, halfwidth_frac = 0.08, nbin = 90) {
  if (!is.finite(kstar) || kstar <= 0) return(list(theta = numeric(0), F = numeric(0)))
  r <- sqrt(k2)
  ring_lo <- (1 - halfwidth_frac) * kstar
  ring_hi <- (1 + halfwidth_frac) * kstar
  sel <- (r >= ring_lo) & (r <= ring_hi)
  if (!any(sel)) return(list(theta = numeric(0), F = numeric(0)))
  
  # power on ring
  Uh <- Uhat[sel]; kx_s <- kx[sel]; ky_s <- ky[sel]
  pow <- Mod(Uh)^2
  theta <- atan2(ky_s, kx_s)  # (-pi, pi]
  # bin theta
  edges <- seq(-pi, pi, length.out = nbin + 1)
  bins <- cut(theta, breaks = edges, include.lowest = TRUE, labels = FALSE)
  Ftheta <- numeric(nbin)
  for (i in seq_len(nbin)) {
    Ftheta[i] <- sum(pow[bins == i])
  }
  th_centers <- 0.5 * (edges[-1] + edges[-length(edges)])
  list(theta = th_centers, F = Ftheta)
}

fourier_mode_score <- function(theta, Ftheta, m) {
  # normalized |sum F(theta) * e^{-im theta}|
  if (length(theta) == 0) return(0)
  z <- sum(Ftheta * exp(-1i * m * theta))
  denom <- sum(abs(Ftheta)) + 1e-12
  as.numeric(Mod(z) / denom)
}

anisotropy_index <- function(Ftheta) {
  if (length(Ftheta) == 0) return(NA_real_)
  (max(Ftheta) - mean(Ftheta)) / (sum(Ftheta) + 1e-12)
}

field_skewness <- function(u) {
  v <- as.numeric(u)
  v <- v - mean(v)
  s2 <- mean(v^2); s3 <- mean(v^3)
  if (s2 <= 1e-16) return(0)
  s3 / (s2^(3/2))
}

classify_pattern <- function(S2, H6, A, skew, rules) {
  if (!is.finite(S2)) S2 <- 0
  if (!is.finite(H6)) H6 <- 0
  if (!is.finite(A))  A  <- 1
  if (!is.finite(skew)) skew <- 0
  
  if (S2 > rules$stripe_S2_thresh && H6 < rules$stripe_H6_max) {
    return("stripes")
  }
  if (H6 > rules$hex_H6_min) {
    if (skew > rules$spots_skew_min)  return("hex_up")
    if (skew < -rules$spots_skew_min) return("hex_down")
    return("hexagons")
  }
  if (A < rules$iso_A_max && abs(skew) > rules$spots_skew_min) {
    return(ifelse(skew > 0, "spots_up", "spots_down"))
  }
  return("labyrinth")
}

# ------------------ SH TIME-STEPPER (IMEX) ----------------
step_SH <- function(u, rparam, gparam, dt) {
  # Swift–Hohenberg: u_t = r u - (1 + ∇^2)^2 u - g u^3 + eta*sqrt(dt)*ξ
  Lhat <- rparam - (1 - k2)^2
  Nphys <- -gparam * (u^3) + eta * sqrt(dt) * matrix(rnorm(nx*ny), nrow = ny, ncol = nx)
  rhs_hat <- fft2(u + dt * Nphys)
  unew_hat <- rhs_hat / (1 - dt * Lhat)
  u <- ifft2(unew_hat)
  if (dealias) {
    uh <- fft2(u); uh <- uh * mask; u <- ifft2(uh)
  }
  u
}

simulate_SH <- function(rparam, gparam, steps, dt, u0 = NULL, save_png = FALSE, tag = "") {
  if (is.null(u0)) {
    u <- matrix(0.01 * rnorm(nx*ny), nrow = ny, ncol = nx)
  } else {
    u <- u0
  }
  frames_dir <- NULL
  if (save_png) {
    frames_dir <- sprintf("sh_frames_r%.3f_g%.3f", rparam, gparam)
    if (!dir.exists(frames_dir)) dir.create(frames_dir, recursive = TRUE)
  }
  if (save_png) {
    df <- as.data.frame(as.table(u)); colnames(df) <- c("y","x","val")
    p <- ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
      scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
      ggtitle(sprintf("t=%.2f  (r=%.2f, g=%.2f)", 0.0, rparam, gparam))
    ggsave(file.path(frames_dir, sprintf("frame_%05d.png", 0)), p, width = 4.5, height = 4.5, dpi = 120)
  }
  for (s in 1:steps) {
    u <- step_SH(u, rparam, gparam, dt)
    if (save_png && (s %% sample_every == 0)) {
      df <- as.data.frame(as.table(u)); colnames(df) <- c("y","x","val")
      p <- ggplot(df, aes(x = x, y = y, fill = val)) +
        geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
        scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026")) +
        ggtitle(sprintf("t=%.2f  (r=%.2f, g=%.2f)", s*dt, rparam, gparam))
      ggsave(file.path(frames_dir, sprintf("frame_%05d.png", s)), p, width = 4.5, height = 4.5, dpi = 120)
    }
  }
  u
}

# ------------------------- SWEEP -------------------------
results <- list()
thumb_records <- list()

for (rr in RG_R) {
  for (ggg in RG_G) {
    u_final <- simulate_SH(rr, ggg, steps = steps, dt = dt, u0 = NULL, save_png = save_samples)
    # Diagnostics
    uh <- fft2(u_final); S <- Mod(uh)^2
    dk <- dominant_k(S)
    kstar <- dk$kstar
    ring <- angular_spectrum(uh, kstar, halfwidth_frac = ring_halfwidth_frac, nbin = n_theta_bins)
    S2 <- fourier_mode_score(ring$theta, ring$F, m = 2)
    H6 <- fourier_mode_score(ring$theta, ring$F, m = 6)
    A  <- anisotropy_index(ring$F)
    skew <- field_skewness(u_final)
    Lscale <- if (is.finite(kstar) && kstar > 0) 2*pi/kstar else NA_real_
    cls <- classify_pattern(S2, H6, A, skew, class_rules)
    
    results[[length(results)+1]] <- data.frame(
      r = rr, g = ggg,
      kstar = kstar, L = Lscale,
      S2 = S2, H6 = H6, A = A, skew = skew,
      class = cls
    )
    
    # Thumbnail
    df <- as.data.frame(as.table(u_final)); colnames(df) <- c("y","x","val")
    p <- ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_raster() + coord_fixed() + theme_void() + theme(legend.position = "none") +
      scale_fill_gradientn(colors=c("#313695","#74add1","#e0f3f8","#fee090","#f46d43","#a50026"))
    thumb_fn <- file.path(thumb_dir, sprintf("thumb_r%.3f_g%.3f.png", rr, ggg))
    ggsave(thumb_fn, p, width = 2.2, height = 2.2, dpi = 120)
    thumb_records[[length(thumb_records)+1]] <- data.frame(r = rr, g = ggg, thumb = thumb_fn, class = cls)
    message(sprintf("[DONE] r=%.3f g=%.3f  => class=%s  (k*=%.3f, S2=%.2f, H6=%.2f, A=%.2f, skew=%.2f)",
                    rr, ggg, cls, kstar, S2, H6, A, skew))
  }
}

res_df <- do.call(rbind, results)
write.csv(res_df, csv_out, row.names = FALSE)
message("Saved metrics: ", csv_out)

# ---------------------- PHASE MAP PLOT --------------------
# Map classes to colors; place labels on grid
classes <- sort(unique(res_df$class))
pal <- c(
  stripes = "#1f78b4",
  hexagons = "#33a02c",
  hex_up = "#33a02c",
  hex_down = "#b2df8a",
  spots_up = "#e31a1c",
  spots_down = "#fb9a99",
  labyrinth = "#6a3d9a"
)
res_df$color <- pal[res_df$class]

# Arrange into image-like grid (r on x-axis, g on y-axis)
phase_plot <- ggplot(res_df, aes(x = factor(r), y = factor(g), fill = class)) +
  geom_tile(color = "grey30") +
  scale_fill_manual(values = pal) +
  labs(x = "r", y = "g", fill = "class", title = "Swift–Hohenberg Phase Map (2D)") +
  theme_minimal(base_size = 12)

ggsave(phase_png, phase_plot, width = 7, height = 5, dpi = 140)
message("Saved phase map: ", phase_png)

# -------------------- PRINT QUICK SUMMARY -----------------
print(res_df)
