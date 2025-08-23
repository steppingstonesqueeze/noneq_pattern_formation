# Buckling beam Euler style
# app.R — Shiny app for beam buckling: Euler (linear), Elastica (post-buckling), and Pitchfork (normal form)
# Author: ChatGPT (for Girish)
# ---------------------------------------------------------------
# FEATURES
# - Linear Euler buckling: exact critical load P_cr for common boundary conditions via effective-length factor K.
#   * Pinned–Pinned, Clamped–Clamped, Clamped–Pinned, Cantilever (Clamped–Free).
#   * First-mode shape: exact for pinned–pinned; BC-satisfying illustratives for others (annotated).
# - Nonlinear post-buckling (perfect column, pinned–pinned) under displacement control using Elastica:
#   * End-shortening ratio δ = Δ/L; relation δ = 1 - E(m)/K(m), load P = EI * (2K(m)/L)^2.
#   * Shape via integrating θ'' + λ^2 sin θ = 0 with symmetry (RK4), then plotting centerline.
# - Pitchfork diagram (normal form) with imperfection ε: 0 = μ a - a^3 + ε; shows (im)perfect pitchfork & stability.
# - Clean UI with tabs and export options.
#
# NOTES
# - SI units internally: E [Pa], I [m^4], L [m], P [N]. UI converts E [GPa], I [cm^4].
# - Mode shapes for non P–P BCs are illustrative; P_cr values are exact.
# - Assumes slender inextensible rod (no shear deformation, no material nonlinearity).
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  ok <- requireNamespace("shiny", quietly = TRUE); if (!ok) stop("Install 'shiny'")
  ok <- requireNamespace("ggplot2", quietly = TRUE); if (!ok) stop("Install 'ggplot2'")
})

library(shiny)
library(ggplot2)

# ---------- Math helpers ----------
# Elliptic integrals K(m), E(m) with robust fallback
KE_from_m <- function(m) {
  stopifnot(m >= 0 && m < 1)
  # Try pracma::ellipke if available (fast), else numeric integrate fallback
  if (requireNamespace("pracma", quietly = TRUE)) {
    out <- try(pracma::ellipke(m), silent = TRUE)
    if (!inherits(out, "try-error")) {
      # pracma::ellipke returns either a named vector c(K, E) or list with $K, $E depending on version
      if (is.list(out)) {
        K <- out$K; E <- out$E
      } else {
        # choose greater as K (K >= E for m in [0,1))
        K <- max(out[1], out[2]); E <- min(out[1], out[2])
      }
      return(list(K = K, E = E))
    }
  }
  # Fallback: numerical integration (slower, but reliable)
  kintegrand <- function(phi) 1 / sqrt(1 - m * (sin(phi))^2)
  eintegrand <- function(phi) sqrt(1 - m * (sin(phi))^2)
  K <- integrate(kintegrand, 0, pi/2, rel.tol = 1e-9, subdivisions = 400L)$value
  E <- integrate(eintegrand, 0, pi/2, rel.tol = 1e-9, subdivisions = 400L)$value
  list(K = K, E = E)
}

# Find m in [0,1) from end-shortening ratio delta = Δ/L satisfying: delta = 1 - E(m)/K(m)
solve_m_from_delta <- function(delta) {
  if (delta <= 0) return(0)
  if (delta >= 0.999) delta <- 0.999  # guard
  f <- function(m) {
    ke <- KE_from_m(m)
    1 - ke$E/ke$K - delta
  }
  uniroot(f, interval = c(1e-12, 1 - 1e-10), tol = 1e-10)$root
}

# Simple RK4 integrator for elastica (center-to-end), returning x(s), y(s), θ(s)
elastica_halfshape <- function(L, m, n = 600) {
  # Parameters from elliptic theory
  ke <- KE_from_m(m)
  K <- ke$K
  lambda <- 2*K/L
  k <- sqrt(m)
  theta_max <- 2*asin(k)  # angle at midspan
  
  ds <- (L/2) / n
  s <- 0
  theta <- theta_max
  q <- 0    # θ'(s)
  x <- 0; y <- 0
  xs <- numeric(n+1); ys <- numeric(n+1); thetas <- numeric(n+1)
  xs[1] <- x; ys[1] <- y; thetas[1] <- theta
  
  derivs <- function(theta, q) {
    # θ' = q; q' = -λ^2 sin θ; x' = cos θ; y' = sin θ
    c(q, -lambda^2 * sin(theta))
  }
  
  for (i in 1:n) {
    # RK4 for theta and q
    k1 <- derivs(theta, q)
    k2 <- derivs(theta + 0.5*ds*k1[1], q + 0.5*ds*k1[2])
    k3 <- derivs(theta + 0.5*ds*k2[1], q + 0.5*ds*k2[2])
    k4 <- derivs(theta + ds*k3[1], q + ds*k3[2])
    
    theta_next <- theta + ds*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
    q_next     <- q     + ds*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
    
    # Update centerline by integrating x' = cos θ, y' = sin θ (midpoint)
    # Use simple RK4 for x,y too, slaved to θ path
    g <- function(th) c(cos(th), sin(th))
    g1 <- g(theta)
    g2 <- g(theta + 0.5*ds*k1[1])
    g3 <- g(theta + 0.5*ds*k2[1])
    g4 <- g(theta + ds*k3[1])
    x <- x + ds*(g1[1] + 2*g2[1] + 2*g3[1] + g4[1])/6
    y <- y + ds*(g1[2] + 2*g2[2] + 2*g3[2] + g4[2])/6
    
    s <- s + ds
    theta <- theta_next; q <- q_next
    xs[i+1] <- x; ys[i+1] <- y; thetas[i+1] <- theta
  }
  
  data.frame(s = seq(0, L/2, length.out = n+1), x = xs, y = ys, theta = thetas,
             lambda = lambda, K = K, k = k)
}

# Build full shape by mirroring the half-shape (center at s=0). Ends should lie on y=0 for hinged case.
build_full_shape <- function(half_df) {
  # Shift so that the end (last point) is at y = 0 (tiny numeric drift otherwise)
  y_end <- tail(half_df$y, 1)
  half_df$y <- half_df$y - y_end
  
  # Mirror about s=0 (x odd/even? position integrates from center -> end, so mirror by reversing and negating x and y symmetry)
  left <- half_df
  left$x <- -rev(half_df$x)
  left$y <-  rev(half_df$y)    # symmetric w.r.t. midspan (y even)
  left$s <- -rev(half_df$s)
  
  # Combine (ensure mid-point not duplicated)
  full <- rbind(left[-1,], transform(half_df, s = +s))
  # Normalize so that ends are at x = ±(Lx/2). Our integration yields chord length Lx ≈ 2 * max(x)
  Lx <- max(full$x) - min(full$x)
  list(df = full, Lx = Lx)
}

# Generate P(δ) curve (displacement control, perfect P–P) for plotting
curve_P_delta <- function(L, E, I, m_max = 0.999, n = 120) {
  ms <- seq(0, m_max, length.out = n)
  out <- lapply(ms, function(m) {
    ke <- KE_from_m(m)
    K <- ke$K; Eke <- ke$E
    delta <- 1 - Eke/K
    lambda <- 2*K/L
    P <- E * I * lambda^2
    c(delta, P)
  })
  mat <- do.call(rbind, out)
  df <- data.frame(delta = mat[,1], P = mat[,2])
  df <- df[order(df$delta), ]
  df
}

# Pitchfork normal form utilities: 0 = μ a - a^3 + ε
pitchfork_branches <- function(eps = 0, amax = 2, n = 600) {
  # Parametrize by amplitude a (excluding 0 when eps != 0 to avoid blow-up)
  a_vals <- seq(-amax, amax, length.out = n)
  a_vals <- a_vals[abs(a_vals) > 1e-6 | eps == 0]  # keep a=0 only if eps==0
  mu <- (a_vals^3 + eps)/a_vals
  # Stability from potential V = -1/2 μ a^2 + 1/4 a^4 - ε a; at equilibrium, stable if d^2V/da^2 = -μ + 3a^2 > 0
  stab <- ( -mu + 3*a_vals^2 ) > 0
  data.frame(a = a_vals, mu = mu, stable = stab)
}

# Effective length factors and illustrative mode shapes
K_factor_for_BC <- function(bc) {
  switch(bc,
         "Pinned–Pinned" = 1.0,
         "Clamped–Clamped" = 0.5,
         "Clamped–Pinned"  = 0.699,   # standard value
         "Cantilever (Clamped–Free)" = 2.0,
         1.0)
}

mode_shape_x_y <- function(bc, L, n = 400) {
  x <- seq(0, L, length.out = n)
  y <- switch(bc,
              "Pinned–Pinned" = sin(pi * x / L),
              # The following are BC-satisfying illustrative shapes (not exact eigenfunctions of the buckling ODE)
              "Clamped–Clamped" = 1 - cos(2*pi * x / L),  # satisfies y=0 and y'=0 at ends
              "Clamped–Pinned" = (x/L)^2 * (1 - x/L),     # y(0)=0, y'(0)=0, y(L)=0
              "Cantilever (Clamped–Free)" = (x/L)^2 * (3 - 2*x/L),  # y(0)=0,y'(0)=0; free-end qualitative
              sin(pi * x / L))
  # Normalize to unit max deflection
  y <- y / max(abs(y))
  data.frame(x = x, y = y)
}

format_force <- function(P) {
  if (P < 1e3) {
    sprintf("%.3f N", P)
  } else if (P < 1e6) {
    sprintf("%.3f kN", P/1e3)
  } else {
    sprintf("%.3f MN", P/1e6)
  }
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Beam Buckling: Euler, Elastica & Pitchfork"),
  sidebarLayout(
    sidebarPanel(
      h4("Material & Geometry (SI)"),
      numericInput("E_gpa", "Young's modulus E [GPa]", value = 200, min = 1, max = 500, step = 1),
      numericInput("I_cm4", "Area moment I [cm^4]", value = 100, min = 0.01, max = 1e6, step = 1),
      numericInput("L_m", "Length L [m]", value = 2, min = 0.1, max = 20, step = 0.1),
      br(),
      h4("Boundary Condition (linear buckling)"),
      selectInput("bc", "End conditions:",
                  choices = c("Pinned–Pinned", "Clamped–Clamped", "Clamped–Pinned", "Cantilever (Clamped–Free)"),
                  selected = "Pinned–Pinned"),
      br(),
      h4("Elastica (Post-buckling, P–P)"),
      sliderInput("delta_ratio", "End-shortening δ = Δ/L", min = 0, max = 0.45, value = 0.05, step = 0.005),
      sliderInput("elastica_n", "Shape resolution (points)", min = 200, max = 2000, value = 600, step = 50),
      br(),
      h4("Pitchfork (Normal Form)"),
      sliderInput("epsilon", "Imperfection ε", min = -0.05, max = 0.05, value = 0.0, step = 0.001),
      sliderInput("amax", "Amplitude range |a| ≤", min = 0.2, max = 3, value = 1.5, step = 0.1),
      br(),
      helpText("Units: E in GPa (converted to Pa), I in cm^4 (converted to m^4).",
               "Elastica tab assumes perfect pinned–pinned under displacement control.")
    ),
    mainPanel(
      tabsetPanel(id = "tabs", type = "pills",
                  tabPanel("Linear Buckling",
                           br(),
                           uiOutput("pcr_info"),
                           plotOutput("mode_plot", height = "360px"),
                           br(),
                           tableOutput("pcr_table"),
                           helpText("Mode shapes for non P–P cases are illustrative BC-satisfying polynomials/trigs; P_cr is exact.")
                  ),
                  tabPanel("Elastica (Post-buckling)",
                           br(),
                           uiOutput("elastica_info"),
                           plotOutput("elastica_shape", height = "380px"),
                           br(),
                           plotOutput("P_delta_plot", height = "320px")
                  ),
                  tabPanel("Pitchfork Diagram",
                           br(),
                           uiOutput("pitchfork_info"),
                           plotOutput("pitchfork_plot", height = "360px")
                  ),
                  tabPanel("About",
                           br(),
                           tags$div(style = "max-width:800px;",
                                    h4("Physics & Assumptions"),
                                    p("Euler buckling: P_cr = π^2 E I / (K L)^2 with effective-length factor K depending on end conditions."),
                                    tags$ul(
                                      tags$li("Pinned–Pinned: K = 1.0"),
                                      tags$li("Clamped–Clamped: K = 0.5"),
                                      tags$li("Clamped–Pinned: K ≈ 0.699"),
                                      tags$li("Cantilever (Clamped–Free): K = 2.0")
                                    ),
                                    p("Elastica (perfect column, P–P, displacement control): δ = 1 − E(m)/K(m), P = E I (2K(m)/L)^2 using complete elliptic integrals (parameter m = k²)."),
                                    p("Pitchfork: normal form 0 = μ a − a³ + ε shows symmetry-breaking via imperfection ε.")
                           )
                  )
      )
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session) {
  E  <- reactive({ input$E_gpa * 1e9 })     # Pa
  I  <- reactive({ input$I_cm4 * 1e-8 })    # m^4
  L  <- reactive({ input$L_m })             # m
  
  # --- Linear Buckling ---
  pcr_val <- reactive({
    Kfac <- K_factor_for_BC(input$bc)
    pi^2 * E() * I() / ( (Kfac * L())^2 )
  })
  
  output$pcr_info <- renderUI({
    Kfac <- K_factor_for_BC(input$bc)
    HTML(sprintf("<b>Critical load P<sub>cr</sub></b> = %s (K = %.3f)",
                 format_force(pcr_val()), Kfac))
  })
  
  output$mode_plot <- renderPlot({
    df <- mode_shape_x_y(input$bc, L())
    ggplot(df, aes(x = x, y = y)) +
      geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
      geom_line(linewidth = 1) +
      labs(x = "x [m]", y = "Mode shape (unit max)",
           title = sprintf("First buckling mode — %s", input$bc)) +
      theme_minimal(base_size = 13)
  })
  
  output$pcr_table <- renderTable({
    bcs <- c("Pinned–Pinned", "Clamped–Clamped", "Clamped–Pinned", "Cantilever (Clamped–Free)")
    Kf  <- sapply(bcs, K_factor_for_BC)
    Pcr <- pi^2 * E() * I() / ((Kf * L())^2)
    data.frame(BC = bcs, K = Kf, Pcr_N = round(Pcr, 2), Pcr_kN = round(Pcr/1e3, 3))
  })
  
  # --- Elastica (Post-buckling, P–P, displacement control) ---
  elastica_data <- reactive({
    req(input$delta_ratio >= 0, input$delta_ratio < 0.9)
    m <- solve_m_from_delta(input$delta_ratio)
    half <- elastica_halfshape(L(), m, n = input$elastica_n)
    built <- build_full_shape(half)
    df <- built$df
    # scale X so that chord length equals L * (1 - δ)
    Lx_target <- L() * (1 - input$delta_ratio)
    scale <- Lx_target / (max(df$x) - min(df$x))
    df$x <- df$x * scale
    
    # Midspan deflection (relative to chord y=0): center at s=0, which corresponds to middle row
    y_mid <- df$y[which.min(abs(df$s))]
    
    # Load from λ(m)
    ke <- KE_from_m(m)
    K <- ke$K
    lambda <- 2*K/L()
    P <- E() * I() * lambda^2
    
    list(df = df, P = P, m = m, y_mid = y_mid, lambda = lambda)
  })
  
  output$elastica_info <- renderUI({
    ed <- elastica_data()
    HTML(sprintf(
      paste0("<b>Elastica (P–P, displacement control)</b>: δ = %.4f, ",
             "m = %.5f, P = %s, midspan deflection (normalized units) y<sub>mid</sub> = %.4f"),
      input$delta_ratio, ed$m, format_force(ed$P), ed$y_mid))
  })
  
  output$elastica_shape <- renderPlot({
    ed <- elastica_data()
    df <- ed$df
    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
      coord_equal(xlim = c(min(df$x), max(df$x))) +
      labs(x = "x [m]", y = "y [m] (scaled)",
           title = "Post-buckled shape (Elastica), pinned–pinned, displacement control") +
      theme_minimal(base_size = 13)
  })
  
  output$P_delta_plot <- renderPlot({
    df <- curve_P_delta(L(), E(), I())
    # Mark current selection
    ed <- elastica_data()
    sel <- data.frame(delta = input$delta_ratio, P = ed$P)
    
    ggplot(df, aes(x = delta, y = P/1e3)) +
      geom_line(linewidth = 1) +
      geom_point(data = sel, aes(x = delta, y = P/1e3), size = 3) +
      labs(x = "End-shortening ratio δ = Δ/L", y = "Load P [kN]",
           title = "Load–end-shortening curve (perfect P–P)") +
      theme_minimal(base_size = 13)
  })
  
  # --- Pitchfork Diagram (normal form) ---
  output$pitchfork_info <- renderUI({
    HTML(sprintf("Normal form: 0 = μ a − a³ + ε. Here ε = %.4f (imperfection). Stable portions shown solid; unstable dashed.", input$epsilon))
  })
  
  output$pitchfork_plot <- renderPlot({
    df <- pitchfork_branches(eps = input$epsilon, amax = input$amax)
    # Determine stability for plotting style
    ggplot(df, aes(x = mu, y = a)) +
      geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
      geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted") +
      geom_path(data = df[df$stable, ], linewidth = 1.1) +
      geom_path(data = df[!df$stable, ], linewidth = 1.1, linetype = "dashed") +
      labs(x = expression(mu~" (control parameter)"), y = "Amplitude a",
           title = if (abs(input$epsilon) < 1e-10) "Perfect pitchfork (supercritical)" else "Imperfect pitchfork (symmetry broken)") +
      theme_minimal(base_size = 13)
  })
}

# ---------- Run ----------
shinyApp(ui, server)
