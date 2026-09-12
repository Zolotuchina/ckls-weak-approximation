library(ggplot2)
set.seed(42)

moment_order <- 2        # p in E[X_t^p]; exact coefficients available for p = 2, 3, 7
S0    <- 1               # initial value
sigma <- 1               # volatility
alpha <- 5/6             # elasticity; the second-order formulas below are derived for 5/6 only
switch_const  <- 2       # near-zero scheme is used when S < switch_const * sigma^6 * h^3
T_max <- 5               # time horizon

eval_times <- seq(0.2, T_max, by = 0.2)   # times at which E[X_t^p] is estimated
n_eval     <- length(eval_times)

h       <- 0.2           # time step
n_paths <- 1000000       # number of Monte Carlo trajectories

path_cap <- 1e12         # trajectories are capped here to avoid overflow


# Exact moments E[X_t^p] = sum (c/d) * (sigma^2 t)^a * x0^b, one row c(a, b, c, d) per term
moment_coeffs <- list(
  "2" = list(
    c(0,    2,       1,         1),
    c(1,  5/3,       1,         1),
    c(2,  4/3,       5,        18),
    c(3,    1,       5,       243)
  ),
  "3" = list(
    c(0,    3,       1,         1),
    c(1,  8/3,       3,         1),
    c(2,  7/3,      10,         3),
    c(3,    2,     140,        81),
    c(4,  5/3,      35,        81),
    c(5,  4/3,      35,       729),
    c(6,    1,      35,     19683)
  ),
  "7" = list(
    c(0,    7,              1,              1),
    c(1, 20/3,             21,              1),
    c(2, 19/3,            595,              3),
    c(3,    6,          90440,             81),
    c(4, 17/3,         113050,             27),
    c(5, 16/3,        2690590,            243),
    c(6,    5,      139910680,           6561),
    c(7, 14/3,      199872400,           6561),
    c(8, 13/3,     1923771850,          59049),
    c(9,    4,   125045170250,        4782969),
    c(10,11/3,    25009034050,        1594323),
    c(11,10/3,   100036136200,       14348907),
    c(12,   3,   875316191750,      387420489),
    c(13, 8/3,    67332014750,      129140163),
    c(14, 7/3,    96188592500,     1162261467),
    c(15,   2,   269328059000,    31381059609),
    c(16, 5/3,    33666007375,    62762119218),
    c(17, 4/3,     9901766875,   564859072962),
    c(18,   1,     9901766875, 45753584909922)
  )
)

# Exact E[X_t^p] on a vector of times
analytical_moment <- function(p, x0, t_vals) {
  key    <- as.character(p)
  w      <- sigma^2 * t_vals
  table  <- moment_coeffs[[key]]
  result <- rep(0, length(t_vals))
  for (row in table) {
    result <- result + (row[3]/row[4]) * w^row[1] * x0^row[2]
  }
  result
}

# Two-point variable with given mean and variance (first-order scheme)
sample_two_point <- function(m1, variance) {
  sd_val <- sqrt(pmax(variance, 0))
  ifelse(runif(length(m1)) < 0.5, m1 - sd_val, m1 + sd_val)
}

# Three-point variable on {0, z1, z2} matching moments m1..m4 (second-order scheme)
sample_three_point <- function(m1, m2, m3, m4) {
  n   <- length(m1)
  eps <- 1e-12
  
  denom <- m1*m3 - m2^2
  valid <- abs(denom) > eps
  
  s    <- ifelse(valid, (m1*m4 - m2*m3)/denom, m1)
  q    <- ifelse(valid, (m2*m4 - m3^2)/denom,  m1^2)
  disc <- pmax(s^2 - 4*q, 0)
  
  z1 <- (s - sqrt(disc))/2
  z2 <- (s + sqrt(disc))/2
  
  dz  <- pmax(abs(z2 - z1), eps)
  sz1 <- ifelse(abs(z1) > eps, z1, sign(z1 - eps/2)*eps)
  sz2 <- ifelse(abs(z2) > eps, z2, eps)
  
  p1 <- ifelse(valid, (m1*z2 - m2)/(sz1*dz), 0)
  p2 <- ifelse(valid, (m2 - m1*z1)/(sz2*dz), 0)
  p0 <- 1 - p1 - p2
  
  p1  <- pmax(0, pmin(1, p1))
  p2  <- pmax(0, pmin(1, p2))
  p0  <- pmax(0, 1 - p1 - p2)
  tot <- p0 + p1 + p2
  p0  <- p0/tot; p1 <- p1/tot; p2 <- p2/tot
  
  u <- runif(n)
  ifelse(!valid, 0,
         ifelse(u <= p0, 0,
                ifelse(u <= p0 + p1, z1, z2)))
}


# Euler-Maruyama step, reflected at zero
euler_step <- function(S, h) {
  dW <- rnorm(length(S), 0, sqrt(h))
  S_new <- S + sigma*pmax(S,0)^alpha*dW
  pmin(pmax(S_new, 0), path_cap)
}

# Step based on raw moments, used when S is close to zero
ckls_step_near_zero <- function(x, h, order) {
  eps   <- 1e-10
  x_pos <- pmax(x, eps)
  
  if (order == 1) {
    m1   <- x_pos
    m2   <- x_pos^2 + h*sigma^2*x_pos^(2*alpha)
    draw <- sample_two_point(m1, m2 - m1^2)
  } else {
    m1 <- x_pos
    m2 <- (5/18)*h^2*sigma^4*x_pos^(4/3) +
      h*sigma^2*x_pos^(5/3) + x_pos^2
    m3 <- (10/3)*h^2*sigma^4*x_pos^(7/3) +
      3*h*sigma^2*x_pos^(8/3) + x_pos^3
    m4 <- (44/3)*h^2*sigma^4*x_pos^(10/3) +
      6*h*sigma^2*x_pos^(11/3) + x_pos^4
    draw <- sample_three_point(m1, m2, m3, m4)
  }
  
  pmax(ifelse(x < eps, x, draw), 0)
}

# Step based on central moments of S_h - x, used away from zero
ckls_step_regular <- function(x, h, order) {
  x_pos <- pmax(x, 1e-10)
  
  if (order == 1) {
    m2   <- sigma^2*x_pos^(2*alpha)*h
    draw <- sample_two_point(rep(0, length(x)), m2)
  } else {
    m1 <- rep(0, length(x))
    m2 <- (5/18)*h^2*sigma^4*x_pos^(4/3) + h*sigma^2*x_pos^(5/3)
    m3 <- (5/2)*h^2*sigma^4*x_pos^(7/3)
    m4 <- 3*h^2*sigma^4*x_pos^(10/3) + (25/4)*h^3*sigma^6*x_pos^2
    draw <- sample_three_point(m1, m2, m3, m4)
  }
  
  pmax(x_pos + draw, 0)
}

# One CKLS step: picks the near-zero or regular scheme by the threshold
ckls_step <- function(S, h, order) {
  threshold <- switch_const * sigma^6 * h^3
  near_zero <- S < threshold
  
  S_new <- ifelse(near_zero,
                  ckls_step_near_zero(S, h, order),
                  ckls_step_regular(S, h, order))
  S_new <- pmax(S_new, 0)
  pmin(S_new, path_cap)
}

ckls_first_order_step <- function(S, h) ckls_step(S, h, order = 1)
ckls_second_order_step <- function(S, h) ckls_step(S, h, order = 2)


# Simulates n_traj trajectories and returns the empirical E[X_t^p] on eval_times
simulate_moment <- function(n_traj, h, step_fn) {
  n_steps <- round(T_max / h)
  h_adj   <- T_max / n_steps
  S       <- rep(S0, n_traj)
  moment  <- rep(NA_real_, n_eval)
  eval_idx    <- 1L
  
  for (i in seq_len(n_steps)) {
    S <- step_fn(S, h_adj)
    
    bad <- !is.finite(S)
    if (any(bad)) S[bad] <- 0
    
    t_now <- i * h_adj
    if (eval_idx <= n_eval &&
        abs(t_now - eval_times[eval_idx]) < h_adj*0.5 + 1e-10) {
      moment[eval_idx] <- mean(S^moment_order)
      eval_idx <- eval_idx + 1L
    }
  }
  moment
}

# Same, in chunks of 1e5 trajectories to limit memory
simulate_moment_chunked <- function(n_traj, h, step_fn, chunk = 1e5) {
  if (n_traj <= chunk) return(simulate_moment(n_traj, h, step_fn))
  
  n_chunks <- ceiling(n_traj / chunk)
  acc      <- matrix(0, nrow = n_eval, ncol = n_chunks)
  sizes    <- integer(n_chunks)
  for (k in seq_len(n_chunks)) {
    this_n   <- min(chunk, n_traj - (k-1)*chunk)
    sizes[k] <- this_n
    acc[, k] <- simulate_moment(this_n, h, step_fn)
  }
  as.numeric(acc %*% sizes / sum(sizes))
}


message("Euler (pure diffusion) ...")
t0 <- proc.time()
mom_euler <- simulate_moment_chunked(n_paths, h, euler_step)
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

message("CKLS 1st order (pure diffusion) ...")
t0 <- proc.time()
mom_first <- simulate_moment_chunked(n_paths, h, ckls_first_order_step)
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

message("CKLS 2nd order (pure diffusion) ...")
t0 <- proc.time()
mom_second <- simulate_moment_chunked(n_paths, h, ckls_second_order_step)
message(sprintf("  %.1f s", (proc.time()-t0)[3]))

mom_exact <- analytical_moment(moment_order, S0, eval_times)

# Root mean square error against the exact values
rmse <- function(e, r) sqrt(mean((e-r)^2, na.rm=TRUE))

results <- data.frame(
  T            = eval_times,
  Analytical   = round(mom_exact,  2),
  Euler        = round(mom_euler, 2),
  CKLS_1st_ord = round(mom_first, 2),
  CKLS_2nd_ord = round(mom_second, 2)
)

cat("\n=== Moment E[X^", moment_order, "_t] ===\n", sep="")
print(results, row.names = FALSE)

cat("\nRMSE vs analytical true values:\n")
cat(sprintf("  Euler             : %.4e\n", rmse(mom_euler, mom_exact)))
cat(sprintf("  CKLS first order  : %.4e\n", rmse(mom_first, mom_exact)))
cat(sprintf("  CKLS second order : %.4e\n", rmse(mom_second, mom_exact)))


moment_label <- switch(as.character(moment_order),
                       "2" = expression(E~"["~X[t]^2~"]"),
                       "3" = expression(E~"["~X[t]^3~"]"),
                       "7" = expression(E~"["~X[t]^7~"]"),
                       expression(E~"["~X[t]^p~"]"))

method_labels <- c("Analytical values", "Euler",
                   "CKLS first order", "CKLS second order")

df <- data.frame(
  T      = rep(eval_times, 4),
  Value  = c(mom_exact, mom_euler, mom_first, mom_second),
  Method = factor(rep(method_labels, each = n_eval),
                  levels = method_labels)
)

pal <- c("Analytical values" = "#000000",
         "Euler"              = "#E41A1C",
         "CKLS first order"     = "#377EB8",
         "CKLS second order"     = "#4DAF4A")

override_lt <- setNames(c(0, 1, 1, 1), method_labels)
override_sh <- setNames(c(16, NA, NA, NA), method_labels)

p <- ggplot(df, aes(x = T, y = Value, colour = Method)) +
  geom_line(data  = subset(df, Method != "Analytical values"),
            linewidth = 1.0) +
  geom_point(data = subset(df, Method == "Analytical values"),
             size = 2.5, shape = 16) +
  scale_colour_manual(name = NULL, values = pal, breaks = method_labels) +
  guides(colour = guide_legend(
    override.aes = list(
      linetype = override_lt[method_labels],
      shape    = override_sh[method_labels]
    )
  )) +
  labs(title    = sprintf("Weak approximations of the CKLS diffusion, E[X_t^%d]", moment_order),
       subtitle = bquote(S[0]==.(S0) ~","~~ sigma==.(sigma) ~","~~
                           alpha==.(format(alpha, digits=6)) ~","~~
                           h==.(h) ~","~~ n==.(n_paths) ~","~~
                           p==.(moment_order)),
       x        = "Time  t",
       y        = moment_label) +
  theme_bw(base_size = 13) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.text      = element_text(size = 11),
        plot.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(p)

dir.create("figures", showWarnings = FALSE)
out_file <- sprintf("figures/comparison_p%d_n%s.png", moment_order, format(n_paths, scientific = FALSE))
ggsave(out_file, p, width = 9, height = 6, dpi = 150)
message("Saved ", out_file)
