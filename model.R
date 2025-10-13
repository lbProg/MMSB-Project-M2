# Libraries --------------------------------------------------------------

library(ggplot2)

# Functions --------------------------------------------------------------

theta <- function(t) {
  y <- sin(((2 * pi * t) / T_period) - (pi / 2))

  y > 0
}

growth_rate <- function(t, vars) {
  r * theta(t) * (vars$nutrients**2 / (alpha**2 + vars$nutrients**2))
}

grazing <- function(t, vars) {
  g * theta(t) * (1 - exp(-I * vars$phyto**2))
}

nutrients <- function(t, vars) {
  -vars$R * vars$phyto
}

phyto <- function(t, vars) {
  0.001
}

zoo <- function(t, vars) {
  0.001
}

L_ZN <- function(t, vars) {
  theta(t) * sigma + physical_control(-t) * gamma
}

# Model parameters -------------------------------------------------------

t_0 <- 0
t_f <- 365
dt <- 1

time_seq <- seq(t_0, t_f, dt)

# Constants --------------------------------------------------------------

T_period <- 365
r <- 1
alpha <- 0.1
I <- 1.2
g <- 0.5

# Initialize tracked variables -------------------------------------------
zeros <- rep(0, length(time_seq))

vars <- data.frame(
  "time" = seq(t_0, t_f, dt),
  "theta" = zeros,
  "R" = zeros,
  "nutrients" = zeros,
  "zoo" = zeros,
  "phyto" = zeros,
  "detritus" = zeros
)

vars$nutrients[1] <- 0
vars$zoo[1] <- 1
vars$phyto[1] <- 1
vars$detritus[1] <- 1

# Main model loop --------------------------------------------------------

for (i in 2:length(time_seq)) {
  t <- vars$time[i]
  vars_before <- vars[i - 1, ]

  vars$theta[i] <- theta(t)
  vars$R[i] <- growth_rate(t, vars_before)
  vars$G[i] <- grazing(t, vars_before)

  dNdt <- nutrients(t, vars_before)
  dZdt <- zoo(t, vars_before)
  dPdt <- phyto(t, vars_before)

  vars$nutrients[i] <- vars$nutrients[i - 1] + dNdt * dt
  vars$zoo[i] <- vars$zoo[i - 1] + dZdt * dt
  vars$phyto[i] <- vars$phyto[i - 1] + dPdt * dt
  # vars$detritus[i] <- vars$detritus[i - 1] + detritus(t) * dt
}

# Plot results -----------------------------------------------------------

par(mfrow = c(2, 2))
plot(vars$time, vars$theta, type = "l", col = "black")
plot(vars$time, vars$nutrients, type = "l", col = "red")

plot(vars$time, vars$R, type = "l", col = "black")
plot(vars$time, vars$G, type = "l", col = "black")