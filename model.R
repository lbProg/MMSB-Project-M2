# Libraries --------------------------------------------------------------

library(ggplot2)

# Functions --------------------------------------------------------------

physical_control <- function(t) {
  y <- sin(((2 * pi * t) / T_period) - (pi / 2))

  y > 0
}

growth_rate <- function(t, vars) {
  r * vars$theta * (vars$nutrients**2 / (alpha**2 + vars$nutrients**2))
}

grazing <- function(t, vars) {
  g * vars$theta * (1 - exp(-I * vars$phytoplankton**2))
}

nutrients <- function(t) {
  0.001
}

phytoplankton <- function(t) {
  0.001
}

zooplankton <- function(t) {
  0.001
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
  "zooplankton" = zeros,
  "phytoplankton" = zeros,
  "detritus" = zeros
)

vars$nutrients[1] <- 0
vars$zooplankton[1] <- 1
vars$phytoplankton[1] <- 1
vars$detritus[1] <- 1

# Main model loop --------------------------------------------------------

for (i in 2:length(time_seq)) {
  t <- vars$time[i]
  vars_before <- vars[i - 1, ]

  vars$theta[i] <- physical_control(t)
  vars$R[i] <- growth_rate(t, vars_before)
  vars$G[i] <- grazing(t, vars_before)

  dNdt <- nutrients(t)
  dZdt <- zooplankton(t)
  dPdt <- phytoplankton(t)

  vars$nutrients[i] <- vars$nutrients[i - 1] + nutrients(t) * dt
  vars$zooplankton[i] <- vars$zooplankton[i - 1] + zooplankton(t) * dt
  vars$phytoplankton[i] <- vars$phytoplankton[i - 1] + phytoplankton(t) * dt
  # vars$detritus[i] <- vars$detritus[i - 1] + detritus(t) * dt
}

# Plot results -----------------------------------------------------------

par(mfrow = c(2, 2))
plot(vars$time, vars$theta, type = "l", col = "black")
plot(vars$time, vars$nutrients, type = "l", col = "red")

plot(vars$time, vars$R, type = "l", col = "black")
plot(vars$time, vars$G, type = "l", col = "black")