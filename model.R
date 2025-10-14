# Libraries --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)

# Functions --------------------------------------------------------------

y <- function(t) {
  sin(((2 * pi * t) / T_period) - (pi / 2))
}

theta <- function(y) {
  y > 0
}

A_mix <- function(t) {
  theta(-y(t)) * a
}

growth_rate <- function(t, vars) {
  r * theta(y(t)) * (vars$nutrients**2 / (alpha**2 + vars$nutrients**2))
}

grazing <- function(t, vars) {
  g * theta(y(t)) * (1 - exp(-I * vars$phyto**2))
}

nutrients <- function(t, vars) {
  -vars$R * vars$phyto +
    vars$L_ZN * (vars$zoo - Z_0) +
    vars$L_PN * (vars$phyto - P_0) +
    L_D * vars$detritus +
    A_mix(t) * (vars$detritus - vars$nutrients) +
    SN_ext
}

zoo <- function(t, vars) {
  vars$G * vars$zoo - vars$L_Z * (vars$zoo - Z_0)
}

phyto <- function(t, vars) {
  vars$R * vars$phyto - vars$L_P * (vars$phyto - P_0) - vars$G * vars$zoo
}

detritus <- function(t, vars) {
  vars$L_ZD * (vars$zoo - Z_0) +
  vars$L_PD * (vars$phyto - P_0) -
  L_D * vars$detritus - A_mix(t) * (vars$detritus - vars$nutrients) + SD_ext
}

compute_state_vars <- function(t, vars) {
  c(nutrients(t, vars), zoo(t, vars), phyto(t, vars), detritus(t, vars))
}

L_PN <- function(t) {
  theta(y(t)) * eta + theta(-y(t)) * lambda
}

L_PD <- function(t) {
  theta(y(t)) * delta + theta(-y(t)) * lambda
}

L_ZN <- function(t, vars) {
  theta(y(t)) * nu + theta(-y(t)) * gamma
}

L_ZD <- function(t, vars) {
  theta(y(t)) * sigma + theta(-y(t)) * gamma
}

# Model parameters -------------------------------------------------------

t_0 <- 0
t_f <- 365
dt <- 0.1

time_seq <- seq(t_0, t_f, dt)

pb <- txtProgressBar(
  min = 0,
  max = length(time_seq),
  style = 3,
  width = 50,
  char = "="
)

# Constants --------------------------------------------------------------

T_period <- 365
r <- 1
alpha <- 0.1
I <- 1.2
g <- 0.5
a <- 0.2

gamma <- 0.05
delta <- 0.02
eta <- 0.01
lambda <- 0.05
nu <- 0.01
sigma <- 0.02

Z_0 <- 0.01
P_0 <- 0.01

L_D <- 0

SN_ext <- 0
SD_ext <- 0

# Initialize tracked variables -------------------------------------------
zeros <- rep(0, length(time_seq))

vars <- data.frame(
  "time" = seq(t_0, t_f, dt),
  "y" = zeros,
  "theta" = zeros,
  "R" = zeros,
  "G" = zeros,
  "L_PN" = zeros,
  "L_PD" = zeros,
  "L_ZN" = zeros,
  "L_ZD" = zeros,
  "L_Z" = zeros,
  "L_D" = zeros,
  "L_P" = zeros,
  "nutrients" = zeros,
  "zoo" = zeros,
  "phyto" = zeros,
  "detritus" = zeros
)

vars$nutrients[1] <- 0.99
vars$zoo[1] <- Z_0
vars$phyto[1] <- P_0
vars$detritus[1] <- 0.99

# Main model loop --------------------------------------------------------

for (i in 2:length(time_seq)) {
  t <- vars$time[i]
  vars_before <- vars[i - 1, ]

  vars$theta[i] <- theta(y(t))
  vars$R[i] <- growth_rate(t, vars_before)
  vars$G[i] <- grazing(t, vars_before)

  vars$L_PN[i] <- L_PN(t)
  vars$L_PD[i] <- L_PD(t)
  vars$L_ZN[i] <- L_ZN(t)
  vars$L_ZD[i] <- L_ZD(t)

  vars$L_Z[i] <- vars$L_ZN[i] + vars$L_ZD[i]
  vars$L_P[i] <- vars$L_PN[i] + vars$L_PD[i]

  K1 <- compute_state_vars(t, vars_before)
  K2 <- compute_state_vars(t + dt / 2, vars_before + dt / 2 * K1)
  K3 <- compute_state_vars(t + dt / 2, vars_before + dt / 2 * K2)
  K4 <- compute_state_vars(t + dt, vars_before + dt * K3)

  K <- compute_state_vars(t, vars_before) + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4)

  dNdt <- K[1]
  dZdt <- K[2]
  dPdt <- K[3]
  dDdt <- K[4]

  # dNdt <- 0
  # dZdt <- 0
  # dPdt <- 0
  # dDdt <- 0

  vars$nutrients[i] <- max(vars$nutrients[i - 1] + dNdt * dt, 0)
  vars$zoo[i] <- max(vars$zoo[i - 1] + dZdt * dt, 0)
  vars$phyto[i] <- max(vars$phyto[i - 1] + dPdt * dt, 0)
  vars$detritus[i] <- max(vars$detritus[i - 1] + dDdt * dt, 0)

  setTxtProgressBar(pb, i)
}

close(pb)

# Plot results -----------------------------------------------------------

# par(mfrow = c(3, 2))

# plot(vars$time, vars$theta, type = "l", col = "black")
# plot(vars$time, vars$R, type = "l", col = "black")
# plot(vars$time, vars$G, type = "l", col = "black")

# plot(vars$time, vars$nutrients, type = "l", col = "red", ylim = c(0, 2))
# plot(vars$time, vars$zoo, type = "l", col = "blue")
# plot(vars$time, vars$phyto, type = "l", col = "green")
# plot(vars$time, vars$detritus, type = "l", col = "yellow")

# plot(vars$time, vars$L_Z, type = "l")
# plot(vars$time, vars$L_P, type = "l")

vars |>
  select(time, nutrients, zoo, phyto, detritus) |>
  pivot_longer(nutrients:detritus, names_to = "variable", values_to = "value") |>
  ggplot(aes(x = time, y = value, color = variable)) +
  geom_line() +
  theme_bw()

# plot(vars$time, vars$L_PN, type = "l")
# lines(vars$time, vars$L_PD, type = "l")
# lines(vars$time, vars$L_ZN, type = "l")
# lines(vars$time, vars$L_ZD, type = "l")
