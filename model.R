# Libraries --------------------------------------------------------------

library(ggplot2)

# Functions --------------------------------------------------------------

theta <- function(t) {
  y <- sin(((2 * pi * t) / T_period) - (pi / 2))

  y > 0
}

A_mix <- function(t) {
  theta(-t) * a
}

growth_rate <- function(t, vars) {
  r * theta(t) * (vars$nutrients**2 / (alpha**2 + vars$nutrients**2))
}

grazing <- function(t, vars) {
  g * theta(t) * (1 - exp(-I * vars$phyto**2))
}

nutrients <- function(t, vars) {
  -vars$R * vars$phyto *
    vars$L_ZN * (vars$zoo - Z_0) +
    vars$L_PN * (vars$phyto - P_0) +
    L_D * vars$detritus +
    A_mix(t) * (vars$detritus - vars$nutrients) +
    SN_ext
}

phyto <- function(t, vars) {
  vars$R * vars$phyto - vars$L_P * (vars$phyto - P_0) - vars$G * vars$zoo
}

zoo <- function(t, vars) {
  vars$G * vars$zoo - vars$L_Z * (vars$zoo - Z_0)
}

detritus <- function(t, vars) {
  vars$L_ZD * (vars$zoo - Z_0) +
    vars$L_PD * (vars$phyto - P_0) -
    L_D * vars$detritus - A_mix(t) * (vars$detritus - vars$nutrients) + SD_ext
}

L_PN <- function(t) {
  theta(t) * eta + theta(-t) * lambda
}

L_PD <- function(t) {
  theta(t) * delta + theta(-t) * lambda
}

L_ZN <- function(t, vars) {
  theta(t) * nu + theta(-t) * gamma
}

L_ZD <- function(t, vars) {
  theta(t) * sigma + theta(-t) * gamma
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

vars$nutrients[1] <- 1
vars$zoo[1] <- Z_0
vars$phyto[1] <- P_0
vars$detritus[1] <- 1

# Main model loop --------------------------------------------------------

for (i in 2:length(time_seq)) {
  t <- vars$time[i]
  vars_before <- vars[i - 1, ]

  vars$theta[i] <- theta(t)
  vars$R[i] <- growth_rate(t, vars_before)
  vars$G[i] <- grazing(t, vars_before)

  vars$L_PN[i] <- L_PN(t)
  vars$L_PD[i] <- L_PD(t)
  vars$L_ZN[i] <- L_ZN(t)
  vars$L_ZD[i] <- L_ZD(t)

  vars$L_Z[i] <- vars$L_ZN[i] + vars$L_ZD[i]
  vars$L_P[i] <- vars$L_PN[i] + vars$L_PD[i]

  dNdt <- nutrients(t, vars_before)
  dZdt <- zoo(t, vars_before)
  dPdt <- phyto(t, vars_before)
  dDdt <- detritus(t, vars_before)

  vars$nutrients[i] <- vars$nutrients[i - 1] + dNdt * dt
  vars$zoo[i] <- vars$zoo[i - 1] + dZdt * dt
  vars$phyto[i] <- vars$phyto[i - 1] + dPdt * dt
  vars$detritus[i] <- vars$detritus[i - 1] + dDdt * dt
}

# Plot results -----------------------------------------------------------

par(mfrow = c(3, 2))

plot(vars$time, vars$theta, type = "l", col = "black")
plot(vars$time, vars$R, type = "l", col = "black")
plot(vars$time, vars$G, type = "l", col = "black")

plot(vars$time, vars$nutrients, type = "l", col = "red", ylim = c(0, 2))
lines(vars$time, vars$zoo, type = "l", col = "blue")
lines(vars$time, vars$phyto, type = "l", col = "green")
lines(vars$time, vars$detritus, type = "l", col = "yellow")

plot(vars$time, vars$L_PN, type = "l")
lines(vars$time, vars$L_PD, type = "l")
lines(vars$time, vars$L_ZN, type = "l")
lines(vars$time, vars$L_ZD, type = "l")
