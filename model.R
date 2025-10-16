# Libraries --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)

# Shape palette ----------------------------------------------------------

state_vars_palette <- c(
  "nutrients" = "solid",
  "zoo" = "dotted",
  "phyto" = "dashed",
  "detritus" = "twodash"
)

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

L_PD <- function(t, vars) {
  theta(y(t)) * vars$delta + theta(-y(t)) * lambda
}

L_ZN <- function(t) {
  theta(y(t)) * nu + theta(-y(t)) * gamma
}

L_ZD <- function(t) {
  theta(y(t)) * sigma + theta(-y(t)) * gamma
}

sinking_coef <- function(t, type) {
  if (type == 0) {
    r <- delta
  } else if (type == 1) {
    r <- delta_1 * exp(-((t - T_s)**2 / C_s))
  } else {
    r <- delta_2 * 0.5 * (1 + cos((2 * pi * t / T_period) - (pi / 2)))
  }

  r
}

sinking_rate <- function(t, vars) {
  theta(y(t)) * vars$delta
}

# Constants --------------------------------------------------------------

T_period <- 365
r <- 1
alpha <- 0.1
I <- 1.2
g <- 0.5
a <- 0.2

gamma <- 0.05
delta <- 0.02
delta_1 <- 0.035
delta_2 <- 0.02
eta <- 0.01
lambda <- 0.05
nu <- 0.01
sigma <- 0.02

Z_0 <- 0.01
P_0 <- 0.01

L_D <- 0

SN_ext <- 0
SD_ext <- 0

T_s <- 93
C_s <- 3600

# Initialize tracked variables -------------------------------------------

init_vars <- function(length, dt) {
  t_0 <- 0
  t_f <- length
  dt <- dt

  time_seq <- seq(t_0, t_f, dt)

  zeros <- rep(0, length(time_seq))

  vars <- data.frame(
    "time" = time_seq,
    "y" = zeros,
    "theta" = zeros,
    "R" = zeros,
    "G" = zeros,
    "delta" = zeros,
    "sinking_rate" = zeros,
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

  vars
}

# Main model loop function --------------------------------------------------------

run_model <- function(dt, vars, sinking_rate_type) {
  time_seq <- vars$time

  pb <- txtProgressBar(
    min = 0,
    max = length(time_seq),
    style = 3,
    width = 50,
    char = "="
  )
  
  for (i in 2:length(time_seq)) {
    t <- vars$time[i]
    vars_before <- vars[i - 1, ]

    vars$theta[i] <- theta(y(t))
    vars$R[i] <- growth_rate(t, vars_before)
    vars$G[i] <- grazing(t, vars_before)

    vars$L_PN[i] <- L_PN(t)
    vars$L_PD[i] <- L_PD(t, vars_before)
    vars$L_ZN[i] <- L_ZN(t)
    vars$L_ZD[i] <- L_ZD(t)

    vars$L_Z[i] <- vars$L_ZN[i] + vars$L_ZD[i]
    vars$L_P[i] <- vars$L_PN[i] + vars$L_PD[i]

    vars$delta[i] <- sinking_coef(t, sinking_rate_type)
    vars$sinking_rate[i] <- sinking_rate(t, vars_before)

    K <- compute_state_vars(t, vars_before)

    dNdt <- K[1]
    dZdt <- K[2]
    dPdt <- K[3]
    dDdt <- K[4]

    vars$nutrients[i] <- max(vars$nutrients[i - 1] + dNdt * dt, 0)
    vars$zoo[i] <- max(vars$zoo[i - 1] + dZdt * dt, 0)
    vars$phyto[i] <- max(vars$phyto[i - 1] + dPdt * dt, 0)
    vars$detritus[i] <- max(vars$detritus[i - 1] + dDdt * dt, 0)

    setTxtProgressBar(pb, i)
  }

  close(pb)

  vars
}

# Run the models ---------------------------------------------------------

length <- 400
dt <- 0.05

vars <- init_vars(length, dt)
model_base <- run_model(dt, vars, sinking_rate_type = 0)

vars <- init_vars(length, dt)
model_sinking_1 <- run_model(dt, vars, sinking_rate_type = 1)

vars <- init_vars(length, dt)
model_sinking_2 <- run_model(dt, vars, sinking_rate_type = 2)

# Plot results -----------------------------------------------------------

plot_model <- function(model) {
  model |>
    select(time, nutrients, zoo, phyto, detritus) |>
    pivot_longer(nutrients:detritus, names_to = "variable", values_to = "value") |>
    ggplot(aes(x = time, y = value, linetype = variable)) +
    geom_line() +
    scale_linetype_manual(values = state_vars_palette) +
    scale_x_continuous(limits = c(0, max(time)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2)) +
    labs(x = "NPZD", y = "t (days)", linetype = "Variable") +
    theme_test()
}

plot_model(model_base)

sinking_data <- data.frame(
  time = model_sinking_1$time,
  type_1 = model_sinking_1$sinking_rate,
  type_2 = model_sinking_2$sinking_rate
) |>
  pivot_longer(2:3, names_to = "type", values_to = "value")

ggplot(sinking_data, aes(x = time, y = value, linetype = type)) +
  geom_line() +
  theme_test()

plot_model(model_sinking_1)
