# Libraries --------------------------------------------------------------

library(ggplot2)

# Functions --------------------------------------------------------------

physical_control <- function(t) {
  y <- sin(((2 * pi * t) / T_period) - (pi / 2))

  y > 0
}

# Model parameters -------------------------------------------------------

t_0 <- 0
t_f <- 365
dt <- 1

# Constants --------------------------------------------------------------

T_period <- 365

# Initialize tracked variables -------------------------------------------

time_seq <- seq(t_0, t_f, dt)

theta <- rep(0, length(time_seq))

# Main model loop --------------------------------------------------------

for (i in 1:length(time_seq)) {
  t <- time_seq[i]

  theta[i] <- physical_control(t)
}

# Plot results -----------------------------------------------------------

plot(time_seq, theta, type = "l")
