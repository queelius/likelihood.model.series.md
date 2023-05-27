# Load required packages
library(MASS)
library(mvtnorm)
library(algebraic.mle)

rate <- c(3,4,5)
rate.hat <- md_mle_exp_series_C1_C2_C3(exp_series_md_1,theta0=rate,eps=1e-50,eta=.005)
# Set the mean and covariance matrix for the 3D multivariate normal distribution
mu <- point(rate.hat)
sigma <- vcov(rate.hat)
# Define a grid for the contour plots
x_grid <- seq(1.2, 6.5, length.out = 200)
y_grid <- seq(1.2, 6.5, length.out = 200)
z_grid <- seq(1.2, 6.5, length.out = 200)

# Create a grid for each pair of variables
xy_grid <- expand.grid(x = x_grid, y = y_grid)
xz_grid <- expand.grid(x = x_grid, z = z_grid)
yz_grid <- expand.grid(y = y_grid, z = z_grid)

# Compute the probability density function (PDF) for each grid
xy_pdf <- matrix(dmvnorm(xy_grid, mean = mu[1:2], sigma = sigma[1:2, 1:2]), nrow = length(x_grid))
xz_pdf <- matrix(dmvnorm(xz_grid, mean = mu[c(1, 3)], sigma = sigma[c(1, 3), c(1, 3)]), nrow = length(x_grid))
yz_pdf <- matrix(dmvnorm(yz_grid, mean = mu[2:3], sigma = sigma[2:3, 2:3]), nrow = length(y_grid))


# Save the plots to PNG files
png("contour_xy_exp_series_1.png")
contour(x_grid, y_grid, xy_pdf, main = expression(paste("Contour plot (", hat(lambda[1]), "-", hat(lambda[2]), ")")), xlab = expression(hat(lambda[1])), ylab = expression(hat(lambda[2])))
dev.off()

png("contour_xz_exp_series_1.png")
contour(x_grid, z_grid, xz_pdf, main = expression(paste("Contour plot (", hat(lambda[1]), "-", hat(lambda[3]), ")")), xlab = expression(hat(lambda[1])), ylab = expression(hat(lambda[3])))
dev.off()
# Create contour plots for each pair of variables
png("contour_yz_exp_series_1.png")
contour(y_grid, z_grid, yz_pdf, main = expression(paste("Contour plot (", hat(lambda[2]), "-", hat(lambda[3]), ")")), xlab = expression(hat(lambda[2])), ylab = expression(hat(lambda[3])))
dev.off()
