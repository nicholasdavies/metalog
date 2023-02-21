
data = metalog:::xtabulate(quakes$long)
x = data$x
y_increase = (cumsum(data$n) - 0.5) / sum(data$n)
y_increase = log(y_increase/(1-y_increase))
plot(x, y_increase)

# Create `I-Spline` basis.
ispline <- splines2::iSpline(x, knots = quantile(x, 1:5/6), degree = 2, intercept = TRUE)

# Add constant to `I-Spline` basis.
ispline <- cbind(1, ispline)


# Creating `D` matrix and `d` vector.
d_mat <- crossprod(ispline, ispline)
d_vec <- crossprod(ispline, y_increase)

# Setting equality constraints s.t., all paramters are >= 0.
a_mat <- diag(1, ncol(ispline))
b_vec <- rep(0, ncol(ispline))

# Freeing the first parameter.
a_mat[1, 1] <- 0

# Solving the problem.
alpha_non_negative <- quadprog::solve.QP(Dmat = d_mat, dvec = d_vec, Amat = t(a_mat), bvec = b_vec)$solution

# Getting fitted values (i.e., the spline).
spline_non_decreasing <- ispline %*% alpha_non_negative

# Fit monotone non-increasing spline for decreasing toy data.
y_decrease <- rev(y_increase)

# Creating `D` matrix and `d` vector.
d_mat <- crossprod(ispline, ispline)
d_vec <- crossprod(ispline, y_decrease)

# Setting equality constraints s.t., all parameters are <= 0.
a_mat <- diag(-1, ncol(ispline))
b_vec <- rep(0, ncol(ispline))

# Freeing the first parameter.
a_mat[1, 1] <- 0

# Solving the problem.
alpha_non_positive <- quadprog::solve.QP(Dmat = d_mat, dvec = d_vec, Amat = t(a_mat), bvec = b_vec)$solution

# Getting fitted values (i.e., the spline).
spline_non_increasing <- ispline %*% alpha_non_positive



plot(x, y_increase, pch = 19, cex = 1.5, col = "darkgray", main = paste0("Non-decreasing | SSQ = ", round(sum((y_increase - spline_non_decreasing)^2), 5)), cex.main = 2)
lines(x, spline_non_decreasing, lwd = 2, col = "rosybrown")
points(x, spline_non_decreasing, pch = 19, cex = 1.5, col = "steelblue")

plot(x, y_decrease, pch = 19, cex = 1.5, col = "darkgray", main = paste0("Non-increasing | SSQ = ", round(sum((y_decrease - spline_non_increasing)^2), 5)), cex.main = 2)
lines(x, spline_non_increasing, lwd = 2, col = "rosybrown")
points(x, spline_non_increasing, pch = 19, cex = 1.5, col = "steelblue")


y2 = exp(spline_non_decreasing) / (1 + exp(spline_non_decreasing))
y0 = exp(y_increase) / (1 + exp(y_increase))
plot(x, y2, type = "l")
lines(x, y0, col = "red")

y_increase = (cumsum(data$n) - 0.5) / sum(data$n)
y_increase = log(y_increase/(1-y_increase))
