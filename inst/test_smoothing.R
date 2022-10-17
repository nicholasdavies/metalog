# An alternative approach to fitting empirical data with a smoothing spline

# The data
tab = metalog:::xtabulate(quakes$long)
x = tab$x
y = cumsum(tab$n)
yn = (y - 0.5) / (max(y))
yy = log(yn/(1-yn))

# x = 1:100
# y = c(rep(0,50), rep(1,50))

# Fit the spline
n = length(x);
h = diff(x);

Delta = matrix(0, nrow = n - 2, ncol = n);
W = matrix(0, nrow = n - 2, ncol = n - 2);
for (i in 1:(n-2)) {
    Delta[i, i] = 1/h[i];
    Delta[i, i+1] = -1/h[i] - 1/h[i+1];
    Delta[i, i+2] = 1/h[i+1];

    W[i-1, i] = h[i]/6;
    W[i, i-1] = h[i]/6;
    W[i, i] = (h[i] + h[i+1])/3;
}
A = t(Delta) %*% solve(W) %*% Delta

# lambda = 100 * dbeta(seq(0, 1, length.out = ncol(A)), 2, 2)
# lambda = c(rep(100, 100), rep(0, 405), rep(100, 100))
# lambda = runif(605, 0, 100)
# lambda = max(lambda) - lambda
lambda = 50
S = solve(diag(n) + lambda*A);
m = S %*% yy;
df = sum(diag(S))
df
plot(x, m[, 1], "l")
points(x, yy)


model = lm(y ~ poly(x,21), data = data.frame(x = x, y = y))
plot(model)
plot(predict(model, newdata = data.frame(x = seq(1, 100, by = 1))), type = "l")



d = data.frame(x, y)
fit.lo = loess(y ~ x, data = d, enp.target = 10)
plot(predict(fit.lo, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000))))

