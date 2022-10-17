# library(CVXR)
library(metalog)
library(profvis)
y <- cdiac$annual
profvis({
    beta <- Variable(m)
    obj <- Minimize(0.5 * sum((y - beta)^2) + lambda * sum(pos(diff(beta))))
    prob <- Problem(obj)
    soln <- solve(prob, solver = "ECOS")
    betaHat <- soln$getValue(beta)
})

# Data etc
tab = metalog:::xtabulate(quakes$long);
x = tab$x;
y = cumsum(tab$n) / (sum(tab$n) + 1);
m = length(x);
n = 21;

# Formulation of minimization problem

# Build required matrices
Y = metalog:::Y_matrix(m, y, NULL, n);
# p = seq(0.00001, 0.99999, 0.00001) # Extreme!
p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
B = metalog:::LP_feas_coeffs(n, length(p), log(p / (1 - p)), p - 0.5, 1 / (p * (1 - p)))

# Specify problem
system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(sum((x - Y %*% a)^2))
    constraint1 = B %*% a >= 0
    problem = CVXR::Problem(objective, constraints = list(constraint1))
    result = CVXR::solve(problem, solver = "ECOS")
})


system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(sum((x - Y %*% a)^2))
    constraint1 = B %*% a >= 0
    problem = CVXR::Problem(objective, constraints = list(constraint1))
    prob_data = CVXR::get_problem_data(problem, solver = "ECOS")
    if (packageVersion("CVXR") > "0.99-7") {
        ECOS_dims = CVXR::ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
    } else {
        ECOS_dims = prob_data$data[["dims"]]
    }
    solver_output = ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                           G = prob_data$data[["G"]],
                                           h = prob_data$data[["h"]],
                                           dims = ECOS_dims,
                                           A = prob_data$data[["A"]],
                                           b = prob_data$data[["b"]])
    if (packageVersion("CVXR") > "0.99-7") {
        direct_soln = CVXR::unpack_results(problem, solver_output, prob_data$chain, prob_data$inverse_data)
    } else {
        direct_soln = CVXR::unpack_results(problem, "ECOS", solver_output)
    }
})

solution = result$getValue(a)[,1]
mlCV = mlog(solution, "coeff", check_validity = TRUE)
plot(mlCV)
qmlog(y, mlCV)


system.time(
    mlLP <- mlog(quakes$long, "data", n = n, method = "lp", check_validity = FALSE)
) # 0.4 s
solution
plot(mlLP)


curve(pmlog(x, mlCV), from = 160, to = 200, n = 2000, col = "red")
curve(pmlog(x, mlLP), from = 160, to = 200, n = 2000, col = "black", add = TRUE)


summary(rmlog(10000000, mlCV))
moment(mlCV)

x = quakes$long; y = "data"; mean = NULL; n = 21
x = c(-1, 0, 1); y = c(0.1, 0.5, 0.9); mean = NULL; n = 3
x = c(-1, 0, 1); y = c(0.1, 0.5, 0.9); mean = 0.01; n = 4

system.time(
    mlOLS <- mlog(x, y, n = n, mean = mean, method = "ols", check_validity = FALSE)
) # 0.004s ; 0.001

system.time(
    mlLP <- mlog(x, y, n = n, mean = mean, method = "lp", check_validity = FALSE)
) # 2.240s ; 0.003

system.time(
    mlNL <- mlog(x, y, n = n, mean = mean, method = "nl", check_validity = FALSE)
) # 10.674s; 0.014

system.time(
    mlCVX <- mlog(x, y, n = n, mean = mean, method = "cvx", check_validity = FALSE)
) # 0.810s; 0.330


plot(mlOLS)
plot(mlLP)
plot(mlNL)
plot(mlCVX)

qmlog(c(0.1, 0.5, 0.9), mlNL)
moment(mlNL)
qmlog(c(0.1, 0.5, 0.9), mlCVX)
moment(mlCVX)

plot(mlNL)
plot(mlCVX)





# Attempt to extend method to allow fitting the mean for bounded distributions.
# Let's start with an "easy" problem, an unbounded distribution, with 5% quantile
# at 1 and mean of 2.

x = 1;
y = 0.05;
mean = 2;
n = 2;
m = 1;

# Formulation of minimization problem

# Quantile fitting
Y = metalog:::Y_matrix(m, y, NULL, n);

# Mean fitting
# Simpson's rule: integral of f(x) from a to b, with n (even number) sub intervals of size h
# ~= (h/3)[ f(a) + 4 SUM_{j,1,n/2} f(x_{2j-1}) + 2 SUM_{j,1,n/2-1} f(x_{2j}) + f(b) ]
# so basically it's h/3 times the evaluated function from 0 to n with coefficients
# 1, 4, 2, 4, 2, 4, 2, 4, 2, 4, 1
q = seq(1e-4, 1-1e-4, length.out = 257); # Q for quadrature
h = q[2] - q[1];
simpsons = c(1, rep_len(c(4, 2), 255), 1);
Q = metalog:::Y_matrix(length(q), q, NULL, n) * simpsons;

# Constraint
p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
B = metalog:::LP_feas_coeffs(n, length(p), log(p / (1 - p)), p - 0.5, 1 / (p * (1 - p)))

# Specify problem
system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(sum((x - Y %*% a)^2) + (mean - h/3 * sum(Q %*% a))^2)
    constraint1 = B %*% a >= 0
    problem = CVXR::Problem(objective, constraints = list(constraint1))
    result = CVXR::solve(problem, solver = "ECOS")
})

result$num_iters
result$getValue(a)[,1]
solution = result$getValue(a)[,1]
mlCV = mlog(solution, "coeff", check_validity = TRUE)
plot(mlCV)
moment(mlCV) - 2



# Attempt 2, but using tanh-sinh quadrature
# Quantile fitting
Y = metalog:::Y_matrix(m, y, NULL, n);

# Mean fitting
# Tanh-sinh quadrature
N = 64
h = 1/16;
j = 0:N;
xj = 0.5/(exp(0.5 * pi * sinh(j * h)) * cosh(0.5 * pi * sinh(j * h)));
wj = (0.25 * pi * h * cosh(j * h)) / (cosh(0.5 * pi * sinh(j * h))^2);
Q = rbind(
    rev(wj) * metalog:::Y_matrix (length(j),     rev(xj), NULL, n),
    wj[-1]  * metalog:::Yc_matrix(length(j) - 1,  xj[-1], NULL, n)
)

# Constraint
p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
B = metalog:::LP_feas_coeffs(n, length(p), log(p / (1 - p)), p - 0.5, 1 / (p * (1 - p)))

# Specify problem
system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(sum((x - Y %*% a)^2) + (mean - sum(Q %*% a))^2)
    constraint1 = B %*% a >= 0
    problem = CVXR::Problem(objective, constraints = list(constraint1))
    result = CVXR::solve(problem)
})

result$getValue(a)[,1]
result$num_iters

solution = result$getValue(a)[,1]
mlCV = mlog(solution, "coeff", check_validity = TRUE)
plot(mlCV)
moment(mlCV) - 2




# Attempt 3: This time with a lower-bounded distribution
x = log(c(10, 12));

x = c(10, 12);
y = c(0.05, 0.95);
mean = 11;
n = 3;
m = 2;

# Quantile fitting
Y = metalog:::Y_matrix(m, y, NULL, n);

# Mean fitting
# Tanh-sinh quadrature
N = 64
h = 1/16;
j = 0:N;
xj = 0.5/(exp(0.5 * pi * sinh(j * h)) * cosh(0.5 * pi * sinh(j * h)));
wj = (0.25 * pi * h * cosh(j * h)) / (cosh(0.5 * pi * sinh(j * h))^2);
W = c(rev(wj), wj[-1])
Q = rbind(
    metalog:::Y_matrix (length(j),     rev(xj), NULL, n),
    metalog:::Yc_matrix(length(j) - 1,  xj[-1], NULL, n)
)

# Constraint
p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
B = metalog:::LP_feas_coeffs(n, length(p), log(p / (1 - p)), p - 0.5, 1 / (p * (1 - p)))

# Specify problem
system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(sum((x - Y %*% a)^2) + (mean - sum(W * exp(Q %*% a)))^2)
    constraint1 = B %*% a >= 0
    problem = CVXR::Problem(objective, constraints = list(constraint1))
    result = CVXR::solve(problem, solver = "ECOS", verbose = FALSE, num_iter = 1000)
})

result$getValue(a)[,1]
solution = result$getValue(a)[,1]
mlCV = mlog(solution, "coeff", check_validity = TRUE)
plot(mlCV)
moment(mlCV) - 11


# Smallest possible similar problem...
# Attempt 3: This time with a lower-bounded distribution
# Mean fitting
# Tanh-sinh quadrature
n = 1;
xj = 0.5;
wj = 1;
Q = 1;

wj = CVXR::Parameter(pos = TRUE, value = 1)

# Specify problem
system.time({
    a = CVXR::Variable(n)
    objective = CVXR::Minimize(exp(wj * Q * a)^2)
    problem = CVXR::Problem(objective)
    result = CVXR::solve(problem)
})

result$getValue(a)[,1]
result$num_iters

solution = result$getValue(a)[,1]
mlCV = mlog(solution, "coeff", check_validity = TRUE)
plot(mlCV)
moment(mlCV) - 2




# Tanh-sinh quadrature over the interval [0,1]
tsq = function(f, N = 128, h = 1/32, fc = function(x) f(1 - x))
{
    # Implementation guided by DH Bailey, Tanh-Sinh High-Precision Quadrature,
    # https://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf.
    # References in comments here are to this paper.

    # Values of j for sum in eq. (1), but only from 0 to N
    j = 0:N;

    # Complement of x_j, y_j = 1 - x_j, referred to in section 3
    y = 1/(exp(0.5 * pi * sinh(j * h)) * cosh(0.5 * pi * sinh(j * h)));

    # Weights w_j
    w = (0.5 * pi * h * cosh(j * h)) / (cosh(0.5 * pi * sinh(j * h))^2);

    # Conduct integration in two halves using f and fc
    # Left half, including j = 0
    lh = sum(w * f(y * 0.5));

    # Right half, excluding j = 0
    rh = sum(w[-1] * fc(y[-1] * 0.5));

    return (0.5 * (rh + lh))
}



# Legendre quadrature . . .
rules = gaussquad::legendre.quadrature.rules(100)
a = c(0, 2, 1, 0)

moment(list(a = a, b = c(0, 1)))
gaussquad::legendre.quadrature(metalog:::metalog_M, rules[[100]], 0, 1,
    a = a, n = length(a), bl = 0, bu = 1)

q = seq(1e-4, 1-1e-4, length.out = 101); # Q for quadrature
h = q[2] - q[1];
simpsons = c(1, rep_len(c(4, 2), 99), 1);
Q = metalog:::Y_matrix(length(q), q, NULL, length(a)) * simpsons;
M = Q %*% a
h/3 * sum(M)
eM = exp(M);
M = (0 + 1*eM) / (1 + eM);
h/3 * sum(M)


# Tanh-sinh quadrature
a = c(0, 1.1)
h = 0.1
k = -30:30

# Tanh-sinh quadrature
f = function(x) metalog:::metalog_M(x, a = c(1, 0.9, -8.2), n = 3, bl = -Inf, bu = Inf)
fc = function(x) metalog:::metalog_Mc(x, a = c(1, 0.9, -8.2), n = 3, bl = -Inf, bu = Inf)

# Tanh-sinh quadrature over the interval [0,1]
tsq = function(f, N = 128, h = 1/32, fc = function(x) f(1 - x))
{
    # Implementation guided by DH Bailey, Tanh-Sinh High-Precision Quadrature,
    # https://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf.
    # References in comments here are to this paper.

    # Values of j for sum in eq. (1), but only from 0 to N
    j = 0:N;

    # Complement of x_j, y_j = 1 - x_j, referred to in section 3
    y = 1/(exp(0.5 * pi * sinh(j * h)) * cosh(0.5 * pi * sinh(j * h)));

    # Weights w_j
    w = (0.5 * pi * h * cosh(j * h)) / (cosh(0.5 * pi * sinh(j * h))^2);

    # Conduct integration in two halves using f and fc
    # Left half, including j = 0
    lh = sum(w * f(y * 0.5));

    # Right half, excluding j = 0
    rh = sum(w[-1] * fc(y[-1] * 0.5));

    return (0.5 * (rh + lh))
}

tsq(f = f, fc = fc)
moment(list(a = c(1, 0.9, -8.2), b = c(-Inf, Inf)))

xk = tanh(0.5 * pi * sinh(k*h))
yk = 1/(exp(0.5 * pi * sinh(k*h)) * cosh(0.5 * pi * sinh(k*h)))
xk
yk
wk = (0.5*h*pi*cosh(k*h))/(cosh(0.5*pi*sinh(k*h))^2)
0.5*sum(wk * metalog:::metalog_M((xk + 1) / 2, a = a, n = length(a), bl = 0, bu = 1))

moment(list(a = a, b = c(0, 1)))

plot(mlog(a, "coeff", bounds = c(0, 1)))

