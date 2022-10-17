# Convex optimisation fitting
CVX_fit_metalog = function(x, y, mean, n, bounds)
{
    x = c(x, mean); # CDF points / moments to fit
    m = length(x);  # number of CDF points / moments to fit

    # Build required matrices
    Y = metalog:::Y_matrix(m, y, mean, n);
    # p = seq(0.00001, 0.99999, 0.00001) # Extreme!
    p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
    B = metalog:::LP_feas_coeffs(n, length(p), log(p / (1 - p)), p - 0.5, 1 / (p * (1 - p)));

    # Formulate and solve problem
    a = CVXR::Variable(n); # coefficients to solve for
    objective = CVXR::Minimize(sum((x - Y %*% a)^2)); # minimize squared error
    # objective = CVXR::Minimize(sum((x - Y %*% a)^2) + 0.1*sum(abs(a[3:n]))); # WITH REGULARIZATION
    constraint1 = B %*% a >= 0; # feasibility constraints
    problem = CVXR::Problem(objective, constraints = list(constraint1)); # formulate problem
    result = CVXR::solve(problem, solver = "ECOS"); # solve problem

    result$getValue(a)[,1]
}
