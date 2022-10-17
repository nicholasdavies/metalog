#' @export
print.mlog = function(ml)
{
    # TODO should there in fact be any square brackets here?
    bounds = paste0(if (ml$b[1] == -Inf) "(" else "[", ml$b[1], ", ",
        ml$b[2], if (ml$b[2] == Inf) ")" else "]");
    cat(length(ml$a), "-term metalog distribution with bounds ", bounds,
        " fit using ", ml$method, "\n", sep = "");
    cat("Coefficients:", ml$a, "\n");
    cat("Contains elements:", names(ml), "\n");
}
