// xtabulate.cpp

#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List xtabulate(NumericVector& y)
{
    NumericVector x = clone(y);
    IntegerVector n(x.size());

    if (x.size() > 0)
    {
        std::sort(x.begin(), x.end());

        auto counter = n.begin();
        auto first = x.begin(), last = x.end();
        auto dest = first;

        while (++first != last)
        {
            ++*counter;
	        if (*dest != *first)
	        {
	            *++dest = *first;
	            ++counter;
	        }
        }
        ++*counter;
        x.erase(dest + 1, last);
        n.erase(counter + 1, n.end());
    }

    List L = List::create(_["x"] = x, _["n"] = n);
    return L;
}
