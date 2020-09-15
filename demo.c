/*
    demo.c  -  Don Cross

    Demonstrates usage of Chebyshev polynomials to approximate a function.

    https://github.com/cosinekitty/chebyshev
*/

#include <math.h>
#include <stdio.h>
#include "chebyshev.h"


double TestFunc(const void *context, double t)
{
    (void)context;      /* This function doesn't need a context. */

    /* Demonstrate with a hyperbolic cosine, which models the shape of a hanging chain. */
    return (exp(t) + exp(-t)) / 2.0;
}


int main(int argc, const char *argv[])
{
    const int nsamples = 20;
    int npoly, i;
    double x, t, t1, t2, dt;
    double y_exact, y_approx;
    ChebEncoder encoder;
    double coeff[CHEB_MAX_POLYS];

    /* Process the command line parameters. */

    if (argc != 4)
    {
        fprintf(stderr, "USAGE:  demo npoly t1 t2\n");
        return 1;
    }

    if (1 != sscanf(argv[1], "%d", &npoly))
    {
        fprintf(stderr, "ERROR: Invalid number of polynomials on command line.\n");
        return 1;
    }

    if (1 != sscanf(argv[2], "%lf", &t1))
    {
        fprintf(stderr, "ERROR: Invalid value of t1 on command line.\n");
        return 1;
    }

    if (1 != sscanf(argv[3], "%lf", &t2))
    {
        fprintf(stderr, "ERROR: Invalid value of t2 on command line.\n");
        return 1;
    }

    /* Initialize the Chebyshev encoder. */

    if (ChebInit(&encoder, npoly))
    {
        fprintf(stderr, "ERROR: Invalid number of polynomials: %d\n", npoly);
        return 1;
    }

    /* Generate the coefficients that best approximate the function TestFunc. */

    ChebGenerate(&encoder, TestFunc, NULL, t1, t2, coeff);

    /* Output a table comparing the original function with its Chebyshev approximation. */
    dt = (t2 - t1) / (nsamples - 1);
    for (i=0; i < nsamples; ++i)
    {
        t = t1 + i*dt;
        y_exact = TestFunc(NULL, t);
        x = ChebScale(t1, t2, t);
        y_approx = ChebApprox(npoly, coeff, x);
        printf("t=%10.6lf, y_exact=%10.6lf, y_approx=%10.6lf\n", t, y_exact, y_approx);
    }

    return 0;
}
