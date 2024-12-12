#include <stdio.h>
#include <math.h>

#define EPSILON 1e-6 // Tolerance for convergence

// Define the function f(x)
double f(double x) {
    return pow(2, x) * pow((x - 2), 2) - 1;
}

// Define the derivative f'(x) for Newton's method
double df(double x) {
    double term1 = pow(2, x) * pow((x - 2), 2) * log(2); // Derivative of 2^x
    double term2 = 2 * (x - 2) * pow(2, x);             // Derivative of (x-2)^2
    return term1 + term2;
}

// Bisection method
double bisection(double a, double b) {
    double c;
    if (f(a) * f(b) >= 0) {
        printf("Bisection method fails: no root in the interval [%.6f, %.6f].\n", a, b);
        return NAN;
    }

    while ((b - a) >= EPSILON) {
        c = (a + b) / 2.0; // Midpoint
        if (fabs(f(c)) < EPSILON) { // Root found
            break;
        }
        else if (f(a) * f(c) < 0) {
            b = c;
        }
        else {
            a = c;
        }
    }
    return c;
}

// Newton's method
double newton(double x0) {
    double x1;
    int max_iterations = 1000;
    int iterations = 0;

    while (iterations < max_iterations) {
        double fx = f(x0);
        double dfx = df(x0);

        if (fabs(dfx) < EPSILON) {
            printf("Newton's method fails: derivative near zero.\n");
            return NAN;
        }

        x1 = x0 - fx / dfx;

        if (fabs(x1 - x0) < EPSILON) { // Converged
            return x1;
        }

        x0 = x1;
        iterations++;
    }

    printf("Newton's method failed to converge within %d iterations.\n", max_iterations);
    return NAN;
}

// Secant method
double secant(double x0, double x1) {
    double x2;
    int max_iterations = 1000;
    int iterations = 0;

    while (iterations < max_iterations) {
        double fx0 = f(x0);
        double fx1 = f(x1);

        if (fabs(fx1 - fx0) < EPSILON) {
            printf("Secant method fails: division by zero.\n");
            return NAN;
        }

        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        if (fabs(x2 - x1) < EPSILON) { // Converged
            return x2;
        }

        x0 = x1;
        x1 = x2;
        iterations++;
    }

    printf("Secant method failed to converge within %d iterations.\n", max_iterations);
    return NAN;
}

int main() {
    double a = 1.5, b = 3.0; // Interval for bisection
    double initial_guess = 2.5; // Initial guess for Newton's method
    double x0 = 1.5, x1 = 3.0; // Initial guesses for secant method

    // Bisection method
    double root_bisection = bisection(a, b);
    if (!isnan(root_bisection)) {
        printf("Root found by Bisection method: %.6f\n", root_bisection);
    }

    // Newton's method
    double root_newton = newton(initial_guess);
    if (!isnan(root_newton)) {
        printf("Root found by Newton's method: %.6f\n", root_newton);
    }

    // Secant method
    double root_secant = secant(x0, x1);
    if (!isnan(root_secant)) {
        printf("Root found by Secant method: %.6f\n", root_secant);
    }

    return 0;
}
