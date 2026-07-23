/***********************************************************************
 * Program: carvalho2025_table.cpp
 *
 * Detailed Description:
 * This program is designed to compute and tabulate normalized wave height parameters
 * for shallow-water environments, utilizing the Composite Weibull distribution model.
 * All wave heights are normalized by the root-mean-square wave height (Hrms) to provide
 * dimensionless ratios, which are fundamental for understanding wave statistics in complex
 * shallow foreshore conditions where wave breaking significantly alters traditional
 * wave height distributions (e.g., Rayleigh distribution).
 *
 * The core functionality of this program involves a robust numerical solution using
 * the Newton-Raphson method for a system of two non-linear equations:
 *
 * 1.  **Iterating through Normalized Transitional Wave Height (Htr/Hrms)**:
 * The program systematically varies the normalized transitional wave height (Htr_Hrms)
 * from 0.01 to 3.50 in increments of 0.01. The transitional wave height (Htr) is a
 * critical parameter in the Composite Weibull distribution, marking the point where
 * the wave height distribution transitions from a Rayleigh-like behavior (for smaller waves)
 * to a different, more complex behavior influenced by wave breaking (for larger waves).
 *
 * 2.  **Solving for Normalized Scale Parameters (H1/Hrms and H2/Hrms) Simultaneously**:
 * For each Htr_Hrms value, the program solves a system of two non-linear equations
 * to determine H1_Hrms and H2_Hrms simultaneously. This is a significant improvement
 * over sequential solving, ensuring consistency and accuracy in the relationship
 * between the two Weibull components.
 * The two equations are:
 * a) The normalized Hrms equation, which states that the overall normalized Hrms of the Composite Weibull distribution
 * must precisely equal one.
 * b) The continuity condition, which ensures
 * a smooth transition between the two Weibull components at Htr.
 * The Newton-Raphson method for systems of equations is employed, which involves
 * calculating the Jacobian matrix (matrix of partial derivatives) at each iteration.
 *
 * 3.  **Computing Normalized Quantile Wave Heights (H(1/N)/Hrms)**:
 * Once H1_Hrms and H2_Hrms are successfully determined, the program computes the
 * corresponding normalized quantile wave heights for a predefined set of exceedance
 * probabilities (represented by N values, e.g., N=3 for H1/3, N=10 for H1/10, etc.).
 * These calculations utilize specific formulas,
 * depending on the relationship between Htr_Hrms and the specific quantile being calculated.
 * These quantiles represent the mean of the highest 1/N-part of the wave heights,
 * providing crucial insights into the statistical properties of extreme waves.
 *
 * Output:
 * The program generates a formatted text file named "carvalho2025_table.txt". This file
 * contains a comprehensive table of Htr_Hrms, H1_Hrms, H2_Hrms, and the various
 * H(1/N)/Hrms values, providing a valuable reference for wave height analysis in
 * shallow-water coastal engineering applications.
 *
 * All gamma functions applied, either complete or incomplete, are non-normalized.
 *
 * Compilation Instructions:
 * To compile this program with high optimization levels, you can use a g++ command similar to this:
 *
 * g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 * -o carvalho2025_table.exe carvalho2025_table.cpp
 *
 ***********************************************************************/

// This preprocessor directive is crucial for ensuring that the M_PI constant (representing the mathematical constant Pi)
// is defined and available for use within the <cmath> header. This is particularly relevant when compiling
// on certain systems, such as Windows with Microsoft Visual C++ (MSVC), where it might not be defined by default.
#define _USE_MATH_DEFINES

// Standard C++ Library Includes:
// These lines import necessary functionalities from the C++ Standard Library, providing tools
// for input/output, file handling, mathematical operations, data structures, and formatting.
#include <iostream>     // Provides standard input/output streams (e.g., `std::cout` for console output, `std::cerr` for error messages).
#include <fstream>      // Enables file stream operations, allowing the program to read from and write to files (e.g., `std::ofstream` for output files).
#include <vector>       // Provides the `std::vector` container, a dynamic array that can resize itself.
#include <cmath>        // Contains a wide range of mathematical functions (e.g., `std::sqrt` for square root, `std::log` for natural logarithm, `std::pow` for exponentiation, `std::lgamma` for the natural logarithm of the complete gamma function, `std::tgamma` for the complete gamma function).
#include <limits>       // Provides `std::numeric_limits`, which allows querying properties of numeric types, suchs as the smallest representable long double (`std::numeric_limits<long double>::min()`), crucial for handling very small numbers and avoiding division by zero in numerical algorithms.
#include <iomanip>      // Offers manipulators for output formatting (e.g., `std::fixed` for fixed-point notation, `std::setprecision` to control decimal places, `std::setw` to set field width for formatted output).
#include <algorithm>    // Includes general-purpose algorithms (e.g., `std::max`, `std::swap`), useful for various operations.

// Using Namespace:
// This directive brings all identifiers (like `cout`, `endl`, `vector`, `sqrt`, etc.) from the
// `std` (standard) namespace directly into the current scope. This avoids the need to prefix
// these identifiers with `std::`, making the code more concise.
using namespace std;

// Global Parameters for the Composite Weibull Distribution:
// These constants define the fixed exponents for the two parts of the Composite Weibull distribution.
// They are declared globally as `const long double` because their values are fundamental to the model
// and remain unchanged throughout the program's execution.
const long double k1 = 2.0L;    // Exponent for the first part of the Composite Weibull distribution.
                          // As per the theoretical foundation (Groenendijk, 1998, Section 2.1 and 3.3.2),
                          // the initial part of the wave height distribution in shallow water is assumed
                          // to be Rayleigh-shaped. A Rayleigh distribution is a special case of the Weibull
                          // distribution with an exponent (shape parameter) of 2.0.
const long double k2 = 3.6L;    // Exponent for the second part of the Composite Weibull distribution.
                          // This value was empirically determined through calibration and optimization
                          // processes described in Groenendijk (1998, Section 2.1) and Groenendijk & Van Gent (1998).
                          // It reflects the altered shape of the wave height distribution for larger waves
                          // due to depth-induced breaking.

// Precision for Numerical Solver Convergence:
// This constant defines the tolerance level used in numerical methods (like Newton-Raphson)
// to determine when an iterative solution has converged to an acceptable accuracy.
const long double EPSILON = 1e-12L; // A small value (10^-12) indicating the maximum allowable error or difference
                             // between successive iterations for a solution to be considered converged.
const long double JACOBIAN_DX = 1e-8L; // Small step size for finite difference approximation of Jacobian derivatives.

// Forward Declarations of Functions:
// These declarations inform the C++ compiler about the existence and signature (return type, name, and parameters)
// of functions that are defined later in the source code. This is necessary because some functions might call
// other functions that have not yet been fully defined at the point of the call. This prevents compilation errors.

// Declares a function to compute the unnormalized lower incomplete gamma function γ(a, x) = ∫[0,x] t^(a-1)e^(-t) dt.
long double incomplete_gamma_lower(long double a, long double x);
// Declares a function to compute the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
long double incomplete_gamma_upper(long double a, long double x);
// Declares a function to calculate HN (the wave height with a 1/N exceedance probability).
long double calculate_HN(long double N, long double H1, long double H2, long double k1, long double k2, long double Htr);
// Declares a function to calculate H1/N (the mean of the highest 1/N-th part of wave heights).
long double calculate_H1N(long double N_val, long double H1, long double H2, long double k1, long double k2, long double Htr);

// Functions for the system of non-linear equations
// F1 represents the normalized Hrms equation.
long double F1(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms);
// F2 represents the continuity condition.
long double F2(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms);

// Function to solve the 2x2 linear system Ax = b using Cramer's rule (or simple substitution)
// for the Newton-Raphson update.
void solve_linear_system_2x2(long double J11, long double J12, long double J21, long double J22,
                             long double b1, long double b2, long double &dx1, long double &dx2);

// Declares the Newton-Raphson solver for a system of two equations.
bool newtonRaphsonSystemSolver(long double Htr_Hrms, long double &H1_Hrms, long double &H2_Hrms,
                               long double tol, int maxit);
// Declares a function to provide initial guesses for the numerical solver.
void get_initial_guesses(long double Htr_Hrms, long double &H1_initial, long double &H2_initial);


/**
 * @brief Computes the unnormalized lower incomplete gamma function γ(a, x) = ∫[0,x] t^(a-1)e^(-t) dt.
 *
 * This function is a critical component for calculating moments of Weibull distributions,
 * which are fundamental to the Composite Weibull model. It employs a hybrid approach
 * for numerical stability and accuracy:
 * - For small values of 'x' (specifically, x < a + 1.0), it uses a series expansion.
 * - For larger values of 'x', it utilizes a continued fraction expansion.
 * This adaptive strategy ensures robust and precise computation across different input ranges.
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (upper integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of γ(a, x). Returns `nan("")` (Not-a-Number) if input parameters are invalid
 * or if the series/continued fraction fails to converge within the maximum iterations.
 */
long double incomplete_gamma_lower(long double a, long double x)
{
    const int MAXIT = 500;      // Maximum number of iterations allowed for either the series expansion or the continued fraction.
    const long double LOCAL_EPS = 1e-16L;   // A very small tolerance value (10^-16) used for checking convergence
                                      // within the gamma function calculations, ensuring high numerical precision.

    // Input validation: Essential for preventing mathematical errors and ensuring function robustness.
    if(a <= 0.0L || x < 0.0L)
        return nan("0"); // If 'a' is not positive or 'x' is negative, return NaN as the result is undefined.
    if(x == 0.0L)
        return 0.0L; // By definition, γ(a, 0) is 0 for any positive 'a'.

    // Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
    // `std::lgamma` is used instead of `std::tgamma` followed by `log` for improved numerical stability,
    // especially when 'a' is large, as `lgamma` avoids potential overflow/underflow issues with very large/small gamma values.
    long double gln = lgamma(a);

    // Conditional logic to select the appropriate numerical method (series or continued fraction).
    // This choice is based on the relationship between 'x' and 'a', a common heuristic for optimal performance and accuracy.
    if(x < a + 1.0L) {  // Use series expansion for x < a+1.
        long double ap = a;          // `ap` is a mutable copy of 'a', used to increment in the series terms.
        long double sum = 1.0L / a;   // Initialize the sum with the first term of the series expansion.
        long double del = sum;       // `del` stores the value of the current term being added to the sum.
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) { // Iterate up to MAXIT to compute series terms.
            ap += 1.0L;          // Increment 'a' for the next term in the series (a+1, a+2, ...).
            del *= x / ap;      // Calculate the next term: `del_n = del_{n-1} * (x / (a + n))`.
            sum += del;         // Add the newly calculated term to the cumulative sum.
            // Check for convergence: If the absolute value of the current term (`del`) is extremely small
            // relative to the absolute value of the cumulative sum (`sum`), then the series has converged.
            if(abs(del) < abs(sum) * LOCAL_EPS)
                // Return the final result for γ(a,x).
                return sum * exp(-x + a * log(x) - gln) * tgamma(a);
        }
        return nan("0"); // If the loop finishes without `del` becoming sufficiently small, it means convergence failed.
    } else {  // Use continued fraction for x >= a+1.
        long double b = x + 1.0L - a; // Initialize the 'b' term (B_0) for the continued fraction expansion.
        // Initialize 'c' and 'd' to extremely large/small values. This is a common technique
        // to prevent division by zero or underflow during the initial steps of the continued fraction
        // algorithm, where denominators might otherwise be zero or very close to zero.
        long double c = 1.0L / numeric_limits<long double>::min();
        long double d = 1.0L / b;
        long double h = d; // `h` accumulates the value of the continued fraction.
        for (int i = 1; i <= MAXIT; ++i) { // Iterate up to MAXIT to compute continued fraction terms.
            long double an = -1.0L * i * (i - a); // Calculate the 'a_n' term (numerator) for the current iteration.
            b += 2.0L; // Update the 'b_n' term (denominator) for the current iteration.
            d = an * d + b; // Apply the recurrence relation for D_n (denominator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'd'.
            if(abs(d) < numeric_limits<long double>::min())
                d = numeric_limits<long double>::min();
            c = b + an / c; // Apply the recurrence relation for C_n (numerator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'c'.
            if(abs(c) < numeric_limits<long double>::min())
                c = numeric_limits<long double>::min();
            d = 1.0L / d;    // Invert 'd' for the next step.
            long double del = d * c; // Calculate the current term (delta_n) of the continued fraction.
            h *= del;       // Multiply `h` by the current term to update the accumulated fraction value.
            // Check for convergence: If the current term `del` is very close to 1.0, the continued fraction has converged.
            if(abs(del - 1.0L) < LOCAL_EPS)
                break; // Exit the loop if converged.
        }
        // This part of the calculation results in Q(a,x) before multiplication by Gamma(a).
        long double lnPart = -x + a * log(x) - gln;
        long double Qval_normalized = exp(lnPart) * h;
        // P(a,x) = 1 - Q(a,x).
        long double Pval_normalized = 1.0L - Qval_normalized;
        // It's clamped to 0.0 to prevent very small negative floating-point results due to precision issues.
        Pval_normalized = (Pval_normalized < 0.0L) ? 0.0L : Pval_normalized;
        // Return unnormalized lower incomplete gamma γ(a,x) = P(a,x) * Γ(a).
        return Pval_normalized * tgamma(a);
    }
}


/**
 * @brief Computes the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
 *
 * This function calculates the unnormalized form of the upper incomplete gamma function.
 * It is derived from the complete gamma function Γ(a) and the unnormalized lower incomplete gamma function γ(a, x),
 * using the identity Γ(a,x) = Γ(a) - γ(a,x).
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
 * @throws `std::invalid_argument` if input parameters ('a' or 'x') are outside their valid ranges.
 */
long double incomplete_gamma_upper(long double a, long double x) {
    // Input validation: Ensures 'a' is positive and 'x' is non-negative, as required by the gamma function definitions.
    if (a <= 0.0L) {
        throw invalid_argument("incomplete_gamma_upper: 'a' must be positive.");
    }
    if (x < 0.0L) {
        throw invalid_argument("incomplete_gamma_upper: 'x' must be non-negative.");
    }
    // Calculation: Γ(a,x) = Γ(a) - γ(a,x).
    // `std::tgamma` computes the complete gamma function Γ(a).
    return tgamma(a) - incomplete_gamma_lower(a, x);
}


/**
 * @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
 *
 * This function determines the specific wave height (H) such that the probability of a wave
 * exceeding this height is 1/N. This is a key statistical measure. The calculation depends
 * on whether the target wave height falls into the first or second part of the Composite Weibull
 * distribution, which is separated by the transitional wave height (Htr).
 * The logic directly implements the formulas from Appendix A.2.1 of Groenendijk (1998).
 *
 * @param N The N value (e.g., 3 for H1/3, 100 for H1%). Must be strictly greater than 1.0
 * because `log(N)` is used, and `N=1` would result in `log(1)=0`.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent (shape parameter) of the first Weibull distribution. Must be positive.
 * @param k2 The exponent (shape parameter) of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height, which defines the boundary between the two Weibull parts.
 * @return The calculated value of HN.
 * @throws `std::invalid_argument` if input parameters are invalid (e.g., N <= 1, H1/H2 <= 0, k1/k2 <= 0).
 */
long double calculate_HN(long double N, long double H1, long double H2, long double k1, long double k2, long double Htr) {
    // Input validation: Ensures all parameters are within their mathematically valid ranges.
    if (N <= 1.0L) {
        throw invalid_argument("calculate_HN: N must be greater than 1 for log(N).");
    }
    if (H1 <= 0.0L || H2 <= 0.0L) {
        throw invalid_argument("calculate_HN: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0L || k2 <= 0.0L) {
        throw invalid_argument("calculate_HN: k1 and k2 must be positive.");
    }

    // Calculate a candidate HN assuming it falls within the *first* part of the distribution (H <= Htr).
    // This uses Equation A.5 from Groenendijk (1998): `H_N,1 = H1 * (ln(N))^(1/k1)`.
    long double HN_candidate1 = H1 * pow(log(N), 1.0L / k1);

    // Decision point: Check if the `HN_candidate1` is indeed less than the transitional wave height (Htr).
    // A small `EPSILON` is subtracted from `Htr` to account for floating-point inaccuracies when comparing.
    if (HN_candidate1 < Htr - EPSILON) {
        // If true, it means the wave height for the given exceedance probability is governed by the first Weibull part.
        return HN_candidate1;
    } else {
        // If `HN_candidate1` is not less than `Htr`, it implies that the wave height for the given
        // exceedance probability is determined by the *second* part of the distribution (H > Htr).
        // This uses Equation A.7 from Groenendijk (1998): `H_N = H2 * (ln(N))^(1/k2)`.
        return H2 * pow(log(N), 1.0L / k2);
    }
}

/**
 * @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N) for the Composite Weibull distribution.
 *
 * This function computes a characteristic wave height that represents the average height of the
 * highest N-th fraction of waves in a given wave field. This is a commonly used metric in
 * oceanography and coastal engineering (e.g., H1/3 for significant wave height).
 * The calculation involves integrals of the probability density function and depends on
 * whether the relevant wave heights fall into the first or second part of the Composite Weibull distribution.
 * The implementation follows the detailed derivations in Appendix A.2.2 of Groenendijk (1998).
 *
 * @param N_val The N parameter for H1/N (e.g., 3 for H1/3, 10 for H1/10). Must be strictly greater than 1.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent (shape parameter) of the first Weibull distribution. Must be positive.
 * @param k2 The exponent (shape parameter) of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height.
 * @return The calculated value of H1/N.
 * @throws `std::invalid_argument` if input parameters are invalid (e.g., H1/H2 <= 0, k1/k2 <= 0, N_val <= 1).
 */
long double calculate_H1N(long double N_val, long double H1, long double H2, long double k1, long double k2, long double Htr) {
    // Input validation: Ensures all parameters are within their mathematically valid ranges.
    if (H1 <= 0.0L || H2 <= 0.0L) {
        throw invalid_argument("calculate_H1N: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0L || k2 <= 0.0L) {
        throw invalid_argument("calculate_H1N: k1 and k2 must be positive.");
    }
    if (N_val <= 1.0L) {
        throw invalid_argument("calculate_H1N: N_val must be greater than 1.");
    }

    // First, determine HN (the wave height with 1/N exceedance probability).
    // This is crucial because the integration limits and the specific formula for H1/N
    // depend on where HN falls relative to Htr.
    long double H_N_val = calculate_HN(N_val, H1, H2, k1, k2, Htr);

    // Case 1: H_N_val is smaller than Htr.
    // This implies that the integration for H1/N spans both parts of the Composite Weibull distribution.
    // The formula used here corresponds to Equation A.15 in Groenendijk (1998).
    if (H_N_val < Htr - EPSILON) {
        // Term 1: Contribution from the first Weibull distribution (F1(H)).
        // The 'a' parameter for the incomplete gamma function in this context is (1/k1 + 1).
        long double term1_a = 1.0L / k1 + 1.0L;
        // The 'x' parameter for the first incomplete gamma term is `ln(N_val)`.
        long double term1_x_ln_Nval = log(N_val);
        // The 'x' parameter for the second incomplete gamma term (related to Htr) is `(Htr/H1)^k1`.
        long double term1_x_HtrH1 = pow(Htr / H1, k1);

        // Calculate the two unnormalized upper incomplete gamma terms for the first part of the integral.
        // This represents `Γ[a, ln(N)] - Γ[a, (Htr/H1)^k1]`.
        long double gamma_term1_part1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval);
        long double gamma_term1_part2 = incomplete_gamma_upper(term1_a, term1_x_HtrH1);
        long double gamma_term1 = gamma_term1_part1 - gamma_term1_part2;

        // Term 2: Contribution from the second Weibull distribution (F2(H)).
        // The 'a' parameter for the incomplete gamma function is (1/k2 + 1).
        long double term2_a = 1.0L / k2 + 1.0L;
        // The 'x' parameter for this term is `(Htr/H2)^k2`.
        long double term2_x_HtrH2 = pow(Htr / H2, k2);

        // Calculate the unnormalized upper incomplete gamma term for the second part of the integral.
        // This represents `Γ[a, (Htr/H2)^k2]`.
        long double gamma_term2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2);

        // Combine the terms as per Equation A.15: `N_val * H1 * (gamma_term1) + N_val * H2 * (gamma_term2)`.
        return N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2;
    } else {
        // Case 2: H_N_val is greater than or equal to Htr.
        // This means the integration for H1/N only involves the second part of the Composite Weibull distribution.
        // The formula used here corresponds to Equation A.20 in Groenendijk (1998).
        long double term_a = 1.0L / k2 + 1.0L; // The 'a' parameter for the incomplete gamma function.
        long double term_x = log(N_val);     // The 'x' parameter is `ln(N_val)`.
        // Calculate `N_val * H2 * Γ[a, ln(N_val)]`.
        return N_val * H2 * incomplete_gamma_upper(term_a, term_x);
    }
}

// --- Functions for the System of Non-Linear Equations (Newton-Raphson Matrix Method) ---

/**
 * @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
 * It is derived from Equation 7.11 from Groenendijk (1998) and Equation 9 in Caires & Van Gent (2012),
 * and states that the overall normalized Hrms of the Composite Weibull distribution must precisely equal 1.
 * The equation directly uses the unnormalized lower incomplete gamma function, γ(a,x),
 * and the unnormalized upper incomplete gamma function, Γ(a,x).
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the first function, which should be driven to zero.
 */
long double F1(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms) {
    // Input validation for H1_Hrms and H2_Hrms to prevent issues like log(0) or sqrt(negative).
    // While the solver should ideally keep values positive, these checks add robustness.
    if (H1_Hrms <= 0.0L || H2_Hrms <= 0.0L) {
        // Return a large value to push the solver away from invalid regions.
        return numeric_limits<long double>::max();
    }

    long double arg1 = pow(Htr_Hrms / H1_Hrms, k1);
    long double arg2 = pow(Htr_Hrms / H2_Hrms, k2);

    // Calculate terms using unnormalized incomplete gamma functions as per Equation 9.
    // This directly corresponds to the terms H1^2 * γ(2/k1 + 1, (Htr/H1)^k1)
    // and H2^2 * Γ(2/k2 + 1, (Htr/H2)^k2) in Equation 9.
    long double term1 = H1_Hrms * H1_Hrms * incomplete_gamma_lower(2.0L / k1 + 1.0L, arg1);
    long double term2 = H2_Hrms * H2_Hrms * incomplete_gamma_upper(2.0L / k2 + 1.0L, arg2);

    // Ensure the argument to sqrt is non-negative, though theoretically it should be.
    long double sum_terms = term1 + term2;
    if (sum_terms < 0.0L) sum_terms = 0.0L; // Clamp to zero to avoid NaN from sqrt of negative.

    return sqrt(sum_terms) - 1.0L;
}

/**
 * @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the continuity condition between the two Weibull distributions
 * at the transitional wave height Htr. It is derived from Equation 3.4 in Groenendijk (1998)
 * and Equation 8 in Caires & Van Gent (2012):
 * `(Htr/H1)^k1 = (Htr/H2)^k2`. Rearranging this gives `(Htr/H1)^k1 - (Htr/H2)^k2 = 0`.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the second function, which should be driven to zero.
 */
long double F2(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms) {
    // Input validation for H1_Hrms and H2_Hrms to prevent division by zero or log(0).
    if (H1_Hrms <= 0.0L || H2_Hrms <= 0.0L) {
        // Return a large value to push the solver away from invalid regions.
        return numeric_limits<long double>::max();
    }
    return pow(Htr_Hrms / H1_Hrms, k1) - pow(Htr_Hrms / H2_Hrms, k2);
}

/**
 * @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
 *
 * This function is a helper for the Newton-Raphson method for systems.
 * It takes the Jacobian matrix elements (J11, J12, J21, J22) and the negative
 * function values (-F1, -F2) as the right-hand side, and computes the updates
 * (dx1, dx2) for H1_Hrms and H2_Hrms.
 *
 * @param J11 Element (1,1) of the Jacobian matrix (dF1/dH1).
 * @param J12 Element (1,2) of the Jacobian matrix (dF1/dH2).
 * @param J21 Element (2,1) of the Jacobian matrix (dF2/dH1).
 * @param J22 Element (2,2) of the Jacobian matrix (dF2/dH2).
 * @param b1 Right-hand side for the first equation (-F1).
 * @param b2 Right-hand side for the second equation (-F2).
 * @param dx1 Output: The calculated change for H1_Hrms.
 * @param dx2 Output: The calculated change for H2_Hrms.
 */
void solve_linear_system_2x2(long double J11, long double J12, long double J21, long double J22,
                             long double b1, long double b2, long double &dx1, long double &dx2) {
    long double determinant = J11 * J22 - J12 * J21;

    // Check for a singular or nearly singular Jacobian matrix.
    if (abs(determinant) < numeric_limits<long double>::epsilon() * 100L) {
        // If determinant is too small, the matrix is singular or ill-conditioned.
        // In this case, we cannot reliably solve the system.
        // Set dx1 and dx2 to zero to prevent large, unstable steps.
        // A more robust solver might use a pseudo-inverse or a different fallback.
        dx1 = 0.0L;
        dx2 = 0.0L;
        return;
    }

    // Apply Cramer's rule:
    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
}

/**
 * @brief Provides initial guesses for H1_Hrms and H2_Hrms based on Htr_Hrms.
 *
 * A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
 * This function uses an empirical regression for H1_Hrms, and then derives H2_Hrms from the
 * continuity condition, assuming an initial relationship.
 * The regression for H1_Hrms is based on Groenendijk's findings.
 *
 * @param Htr_Hrms The normalized transitional wave height.
 * @param H1_initial Output: The initial guess for H1_Hrms.
 * @param H2_initial Output: The initial guess for H2_Hrms.
 */
void get_initial_guesses(long double Htr_Hrms, long double &H1_initial, long double &H2_initial) {
    // Empirical regression for H1/Hrms.
    H1_initial = 2.244660800090239E-03 + std::pow(std::tanh(1.918610494219390E+00 * Htr_Hrms), 1.780892753373355E-01) / std::pow(std::tanh(std::sinh(1.009497360864962E+00 * Htr_Hrms)), 9.777939607559606E-01);

    // Empirical regression for H2_initial
    H2_initial = 1.059259665431797 + (0.2059286860468916 * Htr_Hrms) / (1.0 + 3.865701948059343 * std::pow(Htr_Hrms, -3.479682433107255));

    // Ensure initial guesses are positive, as physical parameters cannot be zero or negative.
    if (H1_initial <= 0.0L) H1_initial = numeric_limits<long double>::min(); // Small positive value
    if (H2_initial <= 0.0L) H2_initial = numeric_limits<long double>::min(); // Small positive value
}

/**
 * @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
 *
 * This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
 * of the system of non-linear equations F1 and F2. It iteratively refines the guesses
 * for H1_Hrms and H2_Hrms until the functions F1 and F2 are sufficiently close to zero.
 *
 * @param Htr_Hrms The normalized transitional wave height (constant for this solve).
 * @param H1_Hrms Output: The converged normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms Output: The converged normalized scale parameter of the second Weibull distribution.
 * @param tol The desired tolerance for convergence (maximum absolute value of F1 and F2).
 * @param maxit The maximum number of iterations allowed.
 * @return `true` if the solver successfully converges, `false` otherwise.
 */
bool newtonRaphsonSystemSolver(long double Htr_Hrms, long double &H1_Hrms, long double &H2_Hrms,
                               long double tol, int maxit) {
    // Get initial guesses for H1_Hrms and H2_Hrms.
    get_initial_guesses(Htr_Hrms, H1_Hrms, H2_Hrms);

    for (int iter = 0; iter < maxit; ++iter) {
        // Evaluate the functions at the current guesses.
        long double f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        long double f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        // Check for convergence. If both function values are close to zero, we've converged.
        if (abs(f1_val) < tol && abs(f2_val) < tol) {
            return true;
        }

        // Calculate the Jacobian matrix elements using central finite differences.
        // J11 = dF1/dH1
        long double J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        // J12 = dF1/dH2
        long double J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        // J21 = dF2/dH1
        long double J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        // J22 = dF2/dH2
        long double J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0L * JACOBIAN_DX);

        // Solve the linear system J * dx = -F for dx.
        // Here, dx = [dH1, dH2]^T and F = [f1_val, f2_val]^T.
        long double dH1, dH2;
        solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2);

        // Update the guesses.
        H1_Hrms += dH1;
        H2_Hrms += dH2;

        // Ensure H1_Hrms and H2_Hrms remain positive. If they become non-positive,
        // clamp them to a small positive value to prevent mathematical errors in subsequent iterations.
        if (H1_Hrms <= 0.0L) H1_Hrms = numeric_limits<long double>::min();
        if (H2_Hrms <= 0.0L) H2_Hrms = numeric_limits<long double>::min();
    }

    // If the loop finishes without converging, print a warning.
    cerr << "Newton-Raphson system solver failed to converge for Htr_Hrms = " << Htr_Hrms << " after " << maxit << " iterations.\n";
    return false; // Indicate failure to converge.
}


/**
 * @brief Main function of the program.
 *
 * This is the entry point of the C++ program. It orchestrates the entire computation
 * and output process. The function iterates through a predefined range of normalized
 * transitional wave heights (Htr_Hrms), for each value performing the necessary
 * numerical solves and calculations to determine the corresponding normalized
 * wave parameters of the Composite Weibull distribution. Finally, it writes all
 * computed results to a formatted text file.
 *
 * @return 0 if the program executes successfully, 1 if there is an error (e.g., file opening failure).
 */
int main()
{
    // Open an output file named "carvalho2025_table.txt" for writing.
    // `std::ofstream` is used to create an output file stream.
    ofstream fout("carvalho2025_table.txt");
    // Check if the file was successfully opened. If `fout` is in a bad state (e.g., due to permissions or disk full),
    // it evaluates to `false`.
    if (!fout) {
        cerr << "Error opening output file: carvalho2025_table.txt" << endl; // Print an error message to console.
        return 1; // Return a non-zero exit code to indicate an error.
    }

    // Define the range and step size for the normalized transitional wave height (Htr/Hrms).
    // These constants control the granularity and extent of the generated table.
    const long double Htr_Hrms_start = 0.01L; // The starting value for Htr/Hrms.
    const long double Htr_Hrms_end = 3.50L;   // The ending value for Htr/Hrms.
    const long double Htr_Hrms_step = 0.01L;  // The increment step size for Htr/Hrms.

    // Define an array of integer N values for which the H(1/N)/Hrms quantiles will be calculated.
    // These N values correspond to specific exceedance probabilities or mean-of-highest-N-part values.
    // For example, N=3 corresponds to H1/3 (significant wave height).
    const int N_values[] = {3, 10, 50, 100, 250, 1000};

    // Write the header row to the output file.
    // `std::setw(X)` is a manipulator that sets the field width for the *next* output item to X characters.
    // This ensures that the columns in the output file are neatly aligned.
    fout << setw(9) << "Htr/Hrms"
         << setw(9) << "H1/Hrms"
         << setw(9) << "H2/Hrms"
         << setw(11) << "H1/3/Hrms"    // Normalized significant wave height.
         << setw(11) << "H1/10/Hrms"   // Normalized mean of the highest 1/10th waves.
         << setw(11) << "H1/50/Hrms"   // Normalized mean of the highest 1/50th waves.
         << setw(12) << "H1/100/Hrms"  // Normalized mean of the highest 1/100th waves (often H1%).
         << setw(12) << "H1/250/Hrms"  // Normalized mean of the highest 1/250th waves.
         << setw(13) << "H1/1000/Hrms" // Normalized mean of the highest 1/1000th waves.
         << endl; // Insert a newline character to move to the next line after the header.

    // Main loop: Iterate through the defined range of Htr_Hrms values.
    // The `+ EPSILON` in the loop condition helps to ensure that the `Htr_Hrms_end` value
    // is included in the iteration, compensating for potential floating-point arithmetic inaccuracies.
    for (long double Htr_Hrms = Htr_Hrms_start; Htr_Hrms <= Htr_Hrms_end + EPSILON; Htr_Hrms += Htr_Hrms_step) {
        try { // Begin a try-catch block to gracefully handle potential runtime errors during calculations.
            long double H1_normalized; // Declare a variable to store the computed normalized H1 (H1/Hrms).
            long double H2_normalized; // Declare a variable to store the computed normalized H2 (H2/Hrms).

            // Call the `newtonRaphsonSystemSolver` function to find H1_normalized and H2_normalized simultaneously.
            if (!newtonRaphsonSystemSolver(Htr_Hrms, H1_normalized, H2_normalized, EPSILON, 100)) {
                // If the solver returns `false`, it indicates a failure to converge, so throw a runtime_error.
                throw runtime_error("Newton-Raphson system solver failed to find H1_Hrms and H2_Hrms.");
            }

            // Create a dynamic array (vector) to store the calculated H1/N values for the current Htr_Hrms.
            vector<long double> h1n_values;
            // Iterate through each predefined N value.
            for (int N_val : N_values) {
                // Calculate the H1/N value using the `calculate_H1N` function and add it to the vector.
                // `static_cast<long double>(N_val)` ensures that N_val is treated as a long double in the function call.
                h1n_values.push_back(calculate_H1N(static_cast<long double>(N_val), H1_normalized, H2_normalized, k1, k2, Htr_Hrms));
            }

            // Set output formatting for the current row in the file.
            // `std::fixed` ensures floating-point numbers are printed in fixed-point notation (not scientific).
            // `std::setprecision(5)` sets the number of digits after the decimal point to 5.
            fout << fixed << setprecision(5);
            // Write the primary normalized parameters (Htr/Hrms, H1/Hrms, H2/Hrms) to the file,
            // formatted with specified widths.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << H1_normalized
                 << setw(9) << H2_normalized;

            // Output each calculated H1/N value from the vector to the file,
            // also formatted with specified widths.
            fout << setw(11) << h1n_values[0]  // H1/3/Hrms
                 << setw(11) << h1n_values[1]  // H1/10/Hrms
                 << setw(11) << h1n_values[2]  // H1/50/Hrms
                 << setw(12) << h1n_values[3]  // H1/100/Hrms
                 << setw(12) << h1n_values[4]  // H1/250/Hrms
                 << setw(13) << h1n_values[5]  // H1/1000/Hrms
                 << endl; // Move to the next line in the output file.

        } catch (const invalid_argument& e) { // Catch `std::invalid_argument` exceptions.
            // These typically occur due to invalid input parameters passed to functions (e.g., log of non-positive number).
            cerr << "Input error for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(13) << "ERROR"
                 << endl;
        } catch (const runtime_error& e) { // Catch `std::runtime_error` exceptions.
            // These typically indicate issues during the execution, such as solver failure to converge.
            cerr << "Calculation error for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(13) << "ERROR"
                 << endl;
        } catch (const exception& e) { // Catch any other standard C++ exceptions.
            // This is a general catch-all for unexpected errors.
            cerr << "An unexpected error occurred for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(13) << "ERROR"
                 << endl;
        }
    }

    fout.close(); // Close the output file stream, ensuring all buffered data is written to the file.
    cout << "Calculation complete. Results saved to carvalho2025_table.txt" << endl; // Inform the user via console.

    return 0; // Return 0 to indicate that the program executed successfully without critical errors.
}
