/***********************************************************************
 * Program: iterations.cpp
 *
 * Detailed Description:
 * This program calculates the number of iterations required for a Newton-Raphson
 * solver to find the parameters (H1/Hrms, H2/Hrms) of a Composite Weibull
 * distribution for a range of dimensionless transitional wave heights (Htr/Hrms).
 *
 * The program iterates through Htr/Hrms values from 0.01 to 3.5 (step 0.01).
 * For each value, it solves a system of two non-linear equations using a
 * Newton-Raphson method, with the core calculation logic being identical to
 * the one used in shallow-water-waves_cli.cpp.
 *
 * The detailed output of each calculation, including the initial guess and
 * step-by-step iteration values, is written to "output.txt".
 *
 * At the end of the execution, the program prints a summary to the console,
 * and also appends it to "output.txt", showing the minimum, average, and
 * maximum number of iterations required for convergence across the entire
 * tested range of Htr/Hrms values.
 *
 * Compilation Instructions (example using g++ on Windows):
 * g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 * -o iterations iterations.cpp
 ***********************************************************************/

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <numeric> // For std::accumulate

using namespace std;

// Global Parameters for the Composite Weibull Distribution:
const long double k1 = 2.0L;
const long double k2 = 3.6L;

// Precision for Numerical Solver Convergence:
const long double EPSILON = 1e-12L;
const long double JACOBIAN_DX = 1e-8L;

// Forward Declarations of Functions:
long double incomplete_gamma_lower(long double a, long double x);
long double incomplete_gamma_upper(long double a, long double x);
long double F1(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms);
long double F2(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms);
void solve_linear_system_2x2(long double J11, long double J12, long double J21, long double J22,
                             long double b1, long double b2, long double &dx1, long double &dx2);
void get_initial_guesses(long double Htr_Hrms, long double &H1_initial, long double &H2_initial);
int newtonRaphsonSystemSolver(long double Htr_Hrms, long double &H1_Hrms, long double &H2_Hrms,
                              int maxit, ofstream &outputFile);


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
    const int MAXIT = 500;
    const long double LOCAL_EPS = 1e-16L;

    if(a <= 0.0L || x < 0.0L)
        return nan("0");
    if(x == 0.0L)
        return 0.0L;

    long double gln = lgamma(a);

    if(x < a + 1.0L) {
        long double ap = a;
        long double sum = 1.0L / a;
        long double del = sum;
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) {
            ap += 1.0L;
            del *= x / ap;
            sum += del;
            if(abs(del) < abs(sum) * LOCAL_EPS)
                return sum * exp(-x + a * log(x) - gln) * tgamma(a);
        }
        return nan("0");
    } else {
        long double b = x + 1.0L - a;
        long double c = 1.0L / numeric_limits<long double>::min();
        long double d = 1.0L / b;
        long double h = d;
        for (int i = 1; i <= MAXIT; ++i) {
            long double an = -1.0L * i * (i - a);
            b += 2.0L;
            d = an * d + b;
            if(abs(d) < numeric_limits<long double>::min())
                d = numeric_limits<long double>::min();
            c = b + an / c;
            if(abs(c) < numeric_limits<long double>::min())
                c = numeric_limits<long double>::min();
            d = 1.0L / d;
            long double del = d * c;
            h *= del;
            if(abs(del - 1.0L) < LOCAL_EPS)
                break;
        }
        long double lnPart = -x + a * log(x) - gln;
        long double Qval_normalized = exp(lnPart) * h;
        long double Pval_normalized = 1.0L - Qval_normalized;
        Pval_normalized = (Pval_normalized < 0.0L) ? 0.0L : Pval_normalized;
        return Pval_normalized * tgamma(a);
    }
}


/**
 * @brief Computes the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
 *
 * It is derived from the complete gamma function Γ(a) and the unnormalized lower incomplete gamma function γ(a, x),
 * using the identity Γ(a,x) = Γ(a) - γ(a,x).
 *
 * @param a The 'a' parameter (shape parameter). Must be positive.
 * @param x The 'x' parameter (lower integration limit). Must be non-negative.
 * @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
 * @throws `std::invalid_argument` if input parameters are outside their valid ranges.
 */
long double incomplete_gamma_upper(long double a, long double x) {
    if (a <= 0.0L) {
        throw invalid_argument("incomplete_gamma_upper: 'a' must be positive.");
    }
    if (x < 0.0L) {
        throw invalid_argument("incomplete_gamma_upper: 'x' must be non-negative.");
    }
    return tgamma(a) - incomplete_gamma_lower(a, x);
}


/**
 * @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the first function, which should be driven to zero.
 */
long double F1(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms) {
    if (H1_Hrms <= 0.0L || H2_Hrms <= 0.0L) {
        return numeric_limits<long double>::max();
    }

    long double arg1 = pow(Htr_Hrms / H1_Hrms, k1);
    long double arg2 = pow(Htr_Hrms / H2_Hrms, k2);

    long double term1 = H1_Hrms * H1_Hrms * incomplete_gamma_lower(2.0L / k1 + 1.0L, arg1);
    long double term2 = H2_Hrms * H2_Hrms * incomplete_gamma_upper(2.0L / k2 + 1.0L, arg2);

    long double sum_terms = term1 + term2;
    if (sum_terms < 0.0L) sum_terms = 0.0L;

    return sqrt(sum_terms) - 1.0L;
}

/**
 * @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the continuity condition between the two Weibull distributions
 * at the transitional wave height Htr.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the second function, which should be driven to zero.
 */
long double F2(long double H1_Hrms, long double H2_Hrms, long double Htr_Hrms) {
    if (H1_Hrms <= 0.0L || H2_Hrms <= 0.0L) {
        return numeric_limits<long double>::max();
    }
    return pow(Htr_Hrms / H1_Hrms, k1) - pow(Htr_Hrms / H2_Hrms, k2);
}

/**
 * @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
 *
 * This is a helper function for the Newton-Raphson method.
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

    if (abs(determinant) < numeric_limits<long double>::epsilon() * 100L) {
        dx1 = 0.0L;
        dx2 = 0.0L;
        return;
    }

    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
}

/**
 * @brief Provides empirical initial guesses for H1_Hrms and H2_Hrms based on Htr_Hrms.
 *
 * A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
 *
 * @param Htr_Hrms The normalized transitional wave height.
 * @param H1_initial Output: The initial guess for H1_Hrms.
 * @param H2_initial Output: The initial guess for H2_Hrms.
 */
void get_initial_guesses(long double Htr_Hrms, long double &H1_initial, long double &H2_initial) {
    H1_initial = 2.244660800090239E-03 + std::pow(std::tanh(1.918610494219390E+00 * Htr_Hrms), 1.780892753373355E-01) / std::pow(std::tanh(std::sinh(1.009497360864962E+00 * Htr_Hrms)), 9.777939607559606E-01);
    H2_initial = 1.059259665431797 + (0.2059286860468916 * Htr_Hrms) / (1.0 + 3.865701948059343 * std::pow(Htr_Hrms, -3.479682433107255));

    if (H1_initial <= 0.0L) H1_initial = numeric_limits<long double>::min();
    if (H2_initial <= 0.0L) H2_initial = numeric_limits<long double>::min();
}

/**
 * @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
 *
 * This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
 * of the system of non-linear equations F1 and F2. The core calculation logic is
 * identical to that in shallow-water-waves_cli.cpp.
 *
 * @param Htr_Hrms The normalized transitional wave height (constant for this solve).
 * @param H1_Hrms Input/Output: The initial guess and converged normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms Input/Output: The initial guess and converged normalized scale parameter of the second Weibull distribution.
 * @param maxit The maximum number of iterations allowed.
 * @param outputFile The output file stream to write iteration details to.
 * @return The number of iterations taken to converge, or -1 if the solver fails to converge.
 */
int newtonRaphsonSystemSolver(long double Htr_Hrms, long double &H1_Hrms, long double &H2_Hrms,
                               int maxit, ofstream &outputFile) {
    // The loop now runs up to and including maxit, to allow for the initial guess (iter=0)
    // plus 'maxit' actual iterations.
    for (int iter = 0; iter <= maxit; ++iter) {
        long double f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        long double f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        // Set precision for this specific output block to 12 decimal places.
        outputFile << fixed << setprecision(12);

        // Label the first step (iter=0) as "Initial Guess" and subsequent steps as "Iteration X".
        if (iter == 0) {
            outputFile << "  Init. Guess: H1_Hrms = " << H1_Hrms << ", H2_Hrms = " << H2_Hrms
                       << ", F1 = " << f1_val << ", F2 = " << f2_val << endl;
        } else {
            outputFile << "  Iteration " << iter << ": H1_Hrms = " << H1_Hrms << ", H2_Hrms = " << H2_Hrms
                       << ", F1 = " << f1_val << ", F2 = " << f2_val << endl;
        }

        // Check for convergence. If converged, return the number of iterations (iter).
        if (abs(f1_val) < EPSILON && abs(f2_val) < EPSILON) {
            return iter;
        }

        // If we have reached the maximum number of iterations without converging, exit.
        if (iter == maxit) {
            break;
        }

        // Calculate Jacobian and solve for the next step if not converged.
        long double J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        long double J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        long double J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0L * JACOBIAN_DX);
        long double J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0L * JACOBIAN_DX);

        long double dH1, dH2;
        solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2);

        H1_Hrms += dH1;
        H2_Hrms += dH2;

        if (H1_Hrms <= 0.0L) H1_Hrms = numeric_limits<long double>::min();
        if (H2_Hrms <= 0.0L) H2_Hrms = numeric_limits<long double>::min();
    }

    return -1; // Return -1 if the loop completes without convergence.
}


/**
 * @brief Main function of the program.
 *
 * This function iterates through Htr/Hrms values from 0.01 to 3.5,
 * calls the Newton-Raphson solver for each, and writes the detailed
 * results to "output.txt". Finally, it calculates and prints summary
 * statistics about the number of iterations to the console.
 *
 * @return 0 on success, 1 on file I/O error.
 */
int main()
{
    ofstream outputFile("output.txt");
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open output.txt for writing." << endl;
        return 1;
    }

    outputFile << "Newton-Raphson Solver Iteration Analysis:\n";
    outputFile << "==================================================\n";

    vector<int> iteration_counts;
    const int max_iterations = 100;

    // Loop from 0.01 to 3.5 with a step of 0.01
    // Use <= 3.51L to robustly handle floating-point comparisons
    for (long double Htr_Hrms_val = 0.01L; Htr_Hrms_val <= 3.51L; Htr_Hrms_val += 0.01L) {
        long double H1_normalized;
        long double H2_normalized;

        get_initial_guesses(Htr_Hrms_val, H1_normalized, H2_normalized);

        outputFile << "\nCalculating for Htr/Hrms = " << fixed << setprecision(2) << Htr_Hrms_val << endl;
        outputFile << "--------------------------------------------------\n";

        int iterations = newtonRaphsonSystemSolver(Htr_Hrms_val, H1_normalized, H2_normalized, max_iterations, outputFile);

        outputFile << "--------------------------------------------------\n";
        outputFile << "Summary for Htr/Hrms = " << fixed << setprecision(2) << Htr_Hrms_val << ": ";
        if (iterations != -1) {
            outputFile << "Converged in " << iterations << " iterations." << endl;
            outputFile << "Final H1_Hrms = " << fixed << setprecision(12) << H1_normalized
                       << ", Final H2_Hrms = " << H2_normalized << endl;
            iteration_counts.push_back(iterations);
        } else {
            outputFile << "Failed to converge within " << max_iterations << " iterations." << endl;
        }
    }

    outputFile.close();
    cout << "Detailed results have been written to output.txt" << endl;

    // Calculate and display statistics to the console and output.txt
    if (!iteration_counts.empty()) {
        long long sum_of_iterations = accumulate(iteration_counts.begin(), iteration_counts.end(), 0LL);
        // Explicitly cast iteration_counts.size() to double to resolve the -Wconversion warning.
        double average_iterations = static_cast<double>(sum_of_iterations) / static_cast<double>(iteration_counts.size());
        auto min_iter = *min_element(iteration_counts.begin(), iteration_counts.end());
        auto max_iter = *max_element(iteration_counts.begin(), iteration_counts.end());

        // Re-open file in append mode to add the summary at the end.
        outputFile.open("output.txt", std::ios_base::app);
        if(outputFile.is_open()){
            outputFile << "\n\n--- Iteration Statistics (for Htr/Hrms from 0.01 to 3.5) ---\n";
            outputFile << "Minimum iterations: " << min_iter << endl;
            outputFile << "Maximum iterations: " << max_iter << endl;
            outputFile << "Average iterations: " << fixed << setprecision(2) << average_iterations << endl;
            outputFile << "-----------------------------------------------------------\n";
            outputFile.close();
        }


        cout << "\n--- Iteration Statistics (for Htr/Hrms from 0.01 to 3.5) ---\n";
        cout << "Minimum iterations: " << min_iter << endl;
        cout << "Maximum iterations: " << max_iter << endl;
        cout << "Average iterations: " << fixed << setprecision(2) << average_iterations << endl;
        cout << "-----------------------------------------------------------\n";
    } else {
        cout << "\nNo successful convergences to calculate statistics." << endl;
    }

    return 0;
}
