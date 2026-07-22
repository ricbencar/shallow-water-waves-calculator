/***********************************************************************
 * Program: shallow-water-waves_gui.cpp
 *
 * Detailed Description:
 * This program computes local shallow-foreshore wave-height distribution
 * parameters using a model based on the Composed Weibull distribution.
 *
 * The application performs the following:
 * 1. The user provides Hm0 (local significant spectral wave height), 
 * d (local water depth), and slopeM (beach slope 1:m) via input fields.
 * * 2. Computes the following intermediate values:
 * - Free-surface variance (m0) is calculated directly using m0 = (Hm0/4)^2.
 * - Root-mean-square wave height (Hrms) is then calculated from m0
 * using the equation: Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0).
 * - A dimensional transitional wave height: Htr = (0.35 + 5.8*(1/m)) * d.
 * - The dimensionless transitional parameter: H̃_tr = Htr / Hrms.
 *
 * 3. Calculates the dimensionless wave-height ratios (Hᵢ/Hrms)
 * by solving a system of non-linear equations derived from the Composite Weibull
 * distribution, ensuring the normalized Hrms of the distribution equals one.
 * This involves using a Newton-Raphson matrix method for simultaneous root-finding
 * and functions for unnormalized incomplete gamma calculations. These ratios are then
 * converted to dimensional quantities (in meters) by multiplying with Hrms.
 *
 * 4. A detailed report is then generated (and written to "report.txt") with
 * the input parameters, intermediate values, calculated ratios and computed
 * dimensional wave heights, as well as diagnostic ratios.
 *
 * Overshoot-prevention logic:
 * To prevent the composite Weibull distribution from predicting wave heights
 * greater than the theoretical Rayleigh distribution (which is valid for deep
 * water), the following logic is applied:
 *
 * a) H_tr Threshold Switch: If the dimensionless transitional height (Htr/Hrms)
 * is greater than 2.75, the conditions are considered deep-water dominant.
 * The program bypasses the B&G calculation and uses established Rayleigh
 * distribution values for all H1/N statistics.
 *
 * b) Capping Statistical Parameters: If the B&G calculation is performed, the
 * resulting dimensional values for all H1/N are capped at their
 * theoretical Rayleigh limits (e.g., H1/3 <= Hm0, H1/10 <= 1.273*Hm0, etc.).
 *
 * Compilation Instructions (example using g++ on Windows):
 *
 * g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -municode shallow-water-waves_gui.cpp -o \
 * shallow-water-waves_gui -mwindows -static -static-libgcc -static-libstdc++
 *
 ***********************************************************************/

#define _USE_MATH_DEFINES // Required on some systems for M_PI

#include <windows.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <locale>
#include <codecvt>
#include <limits>
#include <algorithm>
#include <iterator>

// --- Global Constants ---

// Parameters for the Composite Weibull distribution.
constexpr double K1 = 2.0; // Exponent for the first part (Rayleigh-shaped).
constexpr double K2 = 3.6; // Exponent for the second part.

// Precision for Numerical Methods.
constexpr double EPSILON = 1e-12;     // Max allowable error for solver convergence.
constexpr double JACOBIAN_DX = 1e-8;  // Step size for finite difference Jacobian.
constexpr double LOCAL_EPS = 1e-16;   // Tolerance for gamma function convergence.

// --- Data Structures ---

/**
 * @struct WaveAnalysisResults
 * @brief Holds all input and calculated parameters for the wave analysis.
 * This structure simplifies passing data between functions.
 */
struct WaveAnalysisResults {
    // Inputs
    double Hm0;
    double d;
    double slopeM;

    // OVERSHOOT-PREVENTION: Added to track which model is used.
    std::string distribution_type;

    // Calculated Parameters
    double Hrms;
    double m0;
    double tanAlpha;
    double Htr_dim;
    double Htr_tilde;
    double H1_Hrms;
    double H2_Hrms;

    // Dimensionless Wave Heights (H/Hrms)
    double H1_3_Hrms;
    double H1_10_Hrms;
    double H1_50_Hrms;
    double H1_100_Hrms;
    double H1_250_Hrms;
    double H1_1000_Hrms;

    // Dimensional Wave Heights (m)
    double H1_dim;
    double H2_dim;
    double H1_3_dim;
    double H1_10_dim;
    double H1_50_dim;
    double H1_100_dim;
    double H1_250_dim;
    double H1_1000_dim;

    // Diagnostic Ratios
    double ratio_1_10_div_1_3;
    double ratio_1_50_div_1_3;
    double ratio_1_100_div_1_3;
    double ratio_1_250_div_1_3;
    double ratio_1_1000_div_1_3;
};


// --- Numerical Utility Functions ---

/**
 * @brief Computes the unnormalized lower incomplete gamma function γ(a, x).
 * @details This function is critical for calculating moments of Weibull distributions.
 * It uses a hybrid approach for numerical stability: a series expansion for x < a + 1.0,
 * and a continued fraction expansion for larger x.
 * @param a The shape parameter of the incomplete gamma function (must be positive).
 * @param x The upper integration limit of the incomplete gamma function (must be non-negative).
 * @return The computed value of γ(a, x), or NaN if inputs are invalid or convergence fails.
 */
double incomplete_gamma_lower(double a, double x) {
    constexpr int MAXIT = 500;

    if (a <= 0.0 || x < 0.0) return std::nan("");
    if (x == 0.0) return 0.0;

    double gln = std::lgamma(a);

    if (x < a + 1.0) { // Series expansion
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * LOCAL_EPS) {
                return sum * std::exp(-x + a * std::log(x) - gln) * std::tgamma(a);
            }
        }
        return std::nan(""); // Failed to converge
    } else { // Continued fraction
        double b = x + 1.0 - a;
        double c = 1.0 / std::numeric_limits<double>::min();
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= MAXIT; ++i) {
            double an = -1.0 * i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
            c = b + an / c;
            if (std::abs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < LOCAL_EPS) break;
        }
        double p_normalized = 1.0 - (std::exp(-x + a * std::log(x) - gln) * h);
        return std::max(0.0, p_normalized) * std::tgamma(a);
    }
}

/**
 * @brief Computes the unnormalized upper incomplete gamma function Γ(a, x).
 * @details Calculated using the identity Γ(a,x) = Γ(a) - γ(a,x).
 * @param a The shape parameter (must be positive).
 * @param x The lower integration limit (must be non-negative).
 * @return The computed value of Γ(a, x).
 */
static inline double incomplete_gamma_upper(double a, double x) {
    if (a <= 0.0 || x < 0.0) {
        throw std::invalid_argument("Invalid arguments for incomplete_gamma_upper.");
    }
    return std::tgamma(a) - incomplete_gamma_lower(a, x);
}

// --- Wave Height Calculation Functions ---

/**
 * @brief Calculates HN (wave height with 1/N exceedance probability).
 * @param N The N value (e.g., 3 for H1/3). Must be > 1.
 * @param H1 Scale parameter of the first Weibull distribution.
 * @param H2 Scale parameter of the second Weibull distribution.
 * @param Htr Transitional wave height.
 * @return The calculated value of HN.
 */
double calculate_HN(double N, double H1, double H2, double Htr) {
    if (N <= 1.0 || H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("Invalid arguments for calculate_HN.");
    }
    double HN_candidate1 = H1 * std::pow(std::log(N), 1.0 / K1);
    return (HN_candidate1 < Htr - EPSILON) ? HN_candidate1 : H2 * std::pow(std::log(N), 1.0 / K2);
}

/**
 * @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N).
 * @param N_val The N parameter for H1/N (e.g., 3 for H1/3). Must be > 1.
 * @param H1 Scale parameter of the first Weibull distribution.
 * @param H2 Scale parameter of the second Weibull distribution.
 * @param Htr Transitional wave height.
 * @return The calculated value of H1/N.
 */
double calculate_H1N(double N_val, double H1, double H2, double Htr) {
    if (N_val <= 1.0 || H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("Invalid arguments for calculate_H1N.");
    }

    double H_N_val = calculate_HN(N_val, H1, H2, Htr);
    double term1_x_ln_Nval = std::log(N_val);
    double term2_a = 1.0 / K2 + 1.0;

    if (H_N_val < Htr - EPSILON) {
        double term1_a = 1.0 / K1 + 1.0;
        double term1_x_HtrH1 = std::pow(Htr / H1, K1);
        double term2_x_HtrH2 = std::pow(Htr / H2, K2);

        double gamma1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval) - incomplete_gamma_upper(term1_a, term1_x_HtrH1);
        double gamma2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2);
        return N_val * (H1 * gamma1 + H2 * gamma2);
    } else {
        return N_val * H2 * incomplete_gamma_upper(term2_a, term1_x_ln_Nval);
    }
}


// --- Newton-Raphson Solver for H1/Hrms and H2/Hrms ---

/**
 * @brief First non-linear equation F1 for the solver.
 * @details Represents the normalized Hrms constraint: sqrt(Integral) - 1 = 0.
 */
static inline double F1(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) return std::numeric_limits<double>::max();
    double arg1 = std::pow(Htr_Hrms / H1_Hrms, K1);
    double arg2 = std::pow(Htr_Hrms / H2_Hrms, K2);
    double term1 = H1_Hrms * H1_Hrms * incomplete_gamma_lower(2.0 / K1 + 1.0, arg1);
    double term2 = H2_Hrms * H2_Hrms * incomplete_gamma_upper(2.0 / K2 + 1.0, arg2);
    return std::sqrt(std::max(0.0, term1 + term2)) - 1.0;
}

/**
 * @brief Second non-linear equation F2 for the solver.
 * @details Represents the continuity condition at the transitional height Htr.
 */
static inline double F2(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) return std::numeric_limits<double>::max();
    return std::pow(Htr_Hrms / H1_Hrms, K1) - std::pow(Htr_Hrms / H2_Hrms, K2);
}

/**
 * @brief Solves a 2x2 linear system Ax = b using Cramer's rule.
 */
static inline void solve_linear_system_2x2(double J11, double J12, double J21, double J22,
                                           double b1, double b2, double &dx1, double &dx2) {
    double determinant = J11 * J22 - J12 * J21;
    if (std::abs(determinant) < 1e-20) { // Avoid division by near-zero
        dx1 = 0.0;
        dx2 = 0.0;
        return;
    }
    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
}

/**
 * @brief Provides empirical initial guesses for H1_Hrms and H2_Hrms.
 */
static void get_initial_guesses(double Htr_Hrms, double &H1_initial, double &H2_initial) {
    H1_initial = 0.9718670705250743 + 1.115952604282648 * std::pow(Htr_Hrms, -0.7970446117540275) * std::exp(-1.449005086812895 * Htr_Hrms);
    H2_initial = 1.059259665431797 + (0.2059286860468916 * Htr_Hrms) / (1.0 + 3.865701948059343 * std::pow(Htr_Hrms, -3.479682433107255));
    if (H1_initial <= 0.0) H1_initial = LOCAL_EPS;
    if (H2_initial <= 0.0) H2_initial = LOCAL_EPS;
}

/**
 * @brief Solves for H1_Hrms and H2_Hrms using the Newton-Raphson method for systems.
 * @return `true` if the solver converges, `false` otherwise.
 */
bool newtonRaphsonSystemSolver(double Htr_Hrms, double &H1_Hrms, double &H2_Hrms, int maxit) {
    get_initial_guesses(Htr_Hrms, H1_Hrms, H2_Hrms);

    for (int iter = 0; iter < maxit; ++iter) {
        double f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        double f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        if (std::abs(f1_val) < EPSILON && std::abs(f2_val) < EPSILON) {
            return true;
        }

        double J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);

        double dH1, dH2;
        solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2);

        H1_Hrms += dH1;
        H2_Hrms += dH2;

        if (H1_Hrms <= 0.0) H1_Hrms = LOCAL_EPS;
        if (H2_Hrms <= 0.0) H2_Hrms = LOCAL_EPS;
    }

    // In a GUI, this error is reported via the return value and handled in the event loop.
    return false;
}


// --- Core Calculation Logic ---

/**
 * @brief Performs the full wave analysis, populating the results structure.
 * @param results A reference to the results structure to be filled.
 * @return `true` on success, `false` on failure.
 */
bool perform_wave_analysis(WaveAnalysisResults& results) {
    try {
        // Step 1: Calculate primary parameters
        results.m0 = std::pow(results.Hm0 / 4.0, 2.0);
        const double sqrt_m0 = std::sqrt(results.m0);
        results.Hrms = (2.69 + 3.24 * sqrt_m0 / results.d) * sqrt_m0;

        results.tanAlpha = 1.0 / results.slopeM;
        results.Htr_dim = (0.35 + 5.8 * results.tanAlpha) * results.d;
        results.Htr_tilde = (results.Hrms > 0.0) ? (results.Htr_dim / results.Hrms) : 0.0;

        // OVERSHOOT-PREVENTION: Method 1 - H_tr Threshold Switch
        if (results.Htr_tilde > 2.75) {
            results.distribution_type = "Rayleigh";
            
            // Use theoretically exact H(1/N)/Hm0 ratios for a pure Rayleigh distribution.
            results.H1_3_dim    = 1.001075736951740 * results.Hm0;
            results.H1_10_dim   = 1.272734273369137 * results.Hm0;
            results.H1_50_dim   = 1.560113379974762 * results.Hm0;
            results.H1_100_dim  = 1.668233372358517 * results.Hm0;
            results.H1_250_dim  = 1.801017222497626 * results.Hm0;
            results.H1_1000_dim = 1.984835590575388 * results.Hm0;
            
            // Back-calculate dimensionless ratios for the report
            if (results.Hrms > 0.0) {
                results.H1_3_Hrms = results.H1_3_dim / results.Hrms;
                results.H1_10_Hrms = results.H1_10_dim / results.Hrms;
                results.H1_50_Hrms = results.H1_50_dim / results.Hrms;
                results.H1_100_Hrms = results.H1_100_dim / results.Hrms;
                results.H1_250_Hrms = results.H1_250_dim / results.Hrms;
                results.H1_1000_Hrms = results.H1_1000_dim / results.Hrms;
            }
        } else {
            results.distribution_type = "B&G";

            // Step 2: Solve for H1/Hrms and H2/Hrms
            if (!newtonRaphsonSystemSolver(results.Htr_tilde, results.H1_Hrms, results.H2_Hrms, 100)) {
                 throw std::runtime_error("Newton-Raphson solver failed to converge.");
            }

            // Step 3: Calculate H1/N quantiles
            results.H1_3_Hrms = calculate_H1N(3.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_10_Hrms = calculate_H1N(10.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_50_Hrms = calculate_H1N(50.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_100_Hrms = calculate_H1N(100.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_250_Hrms = calculate_H1N(250.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_1000_Hrms = calculate_H1N(1000.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);

            // Step 4: Convert to dimensional values
            results.H1_dim = results.H1_Hrms * results.Hrms;
            results.H2_dim = results.H2_Hrms * results.Hrms;
            results.H1_3_dim = results.H1_3_Hrms * results.Hrms;
            results.H1_10_dim = results.H1_10_Hrms * results.Hrms;
            results.H1_50_dim = results.H1_50_Hrms * results.Hrms;
            results.H1_100_dim = results.H1_100_Hrms * results.Hrms;
            results.H1_250_dim = results.H1_250_Hrms * results.Hrms;
            results.H1_1000_dim = results.H1_1000_Hrms * results.Hrms;

            // OVERSHOOT-PREVENTION: Method 2 - Capping
            // Use theoretically exact H(1/N)/Hm0 ratios for a pure Rayleigh distribution.
            results.H1_3_dim    = std::min(results.H1_3_dim,    1.001075736951740 * results.Hm0);
            results.H1_10_dim   = std::min(results.H1_10_dim,   1.272734273369137 * results.Hm0);
            results.H1_50_dim   = std::min(results.H1_50_dim,   1.560113379974762 * results.Hm0);
            results.H1_100_dim  = std::min(results.H1_100_dim,  1.668233372358517 * results.Hm0);
            results.H1_250_dim  = std::min(results.H1_250_dim,  1.801017222497626 * results.Hm0);
            results.H1_1000_dim = std::min(results.H1_1000_dim, 1.984835590575388 * results.Hm0);
        }

        // Step 5: Calculate diagnostic ratios (always do this last)
        if (results.H1_3_dim > 0.0) { // Use capped dimensional value for check
            // Recalculate dimensionless ratios from potentially capped values for accurate reporting
            double H1_3_Hrms_capped = results.H1_3_dim / results.Hrms;
            double H1_10_Hrms_capped = results.H1_10_dim / results.Hrms;
            double H1_50_Hrms_capped = results.H1_50_dim / results.Hrms;
            double H1_100_Hrms_capped = results.H1_100_dim / results.Hrms;
            double H1_250_Hrms_capped = results.H1_250_dim / results.Hrms;
            double H1_1000_Hrms_capped = results.H1_1000_dim / results.Hrms;

            results.ratio_1_10_div_1_3 = H1_10_Hrms_capped / H1_3_Hrms_capped;
            results.ratio_1_50_div_1_3 = H1_50_Hrms_capped / H1_3_Hrms_capped;
            results.ratio_1_100_div_1_3 = H1_100_Hrms_capped / H1_3_Hrms_capped;
            results.ratio_1_250_div_1_3 = H1_250_Hrms_capped / H1_3_Hrms_capped;
            results.ratio_1_1000_div_1_3 = H1_1000_Hrms_capped / H1_3_Hrms_capped;
        } else {
            results.ratio_1_10_div_1_3 = 0.0;
            results.ratio_1_50_div_1_3 = 0.0;
            results.ratio_1_100_div_1_3 = 0.0;
            results.ratio_1_250_div_1_3 = 0.0;
            results.ratio_1_1000_div_1_3 = 0.0;
        }

    } catch (const std::exception& e) {
        std::wstringstream ss;
        ss << L"ERROR during analysis: " << e.what();
        MessageBoxW(NULL, ss.str().c_str(), L"Calculation Error", MB_ICONERROR | MB_OK);
        return false;
    }
    return true;
}


// --- Reporting and I/O Functions ---

/**
 * @brief Formats the analysis results into a detailed report string.
 * @param results The structure containing all calculated data.
 * @return A formatted std::wstring report.
 */
std::wstring format_report(const WaveAnalysisResults& r) {
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

    ss << L"======================\n"
       << L"   INPUT PARAMETERS\n"
       << L"======================\n"
       << L"Hm0 (m)         : " << r.Hm0 << L"\n"
       << L"d (m)           : " << r.d << L"\n"
       << L"Beach slope (m) : " << r.slopeM << L"   (tan(alpha) = " << r.tanAlpha << L")\n\n"

       // OVERSHOOT-PREVENTION: Report which distribution was used.
       << L"Distribution Used : " << std::wstring(r.distribution_type.begin(), r.distribution_type.end()) << L"\n\n"

       << L"===========================\n"
       << L"   CALCULATED PARAMETERS\n"
       << L"===========================\n"
       << L"Free-surface variance m0 (m^2)   : " << r.m0 << L"\n"
       << L"Mean square wave height Hrms (m) : " << r.Hrms << L"\n"
       << L"Transitional wave height Htr (m) : " << r.Htr_dim << L"\n"
       << L"Dimensionless H~_tr (Htr/Hrms)   : " << r.Htr_tilde << L"\n\n"

       << L"=========================================\n"
       << L"   DIMENSIONLESS WAVE HEIGHTS (H/Hrms)\n"
       << L"=========================================\n"
       << L"H1/Hrms       : " << r.H1_Hrms << L"\n"
       << L"H2/Hrms       : " << r.H2_Hrms << L"\n"
       << L"H1/3 / Hrms   : " << r.H1_3_Hrms << L"\n"
       << L"H1/10 / Hrms  : " << r.H1_10_Hrms << L"\n"
       << L"H1/50 / Hrms  : " << r.H1_50_Hrms << L"\n"
       << L"H1/100 / Hrms : " << r.H1_100_Hrms << L"\n"
       << L"H1/250 / Hrms : " << r.H1_250_Hrms << L"\n"
       << L"H1/1000 /Hrms : " << r.H1_1000_Hrms << L"\n\n"

       << L"==================================\n"
       << L"   DIMENSIONAL WAVE HEIGHTS (m)\n"
       << L"==================================\n"
       << L"H1 (m)        : " << r.H1_dim << L"\n"
       << L"H2 (m)        : " << r.H2_dim << L"\n"
       << L"H1/3 (m)      : " << r.H1_3_dim << L"\n"
       << L"H1/10 (m)     : " << r.H1_10_dim << L"\n"
       << L"H1/50 (m)     : " << r.H1_50_dim << L"\n"
       << L"H1/100 (m)    : " << r.H1_100_dim << L"\n"
       << L"H1/250 (m)    : " << r.H1_250_dim << L"\n"
       << L"H1/1000 (m)   : " << r.H1_1000_dim << L"\n\n"

       << L"=======================\n"
       << L"   DIAGNOSTIC RATIOS\n"
       << L"=======================\n"
       << L"(H1/10)/(H1/3)   : " << r.ratio_1_10_div_1_3 << L"\n"
       << L"(H1/50)/(H1/3)   : " << r.ratio_1_50_div_1_3 << L"\n"
       << L"(H1/100)/(H1/3)  : " << r.ratio_1_100_div_1_3 << L"\n"
       << L"(H1/250)/(H1/3)  : " << r.ratio_1_250_div_1_3 << L"\n"
       << L"(H1/1000)/(H1/3) : " << r.ratio_1_1000_div_1_3 << L"\n\n"

       << L"End of Report\n";

    return ss.str();
}

/**
 * @brief Writes the report string to "report.txt" with UTF-8 encoding.
 */
static void write_report_to_file(const std::wstring &report) {
    std::wofstream ofs("report.txt");
    if (!ofs) {
        MessageBoxW(NULL, L"Could not open report.txt for writing.", L"File Error", MB_ICONERROR | MB_OK);
        return;
    }
    ofs.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>));
    ofs << report;
}

/**
 * @brief Convert newline characters for proper display in the edit control.
 */
static std::wstring fix_newlines_for_edit_control(const std::wstring &text) {
    std::wstring out;
    out.reserve(text.size() + 100);
    for (wchar_t c : text) {
        if (c == L'\n') {
            out.push_back(L'\r');
        }
        out.push_back(c);
    }
    return out;
}

// --- Win32 GUI Specifics ---

// Control IDs
#define IDC_EDIT_HM0 101
#define IDC_EDIT_D 102
#define IDC_EDIT_SLOPE 103
#define IDC_BUTTON_COMPUTE 104
#define IDC_OUTPUT 105

// Global Handles
HWND hEditHm0, hEditD, hEditSlope, hOutput;

/**
 * @brief Main window procedure to handle messages for the GUI.
 */
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    static HFONT hMonoFont = NULL;

    switch (msg) {
    case WM_CREATE:
    {
        CreateWindowW(L"STATIC", L"Hm0 (m):", WS_CHILD | WS_VISIBLE,
                     10, 10, 120, 20, hwnd, NULL, NULL, NULL);
        hEditHm0 = CreateWindowW(L"EDIT", L"2.5", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                140, 10, 100, 20, hwnd, (HMENU)IDC_EDIT_HM0, NULL, NULL);

        CreateWindowW(L"STATIC", L"d (m):", WS_CHILD | WS_VISIBLE,
                     10, 40, 120, 20, hwnd, NULL, NULL, NULL);
        hEditD = CreateWindowW(L"EDIT", L"5.0", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                              140, 40, 100, 20, hwnd, (HMENU)IDC_EDIT_D, NULL, NULL);

        CreateWindowW(L"STATIC", L"Beach slope m:", WS_CHILD | WS_VISIBLE,
                     10, 70, 120, 20, hwnd, NULL, NULL, NULL);
        hEditSlope = CreateWindowW(L"EDIT", L"100", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                  140, 70, 100, 20, hwnd, (HMENU)IDC_EDIT_SLOPE, NULL, NULL);

        CreateWindowW(L"BUTTON", L"Compute", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                     10, 100, 230, 30, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);

        hOutput = CreateWindowW(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                               ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY,
                               10, 140, 760, 440, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

        hMonoFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                               DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                               DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New");
        SendMessageW(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
    }
    break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE) {
            wchar_t buffer[64];
            WaveAnalysisResults results = {};

            GetWindowTextW(hEditHm0, buffer, 63);
            results.Hm0 = _wtof(buffer);
            
            GetWindowTextW(hEditD, buffer, 63);
            results.d = _wtof(buffer);

            GetWindowTextW(hEditSlope, buffer, 63);
            results.slopeM = _wtof(buffer);

            if (results.Hm0 <= 0.0 || results.d <= 0.0 || results.slopeM <= 0.0) {
                MessageBoxW(hwnd, L"All input values must be positive.", L"Input Error", MB_ICONERROR | MB_OK);
                break;
            }
            
            if (perform_wave_analysis(results)) {
                std::wstring report = format_report(results);
                write_report_to_file(report);
                std::wstring guiText = fix_newlines_for_edit_control(report);
                SetWindowTextW(hOutput, guiText.c_str());
            } else {
                // Error message is shown by perform_wave_analysis itself
                SetWindowTextW(hOutput, L"Calculation failed. See error message pop-up.");
            }
        }
        break;

    case WM_DESTROY:
        if (hMonoFont) {
            DeleteObject(hMonoFont);
        }
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProcW(hwnd, msg, wParam, lParam);
    }
    return 0;
}

/**
 * @brief Main entry point for the Windows application.
 */
int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow) {
    const wchar_t CLASS_NAME[] = L"ShallowWaveWindowClass";

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION);

    if (!RegisterClassExW(&wc)) {
        MessageBoxW(NULL, L"Window Registration Failed!", L"Error!", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    HWND hwnd = CreateWindowExW(
        0, CLASS_NAME, L"Shallow Wave Heights GUI",
        WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 640,
        NULL, NULL, hInstance, NULL);

    if (!hwnd) {
        MessageBoxW(NULL, L"Window Creation Failed!", L"Error!", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg = {};
    while (GetMessageW(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }

    return static_cast<int>(msg.wParam);
}
