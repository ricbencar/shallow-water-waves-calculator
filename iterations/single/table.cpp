/*
 * table.cpp - Composite Weibull RMS closure lookup table
 * tau = Htr/Hrms; q is the nonnegative transition exponent.
 * At tau=0, q=0 is the analytic boundary limit (not a Brent root).
 * Uses Newton-Raphson iteration in z=ln(q), initialized by the GEP
 * approximation, and the existing incomplete-gamma algorithms.
 *
 * Output: table_tau_q.txt, two tab-separated columns, fixed decimal notation.
 * Build: g++ -O3 -std=c++17 -Wall -Wextra -pedantic table.cpp -o table
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace cwd {
constexpr double NUMERICAL_TOL = 2.5e-14;
constexpr int MAX_GAMMA_ITERATIONS = 10000;
constexpr int MAX_ROOT_ITERATIONS = 40;
static void require_finite_positive(double value, const char* name) {
    if (!(std::isfinite(value) && value > 0.0)) {
        throw std::runtime_error(std::string(name) + " must be finite and positive.");
    }
}

struct RegularizedGamma {
    double P = 0.0;
    double Q = 1.0;
};

// Same gamma series / continued fraction and numerical settings as the GUI.
static RegularizedGamma regularized_gamma(double a, double x) {
    if (!std::isfinite(a) || a <= 0.0 || !std::isfinite(x) || x < 0.0) {
        throw std::invalid_argument("Invalid incomplete gamma arguments.");
    }
    if (x == 0.0) return {0.0, 1.0};

    const double prefactor = std::exp(a * std::log(x) - x - std::lgamma(a));
    double P = 0.0;
    double Q = 0.0;
    bool converged = false;

    if (x < a + 1.0) {
        double ap = a;
        double delta = 1.0 / a;
        double sum = delta;
        for (int n = 1; n <= MAX_GAMMA_ITERATIONS; ++n) {
            ap += 1.0;
            delta *= x / ap;
            sum += delta;
            if (!std::isfinite(sum)) {
                throw std::runtime_error("Gamma series became non-finite.");
            }
            if (std::abs(delta) <= NUMERICAL_TOL * std::abs(sum)) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Incomplete gamma series did not converge.");
        }
        P = sum * prefactor;
        Q = 1.0 - P;
    } else {
        constexpr double FPMIN = 1e-300;
        double b = x + 1.0 - a;
        double c = 1.0 / FPMIN;
        double d = 1.0 / b;
        double fraction = d;
        for (int n = 1; n <= MAX_GAMMA_ITERATIONS; ++n) {
            const double nd = static_cast<double>(n);
            const double an = -nd * (nd - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < FPMIN) d = std::copysign(FPMIN, d);
            c = b + an / c;
            if (std::abs(c) < FPMIN) c = std::copysign(FPMIN, c);
            d = 1.0 / d;
            const double delta = c * d;
            fraction *= delta;
            if (!std::isfinite(fraction)) {
                throw std::runtime_error("Gamma continued fraction became non-finite.");
            }
            if (std::abs(delta - 1.0) <= NUMERICAL_TOL) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Incomplete gamma continued fraction did not converge.");
        }
        Q = fraction * prefactor;
        P = 1.0 - Q;
    }

    constexpr double ROUNDING_ALLOWANCE = 2e-13;
    if (!std::isfinite(P) || !std::isfinite(Q) ||
        P < -ROUNDING_ALLOWANCE || P > 1.0 + ROUNDING_ALLOWANCE ||
        Q < -ROUNDING_ALLOWANCE || Q > 1.0 + ROUNDING_ALLOWANCE) {
        throw std::runtime_error("Incomplete gamma outside probability range.");
    }
    return {std::clamp(P, 0.0, 1.0), std::clamp(Q, 0.0, 1.0)};
}

static double gamma_lower(double a, double x) {
    return std::tgamma(a) * regularized_gamma(a, x).P;
}

static double gamma_upper(double a, double x) {
    return std::tgamma(a) * regularized_gamma(a, x).Q;
}


// The supplied GEP approximation is used only to initialize Newton's method.
static double initial_guess(double tau) {
    return 1.002054709 * tau * tau *
           std::tanh(std::pow(tau, 1.652259933)) *
           std::exp(-std::tanh(std::tanh(tau)) *
                    std::exp(-2.146235748 * tau));
}

struct NewtonEvaluation {
    double q = 0.0;
    double residual = 0.0;
    double derivative = 0.0;
};

// R(z) = 2*log(tau) + log(L+U), q=exp(z), with
// L=gamma(2,q)/q and U=q^(-5/9)*Gamma(14/9,q).
// The exp(-q) terms in dL/dq and dU/dq cancel exactly, giving
// dR/dz = -(L+(5/9)*U)/(L+U); its value is between -1 and -5/9.
static NewtonEvaluation evaluate(double z, double log_tau) {
    NewtonEvaluation e;
    e.q = std::exp(z);
    require_finite_positive(e.q, "q");
    const double lower = std::exp(-z) * gamma_lower(2.0, e.q);
    const double upper = std::exp(-(5.0 / 9.0) * z) * gamma_upper(14.0 / 9.0, e.q);
    const double kernel = lower + upper;
    require_finite_positive(kernel, "RMS kernel");
    e.residual = 2.0 * log_tau + std::log(kernel);
    e.derivative = -(lower + (5.0 / 9.0) * upper) / kernel;
    if (!std::isfinite(e.residual) || !std::isfinite(e.derivative) ||
        !(e.derivative < 0.0)) {
        throw std::runtime_error("Invalid Newton residual or derivative");
    }
    return e;
}

// Newton-Raphson on z=ln(q): q remains strictly positive at every step.
// At tau=0 the exact closure is singular; use its analytic limit instead.
static double solve_q(double tau) {
    if (tau == 0.0) return 0.0;
    require_finite_positive(tau, "tau");
    const double q0 = initial_guess(tau);
    require_finite_positive(q0, "GEP initial guess");
    double z = std::log(q0);
    const double log_tau = std::log(tau);
    constexpr double ROOT_TOL = 5e-13;
    for (int iteration = 0; iteration <= MAX_ROOT_ITERATIONS; ++iteration) {
        const NewtonEvaluation e = evaluate(z, log_tau);
        if (std::abs(e.residual) <= ROOT_TOL) return e.q;
        if (iteration == MAX_ROOT_ITERATIONS) break;
        const double step = -e.residual / e.derivative;
        z += step;
        if (!std::isfinite(z))
            throw std::runtime_error("Newton update became non-finite");
    }
    throw std::runtime_error("Newton solver failed to converge");
}

} // namespace cwd

int main() {
    try {
        std::ofstream output("table_tau_q.txt", std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot open table_tau_q.txt");
        output << std::fixed << std::setprecision(18);
        output << "tau\tq\n";
        for (int i = 0; i <= 350; ++i) {
            const double tau = static_cast<double>(i) / 100.0;
            const double q = cwd::solve_q(tau);
            output << std::setprecision(2) << tau << '\t'
                   << std::setprecision(18) << q << '\n';
        }
        if (!output) throw std::runtime_error("Error writing table_tau_q.txt");
        std::cout << "Saved 351 tau-q pairs to table_tau_q.txt\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
