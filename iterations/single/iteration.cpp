/***********************************************************************
 * Program: iteration.cpp
 *
 * Composite Weibull RMS closure: scalar Newton iteration analysis.
 *
 * For each tau = Htr/Hrms from 0.00 to 3.50 in increments of 0.01,
 * solve the exact second-moment condition for tau>0, and use the
 * analytic limiting solution for tau=0:
 * q=0, H2/Hrms=1/sqrt(Gamma(14/9)), and H1/Hrms diverges.
 *
 *   tau^2 * [gamma(2,q)/q + q^(-5/9)*Gamma(14/9,q)] = 1,
 *
 * where gamma and Gamma are UNNORMALIZED lower and upper incomplete
 * gamma functions.  gamma(2,q) = 1 - (1+q)*exp(-q), evaluated by a
 * stable gamma series for small q to avoid cancellation.
 *
 * The solver is Newton-Raphson in z = log(q), initialized by the
 * supplied GEP q(tau) approximation; the RMS closure remains exact.
 *
 * For each tau, output.txt contains the initial guess, every Newton step, and final
 * normalized scales, RMS residual, and continuity residual. Summary
 * iteration statistics are also printed on the console.
 *
 * Build: g++ -O3 -std=c++17 -Wall -Wextra -pedantic \
 *            -Wconversion -Wsign-conversion iteration.cpp -o iteration
 ***********************************************************************/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cwd {

constexpr double K1 = 2.0;
constexpr double K2 = 3.6;
constexpr double A_LOWER = 2.0;
constexpr double A_UPPER = 14.0 / 9.0;
constexpr double UPPER_POWER = 5.0 / 9.0;
constexpr double NUMERICAL_TOL = 2.5e-14;
constexpr int MAX_GAMMA_ITERATIONS = 10000;
constexpr int MAX_ROOT_ITERATIONS = 40;
constexpr double CONSTRAINT_TOL = 1e-10;

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

struct MomentEvaluation {
    double z = 0.0;
    double q = 0.0;
    double lower_term = 0.0;
    double upper_term = 0.0;
    double log_residual = 0.0;
    double rms_residual = 0.0;
    double log_derivative = 0.0;
};

// Exact complete scalar RMS closure, evaluated in logarithmic coordinates:
// R(z) = 2*log(tau) + log(gamma(2,q)/q + q^(-5/9)*Gamma(14/9,q)).
static MomentEvaluation evaluate(double z, double log_tau) {
    MomentEvaluation e;
    e.z = z;
    e.q = std::exp(z);
    require_finite_positive(e.q, "Transition exponent q");
    e.lower_term = std::exp(-z) * gamma_lower(A_LOWER, e.q);
    e.upper_term = std::exp(-UPPER_POWER * z) * gamma_upper(A_UPPER, e.q);
    const double kernel = e.lower_term + e.upper_term;
    require_finite_positive(kernel, "RMS closure kernel");
    e.log_residual = 2.0 * log_tau + std::log(kernel);
    // Analytic derivative: the exp(-q) contributions cancel identically.
    e.log_derivative = -(e.lower_term + UPPER_POWER * e.upper_term) / kernel;
    // expm1 is accurate when log_residual is close to zero.
    e.rms_residual = std::expm1(0.5 * e.log_residual);
    if (!std::isfinite(e.log_residual) || !std::isfinite(e.rms_residual) ||
        !std::isfinite(e.log_derivative) || !(e.log_derivative < 0.0)) {
        throw std::runtime_error("Non-finite closure residual.");
    }
    return e;
}

struct Solution {
    int iterations = 0;
    double q = 0.0;
    double H1_over_Hrms = 0.0;
    double H2_over_Hrms = 0.0;
    double log_residual = 0.0;
    double moment_residual = 0.0;
    double continuity_residual = 0.0;
};

// The GEP expression initializes Newton but is never substituted for the
// solution of the complete RMS closure.
static double initial_guess(double tau) {
    return 1.002054709 * tau * tau *
           std::tanh(std::pow(tau, 1.652259933)) *
           std::exp(-std::tanh(std::tanh(tau)) *
                    std::exp(-2.146235748 * tau));
}

static void log_iteration(std::ostream& out, int iteration,
                          const char* step_type, const MomentEvaluation& e,
                          double tau, double delta_z) {
    const double a1 = tau * std::exp(-e.z / K1);
    const double a2 = tau * std::exp(-e.z / K2);
    out << "  " << std::setw(4) << iteration << "  "
        << std::left << std::setw(12) << step_type << std::right
        << std::scientific << std::setprecision(12)
        << "  " << std::setw(19) << e.z
        << "  " << std::setw(19) << e.q
        << "  " << std::setw(19) << e.log_residual
        << "  " << std::setw(19) << e.log_derivative
        << "  " << std::setw(19) << delta_z
        << "  " << std::setw(19) << e.rms_residual
        << "  " << std::setw(19) << a1
        << "  " << std::setw(19) << a2 << '\n';
}

// Pure Newton-Raphson in z=ln(q), preserving q>0. The analytic derivative
// R'(z)=-(L+(5/9)*U)/(L+U) is strictly negative and bounded away from zero.
static Solution solve(double tau, std::ostream& out) {
    if (tau == 0.0) {
        // The untransformed RMS equation is singular at tau=q=0.
        // Its positive-side limit satisfies
        // q ~ [Gamma(14/9)*tau^2]^(9/5), so q -> 0.
        const double gamma_value = std::tgamma(A_UPPER);
        require_finite_positive(gamma_value, "Gamma(14/9)");
        Solution s;
        s.iterations = 0;
        s.q = 0.0;
        s.H1_over_Hrms = std::numeric_limits<double>::infinity();
        s.H2_over_Hrms = 1.0 / std::sqrt(gamma_value);
        out << "  Analytic boundary limit: no Newton iterations.\n"
            << "  q -> 0 as tau -> 0+.\n"
            << "  H1/Hrms diverges; H2/Hrms has a finite limit.\n"
            << "  Closure residuals reported below are limiting values.\n";
        return s;
    }
    require_finite_positive(tau, "Normalized transition tau");
    const double q0 = initial_guess(tau);
    require_finite_positive(q0, "GEP initial guess");
    double z = std::log(q0);
    const double log_tau = std::log(tau);
    constexpr double ROOT_TOL = 5e-13;
    out << std::scientific << std::setprecision(12)
        << "  Initial GEP q = " << q0 << "; initial z = " << z << '\n'
        << "  " << std::setw(4) << "Iter" << "  "
        << std::left << std::setw(12) << "Step" << std::right
        << "  " << std::setw(19) << "z=ln(q)"
        << "  " << std::setw(19) << "q"
        << "  " << std::setw(19) << "R(log moment)"
        << "  " << std::setw(19) << "dR/dz"
        << "  " << std::setw(19) << "delta z"
        << "  " << std::setw(19) << "R(RMS)"
        << "  " << std::setw(19) << "H1/Hrms"
        << "  " << std::setw(19) << "H2/Hrms" << '\n';

    for (int iteration = 0; iteration <= MAX_ROOT_ITERATIONS; ++iteration) {
        const MomentEvaluation e = evaluate(z, log_tau);
        if (iteration == 0)
            log_iteration(out, 0, "GEP initial", e, tau, 0.0);
        if (std::abs(e.log_residual) <= ROOT_TOL) {
            Solution s;
            s.iterations = iteration;
            s.q = e.q;
            s.H1_over_Hrms = tau * std::exp(-z / K1);
            s.H2_over_Hrms = tau * std::exp(-z / K2);
            require_finite_positive(s.H1_over_Hrms, "First CWD scale");
            require_finite_positive(s.H2_over_Hrms, "Second CWD scale");
            s.log_residual = e.log_residual;
            const double m2 =
                s.H1_over_Hrms * s.H1_over_Hrms * gamma_lower(A_LOWER, s.q) +
                s.H2_over_Hrms * s.H2_over_Hrms * gamma_upper(A_UPPER, s.q);
            require_finite_positive(m2, "Reconstructed second moment");
            s.moment_residual = std::sqrt(m2) - 1.0;
            s.continuity_residual =
                std::pow(tau / s.H1_over_Hrms, K1) -
                std::pow(tau / s.H2_over_Hrms, K2);
            if (!std::isfinite(s.moment_residual) ||
                !std::isfinite(s.continuity_residual) ||
                std::abs(s.moment_residual) > CONSTRAINT_TOL ||
                std::abs(s.continuity_residual) > CONSTRAINT_TOL) {
                throw std::runtime_error("CWD RMS or continuity constraint not satisfied.");
            }
            return s;
        }
        if (iteration == MAX_ROOT_ITERATIONS) break;
        const double delta_z = -e.log_residual / e.log_derivative;
        z += delta_z;
        if (!std::isfinite(z))
            throw std::runtime_error("Non-finite Newton update.");
        const MomentEvaluation next = evaluate(z, log_tau);
        log_iteration(out, iteration + 1, "Newton", next, tau, delta_z);
    }
    throw std::runtime_error("Newton solver failed to converge.");
}

}  // namespace cwd

int main() {
    std::ofstream output("output.txt", std::ios::trunc);
    if (!output) {
        std::cerr << "Error: cannot create output.txt\n";
        return 1;
    }
    output.imbue(std::locale::classic());
    std::cout.imbue(std::locale::classic());
    output << "COMPOSITE WEIBULL RMS CLOSURE - NEWTON ITERATION ANALYSIS\n"
           << "=======================================================\n"
           << "Closure: tau^2 * [gamma(2,q)/q + q^(-5/9)*Gamma(14/9,q)] = 1\n"
           << "Unknown for tau>0: z = ln(q); q = exp(z) > 0\n"
           << "Method: Newton-Raphson in z=ln(q), initialized from GEP\n"
           << "Boundary at tau=0: analytic limit q=0; no Newton iteration\n"
           << "Sweep: tau = 0.00 to 3.50, inclusive, step 0.01 (351 cases)\n"
           << "NOTE: Exact CWD closure for tau>0, without a Rayleigh switch.\n";

    std::vector<int> iteration_counts;
    iteration_counts.reserve(350);
    int failures = 0;
    int successes = 0;
    int boundary_cases = 0;
    double worst_moment_residual = 0.0;
    double worst_continuity_residual = 0.0;

    // Integer iteration prevents floating-point drift at the endpoints.
    for (int index = 0; index <= 350; ++index) {
        const double tau = static_cast<double>(index) / 100.0;
        output << "\n\nCalculating for tau = Htr/Hrms = "
               << std::fixed << std::setprecision(2) << tau << '\n'
               << std::string(184, '-') << '\n';
        try {
            const cwd::Solution s = cwd::solve(tau, output);
            output << std::scientific << std::setprecision(15);
            if (index == 0) {
                ++boundary_cases;
                output << "  Boundary limit (0 Newton iterations).\n";
            } else {
                iteration_counts.push_back(s.iterations);
                output << "  Converged in " << s.iterations << " Newton iterations.\n";
            }
            output << "  q                           = " << s.q << '\n';
            if (index == 0) {
                output << "  H1/Hrms                     = divergent (+infinity as tau -> 0+)\n";
            } else {
                output << "  H1/Hrms                     = " << s.H1_over_Hrms << '\n';
            }
            output << "  H2/Hrms                     = " << s.H2_over_Hrms << '\n'
                   << "  Logarithmic root residual   = " << s.log_residual
                   << (index == 0 ? " (limit)\n" : "\n")
                   << "  RMS residual sqrt(m2)-1     = " << s.moment_residual
                   << (index == 0 ? " (limit)\n" : "\n")
                   << "  Continuity exponent residual= " << s.continuity_residual
                   << (index == 0 ? " (limit)\n" : "\n");
            ++successes;
            worst_moment_residual = std::max(worst_moment_residual,
                                             std::abs(s.moment_residual));
            worst_continuity_residual = std::max(worst_continuity_residual,
                                                 std::abs(s.continuity_residual));
        } catch (const std::exception& ex) {
            ++failures;
            output << "  FAILED: " << ex.what() << '\n';
            std::cerr << "tau=" << std::fixed << std::setprecision(2)
                      << tau << ": " << ex.what() << '\n';
        }
    }

    std::ostringstream summary;
    summary.imbue(std::locale::classic());
    summary << "\n--- Newton Statistics (tau = 0.00 to 3.50) ---\n"
            << "Cases evaluated:       351\n"
            << "Successful solutions:  " << successes << '\n'
            << "Boundary-limit cases:  " << boundary_cases << '\n'
            << "Newton-solved cases:    " << iteration_counts.size() << '\n'
            << "Failed solutions:      " << failures << '\n';
    if (!iteration_counts.empty()) {
        const auto minmax = std::minmax_element(iteration_counts.begin(),
                                                iteration_counts.end());
        const long long total = std::accumulate(iteration_counts.begin(),
                                                iteration_counts.end(), 0LL);
        const double average = static_cast<double>(total) /
                               static_cast<double>(iteration_counts.size());
        summary << "Minimum iterations:    " << *minmax.first << '\n'
                << "Average iterations:    " << std::fixed << std::setprecision(2)
                << average << '\n'
                << "Maximum iterations:    " << *minmax.second << '\n'
                << "Maximum |RMS residual|: " << std::scientific
                << std::setprecision(6) << worst_moment_residual << '\n'
                << "Maximum |cont. resid.|: " << worst_continuity_residual << '\n';
    }
    summary << "----------------------------------------------------\n";

    output << summary.str();
    output.flush();
    if (!output) {
        std::cerr << "Error: could not write output.txt completely.\n";
        return 1;
    }
    std::cout << "Detailed iteration results written to output.txt\n"
              << summary.str();
    return failures == 0 ? 0 : 1;
}
