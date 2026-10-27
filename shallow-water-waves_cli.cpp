/*
 * Shallow-Water Individual Wave-Height Statistics
 *
 * Command-line application for local individual zero-crossing wave-height
 * statistics based on the Composite Weibull Distribution of Battjes and
 * Groenendijk. The Rayleigh distribution is used for weakly depth-limited
 * conditions, with independent Rayleigh limits on the CWD mean-high-wave
 * statistics.
 *
 * Inputs: local spectral Hm0 (m), local still-water depth d (m), and
 * foreshore slope expressed as 1:M.
 *
 * Build (Linux/macOS):
 * g++ -O2 -std=c++17 -Wall -Wextra -pedantic \
 *     shallow-water-waves_cli.cpp -o shallow-water-waves_cli
 *
 * Build (MinGW-w64 / Windows console):
 * g++ -O2 -std=c++17 -Wall -Wextra -pedantic \
 *     shallow-water-waves_cli.cpp -o shallow-water-waves_cli.exe
 */

#include <algorithm>
#include <array>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <utility>
#include <cmath>
#include <cstddef>
#include <cwchar>
#include <cwctype>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cwd {

constexpr double K1 = 2.0;
constexpr double K2 = 3.6;
constexpr double RAYLEIGH_SWITCH = 2.75;
constexpr double NUMERICAL_TOL = 2.5e-14;
constexpr int MAX_GAMMA_ITERATIONS = 10000;
constexpr int MAX_ROOT_ITERATIONS = 40;

constexpr std::array<int, 6> N_VALUES{{3, 10, 50, 100, 250, 1000}};
constexpr std::array<double, 6> RAYLEIGH_MEAN_FACTORS{{
    1.001075736951740,
    1.272734273369137,
    1.560113379974762,
    1.668233372358517,
    1.801017222497626,
    1.984835590575388
}};

struct Statistic {
    int N = 0;
    double threshold = 0.0;              // H_N of the selected distribution
    double rayleigh_threshold = 0.0;     // Reference Rayleigh H_N
    double mean_unlimited = 0.0;         // Mean of highest 1/N, selected distribution
    double rayleigh_mean = 0.0;         // Corresponding Rayleigh mean
    double mean_final = 0.0;            // Mean after applicable Rayleigh limit
    double mean_normalized = 0.0;       // mean_final / Hrms
    bool rayleigh_limit_applied = false;
};

struct Analysis {
    double Hm0 = 0.0;
    double depth = 0.0;
    double slope_denominator = 0.0;
    double slope = 0.0;
    double m0 = 0.0;
    double Hrms = 0.0;
    double saturation = 0.0;
    double Htr = 0.0;
    double transition_normalized = 0.0;
    bool uses_cwd = false;
    bool slope_outside_calibration = false;
    bool saturation_outside_observed_range = false;
    double scale1_normalized = std::numeric_limits<double>::quiet_NaN();
    double scale2_normalized = std::numeric_limits<double>::quiet_NaN();
    double transition_exponent = std::numeric_limits<double>::quiet_NaN();
    double transition_exceedance_probability = std::numeric_limits<double>::quiet_NaN();
    double transition_recurrence_count = std::numeric_limits<double>::quiet_NaN();
    double moment_residual = std::numeric_limits<double>::quiet_NaN();
    double continuity_residual = std::numeric_limits<double>::quiet_NaN();
    int solver_iterations = 0;
    std::array<Statistic, 6> statistics{};
};

static void require_finite_positive(double value, const char* name) {
    if (!(std::isfinite(value) && value > 0.0)) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive.");
    }
}

static double checked_positive(double value, const char* name) {
    require_finite_positive(value, name);
    return value;
}

// Regularized gamma pair P(a,x) and Q(a,x), evaluated without subtracting
// nearly equal *unnormalized* gamma functions in the upper tail.
struct RegularizedGamma {
    double P = 0.0;
    double Q = 1.0;
};

static RegularizedGamma regularized_gamma(double a, double x) {
    if (!std::isfinite(a) || a <= 0.0 || !std::isfinite(x) || x < 0.0) {
        throw std::invalid_argument("Invalid argument in incomplete gamma function.");
    }
    if (x == 0.0) return {0.0, 1.0};

    const double prefactor_log = a * std::log(x) - x - std::lgamma(a);
    const double prefactor = std::exp(prefactor_log);
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
        throw std::runtime_error("Incomplete gamma result outside the probability range.");
    }
    return {std::clamp(P, 0.0, 1.0), std::clamp(Q, 0.0, 1.0)};
}

static double gamma_lower(double a, double x) {
    return std::tgamma(a) * regularized_gamma(a, x).P;
}

static double gamma_upper(double a, double x) {
    return std::tgamma(a) * regularized_gamma(a, x).Q;
}

// Stable evaluation of the difference Gamma(a,low)-Gamma(a,high), equal to
// the integral of u^(a-1)*exp(-u) between the two finite bounds.
static double gamma_interval(double a, double low, double high) {
    if (!std::isfinite(a) || a <= 0.0 ||
        !std::isfinite(low) || low < 0.0 ||
        !std::isfinite(high) || high < low) {
        throw std::invalid_argument("Invalid incomplete gamma integration interval.");
    }
    if (low == high) return 0.0;

    constexpr std::array<double, 8> abscissae{{
        0.09501250983763744, 0.28160355077925891,
        0.45801677765722739, 0.61787624440264375,
        0.75540440835500303, 0.86563120238783174,
        0.94457502307323258, 0.98940093499164993
    }};
    constexpr std::array<double, 8> weights{{
        0.18945061045506850, 0.18260341504492359,
        0.16915651939500254, 0.14959598881657673,
        0.12462897125553387, 0.09515851168249278,
        0.06225352393864789, 0.02715245941175409
    }};

    const double mid = low + 0.5 * (high - low);
    const double half_width = 0.5 * (high - low);
    double integral = 0.0;
    for (std::size_t i = 0; i < abscissae.size(); ++i) {
        const double displacement = half_width * abscissae[i];
        const double u1 = mid - displacement;
        const double u2 = mid + displacement;
        const double f1 = std::pow(u1, a - 1.0) * std::exp(-u1);
        const double f2 = std::pow(u2, a - 1.0) * std::exp(-u2);
        integral += weights[i] * (f1 + f2);
    }
    return checked_positive(integral * half_width, "Gamma interval integral");
}

// Composite Weibull RMS closure, solved in logarithmic coordinates.
// q0 is an explicit approximation; the returned q satisfies the complete moment equation.
static double initial_transition_exponent(double t) {
    return 1.002054709 * t * t * std::tanh(std::pow(t, 1.652259933)) *
           std::exp(-std::tanh(std::tanh(t)) * std::exp(-2.146235748 * t));
}

struct LogMomentEvaluation {
    double residual = 0.0;
    double derivative = 0.0;
};

// With L=q^(-2/K1)*gamma(1+2/K1,q) and
// U=q^(-2/K2)*Gamma(1+2/K2,q), R=2*ln(t)+ln(L+U).
// The incomplete-gamma boundary derivatives cancel exactly:
// dR/dln(q) = -((2/K1)*L + (2/K2)*U)/(L+U).
static LogMomentEvaluation evaluate_log_moment(double log_q, double log_t) {
    const double q = std::exp(log_q);
    require_finite_positive(q, "Transition exponent q");
    const double lower = std::exp(-(2.0 / K1) * log_q) *
                         gamma_lower(1.0 + 2.0 / K1, q);
    const double upper = std::exp(-(2.0 / K2) * log_q) *
                         gamma_upper(1.0 + 2.0 / K2, q);
    const double kernel = lower + upper;
    require_finite_positive(kernel, "Second-moment kernel");
    const double residual = 2.0 * log_t + std::log(kernel);
    const double derivative = -((2.0 / K1) * lower + (2.0 / K2) * upper) / kernel;
    if (!std::isfinite(residual) || !std::isfinite(derivative) ||
        !(derivative < 0.0)) {
        throw std::runtime_error("Invalid Newton residual or derivative.");
    }
    return {residual, derivative};
}

static double solve_transition_exponent(double t, int& iterations) {
    require_finite_positive(t, "Normalized transition height");
    const double q0 = initial_transition_exponent(t);
    require_finite_positive(q0, "Newton initial transition exponent");
    const double log_t = std::log(t);
    double z = std::log(q0);
    constexpr double ROOT_TOL = 5.0e-13;

    for (int iteration = 0; iteration <= MAX_ROOT_ITERATIONS; ++iteration) {
        const LogMomentEvaluation e = evaluate_log_moment(z, log_t);
        if (std::abs(e.residual) <= ROOT_TOL) {
            iterations = iteration;
            const double q = std::exp(z);
            require_finite_positive(q, "CWD transition exponent");
            return q;
        }
        if (iteration == MAX_ROOT_ITERATIONS) break;
        const double delta_z = -e.residual / e.derivative;
        z += delta_z;
        if (!std::isfinite(z)) {
            throw std::runtime_error("Non-finite Newton iteration.");
        }
    }
    throw std::runtime_error("CWD Newton solver did not converge.");
}

static double cwd_second_moment(double a1, double a2, double q) {
    return a1 * a1 * gamma_lower(1.0 + 2.0 / K1, q) +
           a2 * a2 * gamma_upper(1.0 + 2.0 / K2, q);
}

static double cwd_threshold_normalized(int N, double a1,
                                       double a2, double q) {
    const double logN = std::log(static_cast<double>(N));
    return (logN < q) ? a1 * std::pow(logN, 1.0 / K1)
                      : a2 * std::pow(logN, 1.0 / K2);
}

static double cwd_mean_highest_normalized(int N, double a1,
                                          double a2, double q) {
    const double logN = std::log(static_cast<double>(N));
    double tail_moment = 0.0;
    if (logN < q) {
        tail_moment = a1 * gamma_interval(1.0 + 1.0 / K1, logN, q) +
                      a2 * gamma_upper(1.0 + 1.0 / K2, q);
    } else {
        tail_moment = a2 * gamma_upper(1.0 + 1.0 / K2, logN);
    }
    return checked_positive(static_cast<double>(N) * tail_moment,
                            "CWD mean-high-wave height");
}

static Analysis calculate(double Hm0, double depth, double slope_denominator) {
    require_finite_positive(Hm0, "Hm0");
    require_finite_positive(depth, "Water depth");
    require_finite_positive(slope_denominator, "Slope denominator M");

    Analysis result;
    result.Hm0 = Hm0;
    result.depth = depth;
    result.slope_denominator = slope_denominator;
    result.slope = checked_positive(1.0 / slope_denominator, "Foreshore slope");
    result.slope_outside_calibration =
        slope_denominator < 20.0 || slope_denominator > 250.0;
    const double sigma = checked_positive(Hm0 / 4.0, "Spectral standard deviation");
    result.m0 = checked_positive(sigma * sigma, "Spectral variance m0");
    result.saturation = checked_positive(sigma / depth, "Degree of saturation");
    result.saturation_outside_observed_range =
        result.saturation < 0.025 || result.saturation > 0.16;
    result.Hrms = checked_positive((2.69 + 3.24 * result.saturation) * sigma,
                                    "Root-mean-square wave height");
    result.Htr = checked_positive((0.35 + 5.8 * result.slope) * depth,
                                  "Transitional wave height");
    result.transition_normalized = checked_positive(result.Htr / result.Hrms,
                                                    "Normalized transition height");
    result.uses_cwd = result.transition_normalized <= RAYLEIGH_SWITCH;

    if (result.uses_cwd) {
        const double t = result.transition_normalized;
        const double q = solve_transition_exponent(t, result.solver_iterations);
        result.transition_exponent = q;
        result.scale1_normalized = checked_positive(
            t * std::exp(-std::log(q) / K1), "First CWD scale");
        result.scale2_normalized = checked_positive(
            t * std::exp(-std::log(q) / K2), "Second CWD scale");
        result.transition_exceedance_probability = std::exp(-q);
        result.transition_recurrence_count = std::exp(q);

        const double m2 = cwd_second_moment(
            result.scale1_normalized, result.scale2_normalized, q);
        result.moment_residual = std::sqrt(m2) - 1.0;
        const double exp1 = std::pow(t / result.scale1_normalized, K1);
        const double exp2 = std::pow(t / result.scale2_normalized, K2);
        result.continuity_residual = exp1 - exp2;

        if (!std::isfinite(m2) ||
            std::abs(result.moment_residual) > 1e-10 ||
            std::abs(result.continuity_residual) > 1e-10) {
            throw std::runtime_error("CWD probability constraints were not satisfied.");
        }
    }

    const double rayleigh_scale = checked_positive(Hm0 / std::sqrt(2.0),
                                                   "Rayleigh height scale");
    for (std::size_t i = 0; i < N_VALUES.size(); ++i) {
        Statistic& stat = result.statistics[i];
        stat.N = N_VALUES[i];
        const double logN = std::log(static_cast<double>(stat.N));
        stat.rayleigh_threshold = checked_positive(
            rayleigh_scale * std::sqrt(logN), "Rayleigh exceedance threshold");
        stat.rayleigh_mean = checked_positive(
            RAYLEIGH_MEAN_FACTORS[i] * Hm0, "Rayleigh high-wave mean");

        if (result.uses_cwd) {
            const double a1 = result.scale1_normalized;
            const double a2 = result.scale2_normalized;
            const double q = result.transition_exponent;
            stat.threshold = checked_positive(
                cwd_threshold_normalized(stat.N, a1, a2, q) * result.Hrms,
                "CWD exceedance threshold");
            stat.mean_unlimited = checked_positive(
                cwd_mean_highest_normalized(stat.N, a1, a2, q) * result.Hrms,
                "CWD high-wave mean");
            stat.mean_final = std::min(stat.mean_unlimited, stat.rayleigh_mean);
            stat.rayleigh_limit_applied =
                stat.mean_unlimited > stat.rayleigh_mean *
                (1.0 + 8.0 * std::numeric_limits<double>::epsilon());
        } else {
            stat.threshold = stat.rayleigh_threshold;
            stat.mean_unlimited = stat.rayleigh_mean;
            stat.mean_final = stat.rayleigh_mean;
        }

        stat.mean_normalized = checked_positive(stat.mean_final / result.Hrms,
                                                "Normalized high-wave mean");
        if (!(stat.mean_unlimited > stat.threshold) ||
            stat.mean_final > stat.rayleigh_mean * (1.0 + 1e-12)) {
            throw std::runtime_error("Inconsistent wave-height statistics.");
        }
        if (i != 0) {
            if (!(stat.threshold > result.statistics[i - 1].threshold) ||
                !(stat.mean_unlimited > result.statistics[i - 1].mean_unlimited) ||
                !(stat.mean_final > result.statistics[i - 1].mean_final)) {
                throw std::runtime_error("Wave-height ordering constraint not satisfied.");
            }
        }
    }
    return result;
}

// All displayed and saved reports use the same fixed-width text layout.
static std::wstring report_number(double value, int decimals = 6,
                                  bool scientific = false) {
    std::wostringstream number;
    number.imbue(std::locale::classic());
    number << (scientific ? std::scientific : std::fixed)
           << std::setprecision(decimals) << value;
    return number.str();
}

static std::wstring report_slope(double denominator) {
    std::wstring value = report_number(denominator, 4);
    const std::size_t decimal = value.find(L'.');
    if (decimal != std::wstring::npos) {
        while (!value.empty() && value.back() == L'0') value.pop_back();
        if (!value.empty() && value.back() == L'.') value.pop_back();
    }
    return L"1:" + value;
}

static void report_field(std::wostringstream& out,
                         const std::wstring& label,
                         const std::wstring& value,
                         bool right_align_value = true) {
    constexpr int LABEL_WIDTH = 42;
    constexpr int VALUE_WIDTH = 16;
    out << L"  " << std::left << std::setw(LABEL_WIDTH) << label << L" : ";
    if (right_align_value) {
        out << std::right << std::setw(VALUE_WIDTH) << value;
    } else {
        out << value;
    }
    out << L'\n';
}

static void report_statistic_row(std::wostringstream& out,
                                 const std::wstring& fraction,
                                 const std::wstring& threshold,
                                 const std::wstring& raw_mean,
                                 const std::wstring& rayleigh_mean,
                                 const std::wstring& final_mean,
                                 const std::wstring& limit) {
    out << L"  " << std::left << std::setw(10) << fraction << L' '
        << std::right << std::setw(11) << threshold << L' '
        << std::setw(12) << raw_mean << L' '
        << std::setw(12) << rayleigh_mean << L' '
        << std::setw(12) << final_mean << L' '
        << std::setw(8) << limit << L'\n';
}

static std::wstring format_report(const Analysis& r) {
    std::wostringstream out;
    out.imbue(std::locale::classic());

    out << L"SHALLOW-WATER INDIVIDUAL WAVE-HEIGHT STATISTICS\n"
        << std::wstring(72, L'=') << L"\n\n"
        << L"INPUT PARAMETERS\n";
    report_field(out, L"Spectral significant wave height Hm0 (m)",
                 report_number(r.Hm0, 4));
    report_field(out, L"Still-water depth d (m)",
                 report_number(r.depth, 4));
    report_field(out, L"Foreshore slope",
                 report_slope(r.slope_denominator));
    report_field(out, L"tan(alpha)", report_number(r.slope, 6));

    out << L"\nLOCAL WAVE PARAMETERS\n";
    report_field(out, L"Spectral variance m0 (m^2)",
                 report_number(r.m0, 9, true));
    report_field(out, L"Individual-wave Hrms (m)",
                 report_number(r.Hrms));
    report_field(out, L"Degree of saturation sqrt(m0)/d",
                 report_number(r.saturation));
    report_field(out, L"Transitional wave height Htr (m)",
                 report_number(r.Htr));
    report_field(out, L"Normalized transition Htr/Hrms",
                 report_number(r.transition_normalized));
    report_field(out, L"Distribution",
                 r.uses_cwd ? L"Composite Weibull (B&G)" : L"Rayleigh", false);

    out << L"\nDISTRIBUTION PARAMETERS\n";
    if (r.uses_cwd) {
        report_field(out, L"Normalized first Weibull scale H1/Hrms",
                     report_number(r.scale1_normalized));
        report_field(out, L"Normalized second Weibull scale H2/Hrms",
                     report_number(r.scale2_normalized));
        report_field(out, L"Dimensional Weibull scale H1 (m)",
                     report_number(r.scale1_normalized * r.Hrms));
        report_field(out, L"Dimensional Weibull scale H2 (m)",
                     report_number(r.scale2_normalized * r.Hrms));
        report_field(out, L"Transition exponent q",
                     report_number(r.transition_exponent));
        report_field(out, L"Probability P(H > Htr)",
                     report_number(r.transition_exceedance_probability, 9));
        report_field(out, L"Equivalent transition recurrence Ntr",
                     report_number(r.transition_recurrence_count) + L" waves");
    } else {
        report_field(out, L"Composite Weibull scales H1 and H2", L"N/A", false);
        report_field(out, L"Composite Weibull transition exponent", L"N/A", false);
    }

    out << L"\nWAVE-HEIGHT STATISTICS (m)\n";
    report_statistic_row(out, L"HIGHEST", L"H_N", L"RAW MEAN",
                         L"RAYLEIGH", L"FINAL MEAN", L"LIMIT");
    report_statistic_row(out, L"FRACTION", L"THRESHOLD", L"H_(1/N)",
                         L"H_(1/N)", L"H_(1/N)", L"APPLIED");
    out << L"  " << std::wstring(70, L'-') << L'\n';
    for (const Statistic& s : r.statistics) {
        report_statistic_row(out,
            L"1/" + std::to_wstring(s.N),
            report_number(s.threshold, 4),
            report_number(s.mean_unlimited, 4),
            report_number(s.rayleigh_mean, 4),
            report_number(s.mean_final, 4),
            s.rayleigh_limit_applied ? L"YES" : L"NO");
    }

    out << L"\nFINAL NORMALIZED MEAN-HIGH-WAVE HEIGHTS\n";
    for (const Statistic& s : r.statistics) {
        report_field(out,
                     L"H_(1/" + std::to_wstring(s.N) + L") / Hrms",
                     report_number(s.mean_normalized, 7));
    }

    out << L"\nRATIOS RELATIVE TO THE HIGHEST ONE THIRD\n";
    for (std::size_t i = 1; i < r.statistics.size(); ++i) {
        const Statistic& s = r.statistics[i];
        report_field(out,
                     L"H_(1/" + std::to_wstring(s.N) + L") / H_(1/3)",
                     report_number(s.mean_final / r.statistics[0].mean_final));
    }

    out << L"\nNUMERICAL DIAGNOSTICS\n";
    if (r.uses_cwd) {
        report_field(out, L"Newton iterations",
                     std::to_wstring(r.solver_iterations));
        report_field(out, L"Second-moment residual",
                     report_number(r.moment_residual, 10, true));
        report_field(out, L"Continuity exponent residual",
                     report_number(r.continuity_residual, 10, true));
    } else {
        report_field(out, L"CWD nonlinear solver", L"N/A", false);
    }

    out << L"\nINTERPRETATION\n"
        << L"  H_N is the height exceeded with probability 1/N.\n"
        << L"  H_(1/N) is the mean height of the highest 1/N fraction.\n"
        << L"  Both are individual-wave statistics, not storm maxima.\n"
        << L"  Mean-high-wave Rayleigh limits are applied separately.\n"
        << L"  H_N always refers to the selected underlying distribution.\n"
        << L"  The program uses local wave conditions; it does not propagate waves.\n";

    if (r.slope_outside_calibration ||
        (r.uses_cwd && r.saturation_outside_observed_range)) {
        out << L"\nAPPLICABILITY NOTICE\n";
        if (r.slope_outside_calibration) {
            out << L"  Foreshore slope is outside the laboratory slope range\n"
                << L"  1:20 to 1:250 of the Composite Weibull formulation.\n"
                << L"  Predictions for an extensive horizontal, constant-depth\n"
                << L"  seabed require a separate applicability assessment.\n";
        }
        if (r.uses_cwd && r.saturation_outside_observed_range) {
            out << L"  The degree of saturation is outside the approximate\n"
                << L"  interval represented by the original foreshore tests.\n"
                << L"  Interpretation requires engineering assessment.\n";
        }
    }

    // Neither the GUI nor report.txt includes an empty line after the final text.
    std::wstring report = out.str();
    while (!report.empty() && (report.back() == L'\n' || report.back() == L'\r')) {
        report.pop_back();
    }
    return report;
}

} // namespace cwd

namespace {

struct Options {
    std::vector<std::string> positional;
    std::string output_path = "report.txt";
    bool write_file = true;
    bool show_help = false;
};

static void print_usage(std::ostream& stream, const char* program) {
    stream << "Usage:\n"
           << "  " << program << " Hm0 depth slope_M [--output path] [--no-file]\n"
           << "  " << program << " [--output path] [--no-file]\n"
           << "  " << program << " --help\n\n"
           << "Hm0       Local spectral significant wave height (m).\n"
           << "depth     Local still-water depth (m).\n"
           << "slope_M   Positive denominator M of a 1:M foreshore slope.\n\n"
           << "With no numeric arguments, the program prompts for inputs.\n"
           << "The report is displayed on the console and saved as report.txt\n"
           << "unless --no-file is specified. --output selects another file.\n";
}

static Options parse_options(int argc, char* argv[]) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") {
            options.show_help = true;
        } else if (argument == "--output" || argument == "-o") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing path after --output.");
            }
            options.output_path = argv[++i];
            if (options.output_path.empty()) {
                throw std::invalid_argument("Output path must not be empty.");
            }
        } else if (argument == "--no-file") {
            options.write_file = false;
        } else if (argument.size() >= 2 && argument[0] == '-' && argument[1] == '-') {
            throw std::invalid_argument("Unrecognized option: " + argument);
        } else {
            options.positional.push_back(argument);
        }
    }
    if (!options.show_help && !options.positional.empty() &&
        options.positional.size() != 3) {
        throw std::invalid_argument("Provide exactly three numerical inputs or none.");
    }
    return options;
}

static double parse_positive_number(std::string raw, const char* name) {
    std::replace(raw.begin(), raw.end(), ',', '.');
    std::istringstream input(raw);
    input.imbue(std::locale::classic());
    double value = 0.0;
    if (!(input >> value)) {
        throw std::invalid_argument(std::string("Invalid numeric value for ") + name + ".");
    }
    input >> std::ws;
    if (!input.eof() || !std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string("Invalid numeric value for ") + name + ".");
    }
    return value;
}

static double prompt_positive_number(const char* prompt, const char* name) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!std::getline(std::cin, line)) {
        throw std::runtime_error(std::string("No input available for ") + name + ".");
    }
    return parse_positive_number(line, name);
}

// Report text is ASCII so the terminal and saved file contain exactly the
// same content as the wide-character report used by the native GUI.
static std::string report_as_utf8(const std::wstring& report) {
    std::string result;
    result.reserve(report.size());
    for (wchar_t character : report) {
        if (static_cast<unsigned int>(character) > 0x7fu) {
            throw std::runtime_error("The report contains an unsupported character.");
        }
        result.push_back(static_cast<char>(character));
    }
    return result;
}

static void write_report(const std::string& path, const std::string& content) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        throw std::runtime_error("Cannot open output file: " + path);
    }
    file.write(content.data(), static_cast<std::streamsize>(content.size()));
    file.flush();
    if (!file) {
        throw std::runtime_error("Cannot write report completely: " + path);
    }
    file.close();
    if (!file) {
        throw std::runtime_error("Cannot close output file: " + path);
    }
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.show_help) {
            print_usage(std::cout, argv[0]);
            return 0;
        }

        double Hm0 = 0.0;
        double depth = 0.0;
        double slope_denominator = 0.0;
        if (options.positional.empty()) {
            Hm0 = prompt_positive_number("Hm0 (m): ", "Hm0");
            depth = prompt_positive_number("Water depth d (m): ", "water depth");
            slope_denominator = prompt_positive_number("Foreshore slope 1:M, enter M: ",
                                                       "slope denominator");
            std::cout << '\n';
        } else {
            Hm0 = parse_positive_number(options.positional[0], "Hm0");
            depth = parse_positive_number(options.positional[1], "water depth");
            slope_denominator = parse_positive_number(options.positional[2],
                                                      "slope denominator");
        }

        const cwd::Analysis data = cwd::calculate(Hm0, depth, slope_denominator);
        const std::string report = report_as_utf8(cwd::format_report(data));
        if (options.write_file) {
            write_report(options.output_path, report);
        }
        std::cout << report << '\n';
        if (!std::cout) {
            throw std::runtime_error("Unable to write report to standard output.");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
