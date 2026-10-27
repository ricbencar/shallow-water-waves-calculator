/*
 * Shallow-Water Individual Wave-Height Statistics
 *
 * Native Windows application for local individual zero-crossing wave-height
 * statistics based on the Composite Weibull Distribution of Battjes and
 * Groenendijk. The Rayleigh distribution is used for weakly depth-limited
 * conditions, with independent Rayleigh limits on the CWD mean-high-wave
 * statistics.
 *
 * Inputs: local spectral Hm0 (m), local still-water depth d (m), and
 * foreshore slope expressed as 1:M.
 *
 * Build with MinGW-w64:
 * g++ -O2 -std=c++17 -Wall -Wextra -pedantic -municode -mwindows \
 *     shallow-water-waves_gui.cpp -o shallow-water-waves_gui.exe \
 *     -static -static-libgcc -static-libstdc++
 */

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstring>
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
                 report_number(r.m0));
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

#ifdef _WIN32

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#ifndef WM_DPICHANGED
#define WM_DPICHANGED 0x02E0
#endif

namespace {

constexpr int IDC_EDIT_HM0 = 1001;
constexpr int IDC_EDIT_DEPTH = 1002;
constexpr int IDC_EDIT_SLOPE = 1003;
constexpr int IDC_COMPUTE = 1004;
constexpr int IDC_OUTPUT = 1005;

HWND hLabelHm0 = nullptr;
HWND hLabelDepth = nullptr;
HWND hLabelSlope = nullptr;
HWND hEditHm0 = nullptr;
HWND hEditDepth = nullptr;
HWND hEditSlope = nullptr;
HWND hCompute = nullptr;
HWND hOutput = nullptr;
HWND hStatus = nullptr;
HFONT hReportFont = nullptr;
HFONT hUiFont = nullptr;
UINT gDpi = 96;

static int dip(int logical_pixels) {
    return MulDiv(logical_pixels, static_cast<int>(gDpi), 96);
}

// Load optional Windows API entry points without incompatible function-pointer
// casts. On supported Win32/Win64 targets, FARPROC and WINAPI function pointers
// have the same representation.
template <typename FunctionPointer>
static FunctionPointer optional_winapi(HMODULE module, const char* name) {
    static_assert(sizeof(FunctionPointer) == sizeof(FARPROC),
                  "Unsupported Windows function-pointer representation.");
    const FARPROC address = module ? GetProcAddress(module, name) : nullptr;
    FunctionPointer function = nullptr;
    std::memcpy(&function, &address, sizeof(function));
    return function;
}

// Enable per-monitor DPI scaling when supported by Windows. Both APIs are
// resolved at runtime so the executable also loads on older Windows systems.
static void enable_dpi_awareness() {
    HMODULE user32 = GetModuleHandleW(L"user32.dll");
    if (!user32) return;
    using SetDpiContext = BOOL (WINAPI*)(HANDLE);
    const auto set_context = optional_winapi<SetDpiContext>(
        user32, "SetProcessDpiAwarenessContext");
    if (set_context) {
        const HANDLE per_monitor_v2 = reinterpret_cast<HANDLE>(static_cast<INT_PTR>(-4));
        if (set_context(per_monitor_v2)) return;
        const HANDLE per_monitor_v1 = reinterpret_cast<HANDLE>(static_cast<INT_PTR>(-3));
        if (set_context(per_monitor_v1)) return;
    }
    using SetDpiAware = BOOL (WINAPI*)();
    const auto set_aware = optional_winapi<SetDpiAware>(
        user32, "SetProcessDPIAware");
    if (set_aware) set_aware();
}

static UINT system_dpi() {
    HMODULE user32 = GetModuleHandleW(L"user32.dll");
    if (user32) {
        using GetSystemDpi = UINT (WINAPI*)();
        const auto get_dpi = optional_winapi<GetSystemDpi>(
            user32, "GetDpiForSystem");
        if (get_dpi) {
            const UINT value = get_dpi();
            if (value >= 48 && value <= 768) return value;
        }
    }
    HDC desktop = GetDC(nullptr);
    const int dpi = desktop ? GetDeviceCaps(desktop, LOGPIXELSX) : 96;
    if (desktop) ReleaseDC(nullptr, desktop);
    return (dpi >= 48 && dpi <= 768) ? static_cast<UINT>(dpi) : 96u;
}

static UINT window_dpi(HWND window) {
    HMODULE user32 = GetModuleHandleW(L"user32.dll");
    if (user32) {
        using GetWindowDpi = UINT (WINAPI*)(HWND);
        const auto get_dpi = optional_winapi<GetWindowDpi>(
            user32, "GetDpiForWindow");
        if (get_dpi) {
            const UINT value = get_dpi(window);
            if (value >= 48 && value <= 768) return value;
        }
    }
    return system_dpi();
}

static void send_font(HWND control, HFONT font) {
    if (control && font) {
        SendMessageW(control, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    }
}

static void create_fonts(UINT dpi) {
    const UINT bounded_dpi = std::max(48u, std::min(dpi, 768u));
    const int ui_height = -MulDiv(17, static_cast<int>(bounded_dpi), 96);
    const int report_height = -MulDiv(19, static_cast<int>(bounded_dpi), 96);
    HFONT new_ui = CreateFontW(ui_height, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                               DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                               CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                               DEFAULT_PITCH | FF_DONTCARE, L"Segoe UI");
    HFONT new_report = CreateFontW(report_height, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                                   DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                                   CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                   FIXED_PITCH | FF_MODERN, L"Consolas");
    if (!new_ui || !new_report) {
        if (new_ui) DeleteObject(new_ui);
        if (new_report) DeleteObject(new_report);
        throw std::runtime_error("Unable to create interface fonts.");
    }
    HFONT old_ui = hUiFont;
    HFONT old_report = hReportFont;
    hUiFont = new_ui;
    hReportFont = new_report;
    gDpi = bounded_dpi;
    send_font(hLabelHm0, hUiFont);
    send_font(hLabelDepth, hUiFont);
    send_font(hLabelSlope, hUiFont);
    send_font(hEditHm0, hUiFont);
    send_font(hEditDepth, hUiFont);
    send_font(hEditSlope, hUiFont);
    send_font(hCompute, hUiFont);
    send_font(hStatus, hUiFont);
    send_font(hOutput, hReportFont);
    if (old_ui) DeleteObject(old_ui);
    if (old_report) DeleteObject(old_report);
}

static std::wstring error_to_wide(const std::exception& e) {
    const std::string message(e.what());
    return std::wstring(message.begin(), message.end());
}

static std::wstring read_text(HWND edit) {
    const int length = GetWindowTextLengthW(edit);
    if (length <= 0 || length > 200) {
        throw std::invalid_argument("An input field is empty or too long.");
    }
    std::vector<wchar_t> buffer(static_cast<std::size_t>(length) + 1u, L'\0');
    GetWindowTextW(edit, buffer.data(), length + 1);
    return std::wstring(buffer.data());
}

static double parse_positive_input(HWND edit, const char* label) {
    std::wstring input = read_text(edit);
    std::replace(input.begin(), input.end(), L',', L'.');
    errno = 0;
    wchar_t* end = nullptr;
    const double value = std::wcstod(input.c_str(), &end);
    if (end == input.c_str() || errno == ERANGE) {
        throw std::invalid_argument(std::string("Invalid numeric value for ") + label + ".");
    }
    while (*end != L'\0' && std::iswspace(static_cast<wint_t>(*end))) ++end;
    if (*end != L'\0' || !std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string("Invalid numeric value for ") + label + ".");
    }
    return value;
}

static std::string utf8_from_wide(const std::wstring& content) {
    if (content.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Report is too large for UTF-8 encoding.");
    }
    const int count = static_cast<int>(content.size());
    const int bytes = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                                          content.data(), count, nullptr, 0,
                                          nullptr, nullptr);
    if (bytes <= 0) throw std::runtime_error("Unable to encode report as UTF-8.");
    std::string output(static_cast<std::size_t>(bytes), '\0');
    const int written = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                                            content.data(), count, &output[0], bytes,
                                            nullptr, nullptr);
    if (written != bytes) throw std::runtime_error("Incomplete UTF-8 encoding.");
    return output;
}

static void write_report(const std::wstring& report) {
    const std::string utf8 = utf8_from_wide(report);
    HANDLE file = CreateFileW(L"report.txt", GENERIC_WRITE, FILE_SHARE_READ,
                              nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE) {
        throw std::runtime_error("Could not open report.txt for writing.");
    }
    if (utf8.size() > static_cast<std::size_t>(std::numeric_limits<DWORD>::max())) {
        CloseHandle(file);
        throw std::runtime_error("The report is too large to write.");
    }
    DWORD written = 0;
    const BOOL success = WriteFile(file, utf8.data(),
                                   static_cast<DWORD>(utf8.size()), &written, nullptr);
    const BOOL closed = CloseHandle(file);
    if (!success || written != static_cast<DWORD>(utf8.size()) || !closed) {
        throw std::runtime_error("Could not write report.txt completely.");
    }
}

static std::wstring to_windows_newlines(const std::wstring& text) {
    std::wstring result;
    result.reserve(text.size() + 128);
    for (wchar_t ch : text) {
        if (ch == L'\n') result.push_back(L'\r');
        result.push_back(ch);
    }
    return result;
}

static HWND add_control(DWORD style, const wchar_t* type, const wchar_t* text,
                        HWND parent, int id = 0, DWORD ex_style = 0) {
    HWND control = CreateWindowExW(ex_style, type, text, style,
                                   0, 0, 1, 1, parent,
                                   reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
                                   GetModuleHandleW(nullptr), nullptr);
    if (!control) throw std::runtime_error("Could not initialize window controls.");
    send_font(control, hUiFont);
    return control;
}

static void calculate_and_display(HWND window) {
    try {
        const double Hm0 = parse_positive_input(hEditHm0, "Hm0");
        const double depth = parse_positive_input(hEditDepth, "water depth");
        const double M = parse_positive_input(hEditSlope, "foreshore slope denominator");
        const cwd::Analysis data = cwd::calculate(Hm0, depth, M);
        const std::wstring report = cwd::format_report(data);
        const std::wstring window_report = to_windows_newlines(report);
        SetWindowTextW(hOutput, window_report.c_str());
        SendMessageW(hOutput, EM_SETSEL, 0, 0);
        SendMessageW(hOutput, EM_SCROLLCARET, 0, 0);
        SetWindowTextW(hStatus, L"Calculation completed.");
        try {
            write_report(report);
            SetWindowTextW(hStatus, L"Calculation completed. report.txt saved to working directory.");
        } catch (const std::exception& e) {
            const std::wstring message = error_to_wide(e);
            SetWindowTextW(hStatus, L"Calculation completed. Report file not saved.");
            MessageBoxW(window, message.c_str(), L"Report File Error", MB_OK | MB_ICONWARNING);
        }
    } catch (const std::exception& e) {
        const std::wstring message = error_to_wide(e);
        SetWindowTextW(hOutput, L"");
        SetWindowTextW(hStatus, L"Calculation failed.");
        MessageBoxW(window, message.c_str(), L"Calculation Error", MB_OK | MB_ICONERROR);
    }
}

static void move_control(HWND control, int x, int y, int width, int height) {
    if (control) {
        SetWindowPos(control, nullptr, x, y, std::max(1, width), std::max(1, height),
                     SWP_NOZORDER | SWP_NOACTIVATE);
    }
}

// The toolbar places labels above the corresponding edit controls. The
// Calculate button shares an input row when space permits. At smaller client
// widths the controls occupy two rows so nothing extends beyond the window.
static void arrange_controls(HWND window) {
    if (!hLabelHm0 || !hLabelDepth || !hLabelSlope || !hEditHm0 || !hEditDepth ||
        !hEditSlope || !hCompute || !hStatus || !hOutput) return;
    RECT client{};
    GetClientRect(window, &client);
    const int client_width = client.right - client.left;
    const int client_height = client.bottom - client.top;
    if (client_width <= 0 || client_height <= 0) return;

    const int margin = dip(18);
    const int gap = dip(14);
    const int top = dip(16);
    const int label_height = dip(25);
    const int label_gap = dip(4);
    const int input_height = dip(38);
    const int button_height = dip(40);
    const int button_width = dip(150);
    const int status_height = dip(26);
    const int available = std::max(dip(100), client_width - 2 * margin);
    const int edit_top = top + label_height + label_gap;
    const int row_bottom = edit_top + input_height;

    int status_y = 0;
    int status_x = margin;
    int status_width = available;
    auto place_field = [&](HWND label, HWND edit, int x, int y, int w) {
        move_control(label, x, y, w, label_height);
        move_control(edit, x, y + label_height + label_gap, w, input_height);
    };

    if (available >= 3 * dip(185) + 3 * gap + button_width) {
        // Three input groups and the action button in one row.
        const int field_width = (available - button_width - 3 * gap) / 3;
        const int x1 = margin;
        const int x2 = x1 + field_width + gap;
        const int x3 = x2 + field_width + gap;
        place_field(hLabelHm0, hEditHm0, x1, top, field_width);
        place_field(hLabelDepth, hEditDepth, x2, top, field_width);
        place_field(hLabelSlope, hEditSlope, x3, top, field_width);
        move_control(hCompute, margin + available - button_width,
                     edit_top - dip(1), button_width, button_height);
        status_y = row_bottom + dip(15);
    } else if (available >= 3 * dip(170) + 2 * gap) {
        // Three inputs above an action/status row.
        const int field_width = (available - 2 * gap) / 3;
        place_field(hLabelHm0, hEditHm0, margin, top, field_width);
        place_field(hLabelDepth, hEditDepth, margin + field_width + gap,
                    top, field_width);
        place_field(hLabelSlope, hEditSlope, margin + 2 * (field_width + gap),
                    top, field_width);
        const int button_y = row_bottom + dip(13);
        move_control(hCompute, margin + available - button_width, button_y,
                     button_width, button_height);
        status_y = button_y + (button_height - status_height) / 2;
        status_width = available - button_width - gap;
    } else if (available >= 2 * dip(140) + gap) {
        // Compact layout: Hm0 and depth above slope and Calculate.
        const int field_width = (available - gap) / 2;
        place_field(hLabelHm0, hEditHm0, margin, top, field_width);
        place_field(hLabelDepth, hEditDepth, margin + field_width + gap,
                    top, field_width);
        const int second_row = row_bottom + dip(13);
        place_field(hLabelSlope, hEditSlope, margin, second_row, field_width);
        move_control(hCompute, margin + field_width + gap,
                     second_row + label_height + label_gap - dip(1),
                     field_width, button_height);
        status_y = second_row + label_height + label_gap + button_height + dip(12);
    } else {
        // The single-column arrangement accommodates very small displays.
        place_field(hLabelHm0, hEditHm0, margin, top, available);
        const int second_row = row_bottom + dip(10);
        place_field(hLabelDepth, hEditDepth, margin, second_row, available);
        const int third_row = second_row + label_height + label_gap + input_height + dip(10);
        place_field(hLabelSlope, hEditSlope, margin, third_row, available);
        const int button_y = third_row + label_height + label_gap + input_height + dip(12);
        move_control(hCompute, margin, button_y, available, button_height);
        status_y = button_y + button_height + dip(12);
    }

    move_control(hStatus, status_x, status_y, status_width, status_height);
    const int report_top = (status_width == available ? status_y + status_height + dip(12) :
                            row_bottom + dip(13) + button_height + dip(13));
    move_control(hOutput, margin, report_top, available,
                 std::max(1, client_height - report_top - margin));
}

static LRESULT CALLBACK window_proc(HWND window, UINT message,
                                    WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_CREATE:
        try {
            create_fonts(window_dpi(window));
            hLabelHm0 = add_control(WS_CHILD | WS_VISIBLE | SS_LEFT,
                                    L"STATIC", L"Hm0 (m)", window);
            hEditHm0 = add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL,
                                   L"EDIT", L"2.5", window, IDC_EDIT_HM0, WS_EX_CLIENTEDGE);
            hLabelDepth = add_control(WS_CHILD | WS_VISIBLE | SS_LEFT,
                                      L"STATIC", L"Depth d (m)", window);
            hEditDepth = add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL,
                                     L"EDIT", L"5.0", window, IDC_EDIT_DEPTH, WS_EX_CLIENTEDGE);
            hLabelSlope = add_control(WS_CHILD | WS_VISIBLE | SS_LEFT,
                                      L"STATIC", L"Slope 1:M", window);
            hEditSlope = add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL,
                                     L"EDIT", L"100", window, IDC_EDIT_SLOPE, WS_EX_CLIENTEDGE);
            hCompute = add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | BS_DEFPUSHBUTTON,
                                   L"BUTTON", L"Calculate", window, IDC_COMPUTE);
            hStatus = add_control(WS_CHILD | WS_VISIBLE | SS_LEFT | SS_ENDELLIPSIS,
                                  L"STATIC", L"Ready.", window);
            hOutput = add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_VSCROLL |
                                   WS_HSCROLL | ES_MULTILINE | ES_AUTOVSCROLL |
                                   ES_AUTOHSCROLL | ES_READONLY,
                                   L"EDIT", L"", window, IDC_OUTPUT, WS_EX_CLIENTEDGE);
            send_font(hOutput, hReportFont);
            const int edit_margin = dip(7);
            const LPARAM margins = MAKELPARAM(edit_margin, edit_margin);
            SendMessageW(hEditHm0, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            SendMessageW(hEditDepth, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            SendMessageW(hEditSlope, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            arrange_controls(window);
        } catch (const std::exception& e) {
            const std::wstring error = error_to_wide(e);
            MessageBoxW(window, error.c_str(), L"Startup Error", MB_OK | MB_ICONERROR);
            return -1;
        }
        return 0;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_COMPUTE && HIWORD(wParam) == BN_CLICKED) {
            calculate_and_display(window);
        }
        return 0;

    case WM_SIZE:
        if (wParam != SIZE_MINIMIZED) arrange_controls(window);
        return 0;

    case WM_DPICHANGED: {
        const RECT* suggested = reinterpret_cast<const RECT*>(lParam);
        try {
            create_fonts(LOWORD(wParam));
            if (suggested) {
                SetWindowPos(window, nullptr, suggested->left, suggested->top,
                             suggested->right - suggested->left,
                             suggested->bottom - suggested->top,
                             SWP_NOZORDER | SWP_NOACTIVATE);
            }
            const int edit_margin = dip(7);
            const LPARAM margins = MAKELPARAM(edit_margin, edit_margin);
            SendMessageW(hEditHm0, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            SendMessageW(hEditDepth, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            SendMessageW(hEditSlope, EM_SETMARGINS, EC_LEFTMARGIN | EC_RIGHTMARGIN, margins);
            arrange_controls(window);
        } catch (const std::exception& e) {
            const std::wstring error = error_to_wide(e);
            MessageBoxW(window, error.c_str(), L"Display Scaling Error", MB_OK | MB_ICONWARNING);
        }
        return 0;
    }

    case WM_GETMINMAXINFO: {
        const UINT dpi = window_dpi(window);
        const int minimum_width = MulDiv(510, static_cast<int>(dpi), 96);
        const int minimum_height = MulDiv(460, static_cast<int>(dpi), 96);
        int allowed_width = minimum_width;
        int allowed_height = minimum_height;
        MONITORINFO monitor{};
        monitor.cbSize = sizeof(monitor);
        const HMONITOR current = MonitorFromWindow(window, MONITOR_DEFAULTTONEAREST);
        if (current && GetMonitorInfoW(current, &monitor)) {
            const int work_width = static_cast<int>(monitor.rcWork.right - monitor.rcWork.left);
            const int work_height = static_cast<int>(monitor.rcWork.bottom - monitor.rcWork.top);
            allowed_width = std::max(320, work_width);
            allowed_height = std::max(320, work_height);
        }
        MINMAXINFO* info = reinterpret_cast<MINMAXINFO*>(lParam);
        info->ptMinTrackSize.x = std::min(minimum_width, allowed_width);
        info->ptMinTrackSize.y = std::min(minimum_height, allowed_height);
        return 0;
    }

    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;

    case WM_NCDESTROY:
        if (hReportFont) { DeleteObject(hReportFont); hReportFont = nullptr; }
        if (hUiFont) { DeleteObject(hUiFont); hUiFont = nullptr; }
        return DefWindowProcW(window, message, wParam, lParam);

    default:
        return DefWindowProcW(window, message, wParam, lParam);
    }
}

} // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show_command) {
    enable_dpi_awareness();
    constexpr wchar_t CLASS_NAME[] = L"ShallowForeshoreWaveStatistics";
    WNDCLASSEXW wc{};
    wc.cbSize = static_cast<UINT>(sizeof(wc));
    wc.lpfnWndProc = window_proc;
    wc.hInstance = instance;
    wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    wc.hIcon = LoadIconW(nullptr, IDI_APPLICATION);
    wc.hIconSm = LoadIconW(nullptr, IDI_APPLICATION);
    wc.hbrBackground = reinterpret_cast<HBRUSH>(COLOR_BTNFACE + 1);
    wc.lpszClassName = CLASS_NAME;

    if (!RegisterClassExW(&wc)) {
        MessageBoxW(nullptr, L"Window registration failed.", L"Startup Error",
                    MB_OK | MB_ICONERROR);
        return 1;
    }

    // Select a comfortable default size without extending outside the
    // available desktop work area (including at high display scaling).
    const UINT startup_dpi = system_dpi();
    const int scale = static_cast<int>(startup_dpi);
    RECT work{};
    if (!SystemParametersInfoW(SPI_GETWORKAREA, 0, &work, 0)) {
        work = {0, 0, GetSystemMetrics(SM_CXSCREEN), GetSystemMetrics(SM_CYSCREEN)};
    }
    const int available_width = std::max(320, static_cast<int>(work.right - work.left));
    const int available_height = std::max(320, static_cast<int>(work.bottom - work.top));
    const int preferred_width = MulDiv(990, scale, 96);
    const int preferred_height = MulDiv(730, scale, 96);
    const int inset = MulDiv(20, scale, 96);
    const int width = std::min(preferred_width, std::max(320, available_width - inset));
    const int height = std::min(preferred_height, std::max(320, available_height - inset));
    const int x = work.left + (available_width - width) / 2;
    const int y = work.top + (available_height - height) / 2;

    HWND window = CreateWindowExW(0, CLASS_NAME,
                                  L"Shallow-Water Individual Wave-Height Statistics",
                                  WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
                                  x, y, width, height,
                                  nullptr, nullptr, instance, nullptr);
    if (!window) {
        MessageBoxW(nullptr, L"Window creation failed.", L"Startup Error",
                    MB_OK | MB_ICONERROR);
        return 1;
    }
    ShowWindow(window, show_command);
    UpdateWindow(window);
    if (hEditHm0) {
        SetFocus(hEditHm0);
        SendMessageW(hEditHm0, EM_SETSEL, 0, -1);
    }

    MSG msg{};
    int message_status = 0;
    while ((message_status = GetMessageW(&msg, nullptr, 0, 0)) > 0) {
        if (msg.message == WM_KEYDOWN && msg.wParam == VK_RETURN &&
            (msg.hwnd == hEditHm0 || msg.hwnd == hEditDepth ||
             msg.hwnd == hEditSlope || msg.hwnd == hCompute)) {
            calculate_and_display(window);
            continue;
        }
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    if (message_status < 0) return 1;
    return static_cast<int>(msg.wParam);
}

#endif // _WIN32
