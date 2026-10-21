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
constexpr int MAX_ROOT_ITERATIONS = 200;

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

// For q=(Htr/H1)^k1=(Htr/H2)^k2, the normalized second moment is
// t^2 * [q^(-2/k1)*gamma(1+2/k1,q)
//       +q^(-2/k2)*Gamma(1+2/k2,q)], where t=Htr/Hrms.
static double logarithmic_moment_residual(double log_q, double log_t) {
    const double q = std::exp(log_q);
    const double lower = std::exp(-(2.0 / K1) * log_q) *
                         gamma_lower(1.0 + 2.0 / K1, q);
    const double upper = std::exp(-(2.0 / K2) * log_q) *
                         gamma_upper(1.0 + 2.0 / K2, q);
    const double bracket = lower + upper;
    require_finite_positive(bracket, "Second-moment kernel");
    return 2.0 * log_t + std::log(bracket);
}

// Brent's bracketed scalar solver, operating on z=log(q) to keep q positive.
// The moment function is strictly decreasing with q.
static double solve_transition_exponent(double t, int& iterations) {
    require_finite_positive(t, "Normalized transition height");
    const double log_t = std::log(t);
    const auto f = [log_t](double z) {
        return logarithmic_moment_residual(z, log_t);
    };

    double a = -700.0;
    double b = std::log(64.0);
    double fa = f(a);
    double fb = f(b);
    if (!(std::isfinite(fa) && std::isfinite(fb) && fa > 0.0 && fb < 0.0)) {
        throw std::runtime_error("CWD scalar solution is outside the representable range.");
    }

    double c = a;
    double fc = fa;
    double step = b - a;
    double previous_step = step;

    for (int iteration = 1; iteration <= MAX_ROOT_ITERATIONS; ++iteration) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            step = b - a;
            previous_step = step;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        const double tolerance = 2.0 * std::numeric_limits<double>::epsilon() *
                                 std::abs(b) + 2.0e-14;
        const double mid = 0.5 * (c - b);
        if (std::abs(mid) <= tolerance || fb == 0.0) {
            iterations = iteration;
            const double q = std::exp(b);
            require_finite_positive(q, "CWD transition exponent");
            if (std::abs(f(b)) > 1e-10) {
                throw std::runtime_error("CWD root residual exceeds tolerance.");
            }
            return q;
        }

        if (std::abs(previous_step) >= tolerance && std::abs(fa) > std::abs(fb)) {
            const double s = fb / fa;
            double p = 0.0;
            double denom = 0.0;
            if (a == c) {
                p = 2.0 * mid * s;
                denom = 1.0 - s;
            } else {
                const double q_ratio = fa / fc;
                const double r_ratio = fb / fc;
                p = s * (2.0 * mid * q_ratio * (q_ratio - r_ratio) -
                         (b - a) * (r_ratio - 1.0));
                denom = (q_ratio - 1.0) * (r_ratio - 1.0) * (s - 1.0);
            }
            if (p > 0.0) denom = -denom;
            p = std::abs(p);
            const double bound1 = 3.0 * mid * denom - std::abs(tolerance * denom);
            const double bound2 = std::abs(previous_step * denom);
            if (denom != 0.0 && 2.0 * p < std::min(bound1, bound2)) {
                previous_step = step;
                step = p / denom;
            } else {
                step = mid;
                previous_step = mid;
            }
        } else {
            step = mid;
            previous_step = mid;
        }

        a = b;
        fa = fb;
        b += (std::abs(step) > tolerance) ? step : std::copysign(tolerance, mid);
        fb = f(b);
        if (!std::isfinite(fb)) {
            throw std::runtime_error("Non-finite scalar root residual.");
        }
    }
    throw std::runtime_error("CWD scalar solver did not converge.");
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
        report_field(out, L"Scalar root iterations",
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

namespace {

constexpr int IDC_EDIT_HM0 = 1001;
constexpr int IDC_EDIT_DEPTH = 1002;
constexpr int IDC_EDIT_SLOPE = 1003;
constexpr int IDC_COMPUTE = 1004;
constexpr int IDC_OUTPUT = 1005;

HWND hEditHm0 = nullptr;
HWND hEditDepth = nullptr;
HWND hEditSlope = nullptr;
HWND hOutput = nullptr;
HWND hStatus = nullptr;
HFONT hReportFont = nullptr;
HFONT hUiFont = nullptr;

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
    // Accept decimal points and decimal commas, without accepting partial input.
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
                        int x, int y, int width, int height, HWND parent, int id = 0,
                        DWORD ex_style = 0) {
    HWND window = CreateWindowExW(ex_style, type, text, style,
                                  x, y, width, height, parent,
                                  reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
                                  GetModuleHandleW(nullptr), nullptr);
    if (!window) throw std::runtime_error("Could not initialize window controls.");
    if (hUiFont) SendMessageW(window, WM_SETFONT, reinterpret_cast<WPARAM>(hUiFont), TRUE);
    return window;
}

static void calculate_and_display(HWND window) {
    try {
        const double Hm0 = parse_positive_input(hEditHm0, "Hm0");
        const double depth = parse_positive_input(hEditDepth, "water depth");
        const double M = parse_positive_input(hEditSlope, "foreshore slope denominator");
        const cwd::Analysis data = cwd::calculate(Hm0, depth, M);
        std::wstring report = cwd::format_report(data);
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

static void arrange_controls(HWND window, int client_width, int client_height) {
    (void)window;
    if (!hOutput || !hStatus) return;
    const int margin = 16;
    const int report_top = 98;
    const int report_width = std::max(150, client_width - 2 * margin);
    const int report_height = std::max(150, client_height - report_top - margin);
    MoveWindow(hOutput, margin, report_top, report_width, report_height, TRUE);
    MoveWindow(hStatus, margin, 70, report_width, 23, TRUE);
}

static LRESULT CALLBACK window_proc(HWND window, UINT message,
                                    WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_CREATE:
        try {
            hUiFont = CreateFontW(-17, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                                  DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                                  CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                  DEFAULT_PITCH | FF_DONTCARE, L"Segoe UI");
            hReportFont = CreateFontW(-19, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                                      DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                                      CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                      FIXED_PITCH | FF_MODERN, L"Consolas");
            add_control(WS_CHILD | WS_VISIBLE, L"STATIC", L"Hm0 (m):",
                        16, 17, 88, 23, window);
            hEditHm0 = add_control(WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP |
                                   ES_AUTOHSCROLL, L"EDIT", L"2.5",
                                   106, 15, 104, 26, window, IDC_EDIT_HM0);
            add_control(WS_CHILD | WS_VISIBLE, L"STATIC", L"Depth d (m):",
                        228, 17, 110, 23, window);
            hEditDepth = add_control(WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP |
                                     ES_AUTOHSCROLL, L"EDIT", L"5.0",
                                     336, 15, 104, 26, window, IDC_EDIT_DEPTH);
            add_control(WS_CHILD | WS_VISIBLE, L"STATIC", L"Slope 1:M:",
                        458, 17, 101, 23, window);
            hEditSlope = add_control(WS_CHILD | WS_VISIBLE | WS_BORDER | WS_TABSTOP |
                                     ES_AUTOHSCROLL, L"EDIT", L"100",
                                     554, 15, 104, 26, window, IDC_EDIT_SLOPE);
            add_control(WS_CHILD | WS_VISIBLE | WS_TABSTOP | BS_DEFPUSHBUTTON,
                        L"BUTTON", L"Compute", 681, 13, 155, 31, window, IDC_COMPUTE);
            hStatus = add_control(WS_CHILD | WS_VISIBLE, L"STATIC", L"Ready.",
                                  16, 70, 850, 24, window);
            hOutput = add_control(WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL |
                                   ES_MULTILINE | ES_AUTOVSCROLL | ES_AUTOHSCROLL |
                                   ES_READONLY, L"EDIT", L"",
                                   16, 98, 920, 620, window, IDC_OUTPUT,
                                   WS_EX_CLIENTEDGE);
            if (hReportFont) SendMessageW(hOutput, WM_SETFONT,
                                          reinterpret_cast<WPARAM>(hReportFont), TRUE);
            RECT client{};
            GetClientRect(window, &client);
            arrange_controls(window, static_cast<int>(client.right),
                             static_cast<int>(client.bottom));
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
        arrange_controls(window, LOWORD(lParam), HIWORD(lParam));
        return 0;

    case WM_GETMINMAXINFO: {
        MINMAXINFO* minmax = reinterpret_cast<MINMAXINFO*>(lParam);
        minmax->ptMinTrackSize.x = 900;
        minmax->ptMinTrackSize.y = 510;
        return 0;
    }

    case WM_DESTROY:
        if (hReportFont) { DeleteObject(hReportFont); hReportFont = nullptr; }
        if (hUiFont) { DeleteObject(hUiFont); hUiFont = nullptr; }
        PostQuitMessage(0);
        return 0;

    default:
        return DefWindowProcW(window, message, wParam, lParam);
    }
}

} // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show_command) {
    constexpr wchar_t CLASS_NAME[] = L"ShallowForeshoreWaveStatistics";
    WNDCLASSEXW wc{};
    wc.cbSize = static_cast<UINT>(sizeof(wc));
    wc.lpfnWndProc = window_proc;
    wc.hInstance = instance;
    wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    wc.hIcon = LoadIconW(nullptr, IDI_APPLICATION);
    wc.hIconSm = LoadIconW(nullptr, IDI_APPLICATION);
    wc.hbrBackground = reinterpret_cast<HBRUSH>(COLOR_WINDOW + 1);
    wc.lpszClassName = CLASS_NAME;

    if (!RegisterClassExW(&wc)) {
        MessageBoxW(nullptr, L"Window registration failed.", L"Startup Error",
                    MB_OK | MB_ICONERROR);
        return 1;
    }

    HWND window = CreateWindowExW(0, CLASS_NAME,
                                  L"Shallow-Water Individual Wave-Height Statistics",
                                  WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                                  930, 700, nullptr, nullptr, instance, nullptr);
    if (!window) {
        MessageBoxW(nullptr, L"Window creation failed.", L"Startup Error",
                    MB_OK | MB_ICONERROR);
        return 1;
    }
    ShowWindow(window, show_command);
    UpdateWindow(window);

    MSG msg{};
    int message_status = 0;
    while ((message_status = GetMessageW(&msg, nullptr, 0, 0)) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    if (message_status < 0) return 1;
    return static_cast<int>(msg.wParam);
}

#endif // _WIN32
