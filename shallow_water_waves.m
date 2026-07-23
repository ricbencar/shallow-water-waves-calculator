function results = shallow_water_waves(HM0, DEPTH, M)
%
% Syntax:
%   results = shallow_water_waves(HM0, DEPTH, M)
%
% Inputs:
%   HM0   - significant spectral wave height, Hm0 [m]
%   DEPTH - local water depth, d [m]
%   M     - foreshore slope denominator for a 1:M slope [-]
%
% Output:
%   results - structure containing the same calculated quantities as the
%             C++ WaveAnalysisResults struct.
%

    if nargin ~= 3
        error('shallow_water_waves requires exactly 3 inputs: HM0, DEPTH, M.');
    end
    if ~isnumeric(HM0) || prod(size(HM0)) ~= 1 || ~is_finite_value(HM0) || HM0 <= 0
        error('HM0 must be a positive finite scalar.');
    end
    if ~isnumeric(DEPTH) || prod(size(DEPTH)) ~= 1 || ~is_finite_value(DEPTH) || DEPTH <= 0
        error('DEPTH must be a positive finite scalar.');
    end
    if ~isnumeric(M) || prod(size(M)) ~= 1 || ~is_finite_value(M) || M <= 0
        error('M must be a positive finite scalar.');
    end

    results = initialize_results_struct();
    results.Hm0 = HM0;
    results.d = DEPTH;
    results.slopeM = M;

    [results, ok] = perform_wave_analysis(results);
    if ~ok
        error('Calculation failed.');
    end

    report = format_report(results);
    write_report_to_file(report);

    if nargout == 0
        fprintf('%s', report);
    end
end


function results = initialize_results_struct()
    results = struct( ...
        'Hm0', 0.0, ...
        'd', 0.0, ...
        'slopeM', 0.0, ...
        'distribution_type', '', ...
        'Hrms', 0.0, ...
        'm0', 0.0, ...
        'tanAlpha', 0.0, ...
        'Htr_dim', 0.0, ...
        'Htr_tilde', 0.0, ...
        'H1_Hrms', 0.0, ...
        'H2_Hrms', 0.0, ...
        'H1_3_Hrms', 0.0, ...
        'H1_10_Hrms', 0.0, ...
        'H1_50_Hrms', 0.0, ...
        'H1_100_Hrms', 0.0, ...
        'H1_250_Hrms', 0.0, ...
        'H1_1000_Hrms', 0.0, ...
        'H1_dim', 0.0, ...
        'H2_dim', 0.0, ...
        'H1_3_dim', 0.0, ...
        'H1_10_dim', 0.0, ...
        'H1_50_dim', 0.0, ...
        'H1_100_dim', 0.0, ...
        'H1_250_dim', 0.0, ...
        'H1_1000_dim', 0.0, ...
        'ratio_1_10_div_1_3', 0.0, ...
        'ratio_1_50_div_1_3', 0.0, ...
        'ratio_1_100_div_1_3', 0.0, ...
        'ratio_1_250_div_1_3', 0.0, ...
        'ratio_1_1000_div_1_3', 0.0);
end


function [results, ok] = perform_wave_analysis(results)
    try
        results.m0 = (results.Hm0 / 4.0)^2;
        sqrt_m0 = sqrt(results.m0);
        results.Hrms = (2.69 + 3.24 * sqrt_m0 / results.d) * sqrt_m0;

        results.tanAlpha = 1.0 / results.slopeM;
        results.Htr_dim = (0.35 + 5.8 * results.tanAlpha) * results.d;
        if results.Hrms > 0.0
            results.Htr_tilde = results.Htr_dim / results.Hrms;
        else
            results.Htr_tilde = 0.0;
        end

        if results.Htr_tilde > 2.75
            results.distribution_type = 'Rayleigh';

            results.H1_3_dim    = 1.001075736951740 * results.Hm0;
            results.H1_10_dim   = 1.272734273369137 * results.Hm0;
            results.H1_50_dim   = 1.560113379974762 * results.Hm0;
            results.H1_100_dim  = 1.668233372358517 * results.Hm0;
            results.H1_250_dim  = 1.801017222497626 * results.Hm0;
            results.H1_1000_dim = 1.984835590575388 * results.Hm0;

            if results.Hrms > 0.0
                results.H1_3_Hrms    = results.H1_3_dim    / results.Hrms;
                results.H1_10_Hrms   = results.H1_10_dim   / results.Hrms;
                results.H1_50_Hrms   = results.H1_50_dim   / results.Hrms;
                results.H1_100_Hrms  = results.H1_100_dim  / results.Hrms;
                results.H1_250_Hrms  = results.H1_250_dim  / results.Hrms;
                results.H1_1000_Hrms = results.H1_1000_dim / results.Hrms;
            end
        else
            results.distribution_type = 'B&G';

            [results.H1_Hrms, results.H2_Hrms, ok_solver] = newtonRaphsonSystemSolver(results.Htr_tilde, 100);
            if ~ok_solver
                error('Newton-Raphson solver failed to converge.');
            end

            results.H1_3_Hrms    = calculate_H1N(3.0,    results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_10_Hrms   = calculate_H1N(10.0,   results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_50_Hrms   = calculate_H1N(50.0,   results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_100_Hrms  = calculate_H1N(100.0,  results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_250_Hrms  = calculate_H1N(250.0,  results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);
            results.H1_1000_Hrms = calculate_H1N(1000.0, results.H1_Hrms, results.H2_Hrms, results.Htr_tilde);

            results.H1_dim = results.H1_Hrms * results.Hrms;
            results.H2_dim = results.H2_Hrms * results.Hrms;
            results.H1_3_dim = results.H1_3_Hrms * results.Hrms;
            results.H1_10_dim = results.H1_10_Hrms * results.Hrms;
            results.H1_50_dim = results.H1_50_Hrms * results.Hrms;
            results.H1_100_dim = results.H1_100_Hrms * results.Hrms;
            results.H1_250_dim = results.H1_250_Hrms * results.Hrms;
            results.H1_1000_dim = results.H1_1000_Hrms * results.Hrms;

            results.H1_3_dim    = min(results.H1_3_dim,    1.001075736951740 * results.Hm0);
            results.H1_10_dim   = min(results.H1_10_dim,   1.272734273369137 * results.Hm0);
            results.H1_50_dim   = min(results.H1_50_dim,   1.560113379974762 * results.Hm0);
            results.H1_100_dim  = min(results.H1_100_dim,  1.668233372358517 * results.Hm0);
            results.H1_250_dim  = min(results.H1_250_dim,  1.801017222497626 * results.Hm0);
            results.H1_1000_dim = min(results.H1_1000_dim, 1.984835590575388 * results.Hm0);
        end

        if results.H1_3_dim > 0.0
            H1_3_Hrms_capped    = results.H1_3_dim    / results.Hrms;
            H1_10_Hrms_capped   = results.H1_10_dim   / results.Hrms;
            H1_50_Hrms_capped   = results.H1_50_dim   / results.Hrms;
            H1_100_Hrms_capped  = results.H1_100_dim  / results.Hrms;
            H1_250_Hrms_capped  = results.H1_250_dim  / results.Hrms;
            H1_1000_Hrms_capped = results.H1_1000_dim / results.Hrms;

            results.ratio_1_10_div_1_3   = H1_10_Hrms_capped   / H1_3_Hrms_capped;
            results.ratio_1_50_div_1_3   = H1_50_Hrms_capped   / H1_3_Hrms_capped;
            results.ratio_1_100_div_1_3  = H1_100_Hrms_capped  / H1_3_Hrms_capped;
            results.ratio_1_250_div_1_3  = H1_250_Hrms_capped  / H1_3_Hrms_capped;
            results.ratio_1_1000_div_1_3 = H1_1000_Hrms_capped / H1_3_Hrms_capped;
        else
            results.ratio_1_10_div_1_3   = 0.0;
            results.ratio_1_50_div_1_3   = 0.0;
            results.ratio_1_100_div_1_3  = 0.0;
            results.ratio_1_250_div_1_3  = 0.0;
            results.ratio_1_1000_div_1_3 = 0.0;
        end

        ok = true;
    catch
        ok = false;
    end
end


function [H1_Hrms, H2_Hrms, ok] = newtonRaphsonSystemSolver(Htr_Hrms, maxit)
    EPSILON = 1e-12;
    JACOBIAN_DX = 1e-8;
    LOCAL_EPS = 1e-16;

    [H1_Hrms, H2_Hrms] = get_initial_guesses(Htr_Hrms);
    ok = false;

    for iter = 0:(maxit - 1)
        f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        if abs(f1_val) < EPSILON && abs(f2_val) < EPSILON
            ok = true;
            return;
        end

        J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);

        [dH1, dH2] = solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val);

        H1_Hrms = H1_Hrms + dH1;
        H2_Hrms = H2_Hrms + dH2;

        if H1_Hrms <= 0.0
            H1_Hrms = LOCAL_EPS;
        end
        if H2_Hrms <= 0.0
            H2_Hrms = LOCAL_EPS;
        end
    end
end


function val = F1(H1_Hrms, H2_Hrms, Htr_Hrms)
    K1 = 2.0;
    K2 = 3.6;
    DBL_MAX = 1.7976931348623157e308;

    if H1_Hrms <= 0.0 || H2_Hrms <= 0.0
        val = DBL_MAX;
        return;
    end

    arg1 = (Htr_Hrms / H1_Hrms)^K1;
    arg2 = (Htr_Hrms / H2_Hrms)^K2;
    term1 = H1_Hrms * H1_Hrms * incomplete_gamma_lower(2.0 / K1 + 1.0, arg1);
    term2 = H2_Hrms * H2_Hrms * incomplete_gamma_upper(2.0 / K2 + 1.0, arg2);
    val = sqrt(max(0.0, term1 + term2)) - 1.0;
end


function val = F2(H1_Hrms, H2_Hrms, Htr_Hrms)
    K1 = 2.0;
    K2 = 3.6;
    DBL_MAX = 1.7976931348623157e308;

    if H1_Hrms <= 0.0 || H2_Hrms <= 0.0
        val = DBL_MAX;
        return;
    end

    val = (Htr_Hrms / H1_Hrms)^K1 - (Htr_Hrms / H2_Hrms)^K2;
end


function [dx1, dx2] = solve_linear_system_2x2(J11, J12, J21, J22, b1, b2)
    determinant = J11 * J22 - J12 * J21;
    if abs(determinant) < 1e-20
        dx1 = 0.0;
        dx2 = 0.0;
        return;
    end

    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
end


function [H1_initial, H2_initial] = get_initial_guesses(Htr_Hrms)
    LOCAL_EPS = 1e-16;

    H1_initial = 2.244660800090239e-03 + ...
                 tanh(1.918610494219390e+00 * Htr_Hrms)^1.780892753373355e-01 / ...
                 tanh(sinh(1.009497360864962e+00 * Htr_Hrms))^9.777939607559606e-01;
    H2_initial = 1.059259665431797 + (0.2059286860468916 * Htr_Hrms) / ...
                 (1.0 + 3.865701948059343 * Htr_Hrms^(-3.479682433107255));

    if H1_initial <= 0.0 || isnan(H1_initial)
        H1_initial = LOCAL_EPS;
    end
    if H2_initial <= 0.0 || isnan(H2_initial)
        H2_initial = LOCAL_EPS;
    end
end


function HN = calculate_HN(N, H1, H2, Htr)
    EPSILON = 1e-12;
    K1 = 2.0;
    K2 = 3.6;

    if N <= 1.0 || H1 <= 0.0 || H2 <= 0.0
        error('Invalid arguments for calculate_HN.');
    end

    HN_candidate1 = H1 * (log(N)^(1.0 / K1));
    if HN_candidate1 < Htr - EPSILON
        HN = HN_candidate1;
    else
        HN = H2 * (log(N)^(1.0 / K2));
    end
end


function H1N = calculate_H1N(N_val, H1, H2, Htr)
    EPSILON = 1e-12;
    K1 = 2.0;
    K2 = 3.6;

    if N_val <= 1.0 || H1 <= 0.0 || H2 <= 0.0
        error('Invalid arguments for calculate_H1N.');
    end

    H_N_val = calculate_HN(N_val, H1, H2, Htr);
    term1_x_ln_Nval = log(N_val);
    term2_a = 1.0 / K2 + 1.0;

    if H_N_val < Htr - EPSILON
        term1_a = 1.0 / K1 + 1.0;
        term1_x_HtrH1 = (Htr / H1)^K1;
        term2_x_HtrH2 = (Htr / H2)^K2;

        gamma1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval) - incomplete_gamma_upper(term1_a, term1_x_HtrH1);
        gamma2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2);
        H1N = N_val * (H1 * gamma1 + H2 * gamma2);
    else
        H1N = N_val * H2 * incomplete_gamma_upper(term2_a, term1_x_ln_Nval);
    end
end


function val = incomplete_gamma_lower(a, x)
    MAXIT = 500;
    LOCAL_EPS = 1e-16;
    DBL_MIN = 2.2250738585072014e-308;

    if a <= 0.0 || x < 0.0
        val = NaN;
        return;
    end
    if x == 0.0
        val = 0.0;
        return;
    end

    gln = gammaln(a);

    if x < a + 1.0
        ap = a;
        sum_val = 1.0 / a;
        del = sum_val;
        for n_iter = 1:MAXIT
            ap = ap + 1.0;
            del = del * x / ap;
            sum_val = sum_val + del;
            if abs(del) < abs(sum_val) * LOCAL_EPS
                val = sum_val * exp(-x + a * log(x) - gln) * gamma(a);
                return;
            end
        end
        val = NaN;
    else
        b = x + 1.0 - a;
        c = 1.0 / DBL_MIN;
        d = 1.0 / b;
        h = d;
        for i = 1:MAXIT
            an = -1.0 * i * (i - a);
            b = b + 2.0;
            d = an * d + b;
            if abs(d) < DBL_MIN
                d = DBL_MIN;
            end
            c = b + an / c;
            if abs(c) < DBL_MIN
                c = DBL_MIN;
            end
            d = 1.0 / d;
            del = d * c;
            h = h * del;
            if abs(del - 1.0) < LOCAL_EPS
                break;
            end
        end
        p_normalized = 1.0 - (exp(-x + a * log(x) - gln) * h);
        val = max(0.0, p_normalized) * gamma(a);
    end
end


function val = incomplete_gamma_upper(a, x)
    if a <= 0.0 || x < 0.0
        error('Invalid arguments for incomplete_gamma_upper.');
    end
    val = gamma(a) - incomplete_gamma_lower(a, x);
end


function report = format_report(r)
    nl = sprintf('\n');
    report = [ ...
        '======================', nl, ...
        '   INPUT PARAMETERS', nl, ...
        '======================', nl, ...
        sprintf('Hm0 (m)         : %.4f%s', r.Hm0, nl), ...
        sprintf('d (m)           : %.4f%s', r.d, nl), ...
        sprintf('Beach slope (m) : %.4f   (tan(alpha) = %.4f)%s', r.slopeM, r.tanAlpha, nl), ...
        nl, ...
        'Distribution Used : ', r.distribution_type, nl, ...
        nl, ...
        '===========================', nl, ...
        '   CALCULATED PARAMETERS', nl, ...
        '===========================', nl, ...
        sprintf('Free-surface variance m0 (m^2)   : %.4f%s', r.m0, nl), ...
        sprintf('Mean square wave height Hrms (m) : %.4f%s', r.Hrms, nl), ...
        sprintf('Transitional wave height Htr (m) : %.4f%s', r.Htr_dim, nl), ...
        sprintf('Dimensionless H~_tr (Htr/Hrms)   : %.4f%s', r.Htr_tilde, nl), ...
        nl, ...
        '=========================================', nl, ...
        '   DIMENSIONLESS WAVE HEIGHTS (H/Hrms)', nl, ...
        '=========================================', nl, ...
        sprintf('H1/Hrms       : %.4f%s', r.H1_Hrms, nl), ...
        sprintf('H2/Hrms       : %.4f%s', r.H2_Hrms, nl), ...
        sprintf('H1/3 / Hrms   : %.4f%s', r.H1_3_Hrms, nl), ...
        sprintf('H1/10 / Hrms  : %.4f%s', r.H1_10_Hrms, nl), ...
        sprintf('H1/50 / Hrms  : %.4f%s', r.H1_50_Hrms, nl), ...
        sprintf('H1/100 / Hrms : %.4f%s', r.H1_100_Hrms, nl), ...
        sprintf('H1/250 / Hrms : %.4f%s', r.H1_250_Hrms, nl), ...
        sprintf('H1/1000 /Hrms : %.4f%s', r.H1_1000_Hrms, nl), ...
        nl, ...
        '==================================', nl, ...
        '   DIMENSIONAL WAVE HEIGHTS (m)', nl, ...
        '==================================', nl, ...
        sprintf('H1 (m)        : %.4f%s', r.H1_dim, nl), ...
        sprintf('H2 (m)        : %.4f%s', r.H2_dim, nl), ...
        sprintf('H1/3 (m)      : %.4f%s', r.H1_3_dim, nl), ...
        sprintf('H1/10 (m)     : %.4f%s', r.H1_10_dim, nl), ...
        sprintf('H1/50 (m)     : %.4f%s', r.H1_50_dim, nl), ...
        sprintf('H1/100 (m)    : %.4f%s', r.H1_100_dim, nl), ...
        sprintf('H1/250 (m)    : %.4f%s', r.H1_250_dim, nl), ...
        sprintf('H1/1000 (m)   : %.4f%s', r.H1_1000_dim, nl), ...
        nl, ...
        '=======================', nl, ...
        '   DIAGNOSTIC RATIOS', nl, ...
        '=======================', nl, ...
        sprintf('(H1/10)/(H1/3)   : %.4f%s', r.ratio_1_10_div_1_3, nl), ...
        sprintf('(H1/50)/(H1/3)   : %.4f%s', r.ratio_1_50_div_1_3, nl), ...
        sprintf('(H1/100)/(H1/3)  : %.4f%s', r.ratio_1_100_div_1_3, nl), ...
        sprintf('(H1/250)/(H1/3)  : %.4f%s', r.ratio_1_250_div_1_3, nl), ...
        sprintf('(H1/1000)/(H1/3) : %.4f%s', r.ratio_1_1000_div_1_3, nl), ...
        nl, ...
        'End of Report', nl ...
    ];
end


function write_report_to_file(report)
    fid = fopen('report.txt', 'w');
    if fid == -1
        error('Could not open report.txt for writing.');
    end
    fprintf(fid, '%s', report);
    fclose(fid);
end


function tf = is_finite_value(x)
    tf = (~isnan(x)) && (~isinf(x));
end
