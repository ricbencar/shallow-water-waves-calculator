function results = shallow_water_waves(HM0, DEPTH, M)
%SHALLOW_WATER_WAVES Local individual-wave statistics on shallow foreshores.
%   R = SHALLOW_WATER_WAVES(HM0, DEPTH, M) calculates the Composite Weibull
%   Distribution (Battjes and Groenendijk) and Rayleigh-limited statistics.
%   HM0 and DEPTH are in metres; M is the denominator of a 1:M slope.
%   R contains numerical results and individual-wave statistics.
%   The report is written to report.txt. Without an output argument it is
%   also displayed in the MATLAB command window.
%
%   Example: R = shallow_water_waves(2.5, 5.0, 100);

    validate_positive_scalar(HM0, 'HM0');
    validate_positive_scalar(DEPTH, 'DEPTH');
    validate_positive_scalar(M, 'M');
    HM0 = double(HM0); DEPTH = double(DEPTH); M = double(M);
    k1 = 2.0; k2 = 3.6;
    ns = [3 10 50 100 250 1000];
    rayleighFactors = [1.001075736951740, 1.272734273369137, ...
        1.560113379974762, 1.668233372358517, ...
        1.801017222497626, 1.984835590575388];
    results = struct();
    results.Hm0 = HM0; results.d = DEPTH; results.slopeM = M;
    results.depth = DEPTH; results.slope_denominator = M;
    results.tanAlpha = 1/M; results.slope = results.tanAlpha;
    sigma = positive(HM0/4, 'Spectral standard deviation');
    results.m0 = positive(sigma*sigma, 'Spectral variance');
    results.saturation = positive(sigma/DEPTH, 'Saturation');
    results.Hrms = positive((2.69+3.24*results.saturation)*sigma, 'Hrms');
    results.Htr_dim = positive((0.35+5.8/M)*DEPTH, 'Transition height');
    results.Htr = results.Htr_dim;
    results.Htr_tilde = positive(results.Htr/results.Hrms, 'Normalized transition');
    results.transition_normalized = results.Htr_tilde;
    results.uses_cwd = results.Htr_tilde <= 2.75;
    results.distribution_type = 'Rayleigh';
    if results.uses_cwd, results.distribution_type = 'Composite Weibull (B&G)'; end
    results.slope_outside_calibration = M < 20 || M > 250;
    results.saturation_outside_observed_range = results.saturation < 0.025 || results.saturation > 0.16;
    results.H1_Hrms = NaN; results.H2_Hrms = NaN;
    results.H1_dim = NaN; results.H2_dim = NaN;
    results.transition_exponent = NaN;
    results.transition_exceedance_probability = NaN;
    results.transition_recurrence_count = NaN;
    results.moment_residual = NaN;
    results.continuity_residual = NaN;
    results.solver_iterations = 0;

    if results.uses_cwd
        t = results.Htr_tilde;
        [q, z, iterations] = solve_transition_exponent(t, k1, k2);
        a1 = positive(t*exp(-z/k1), 'First Weibull scale');
        a2 = positive(t*exp(-z/k2), 'Second Weibull scale');
        results.H1_Hrms = a1; results.H2_Hrms = a2;
        results.H1_dim = a1*results.Hrms;
        results.H2_dim = a2*results.Hrms;
        results.transition_exponent = q;
        results.transition_exceedance_probability = exp(-q);
        results.transition_recurrence_count = exp(q);
        results.solver_iterations = iterations;
        m2 = a1^2*gamma_lower(1+2/k1,q) + a2^2*gamma_upper(1+2/k2,q);
        results.moment_residual = sqrt(m2)-1;
        results.continuity_residual = (t/a1)^k1-(t/a2)^k2;
        if ~finite_value(m2) || abs(results.moment_residual)>1e-10 || ...
                abs(results.continuity_residual)>1e-10
            error('CWD probability constraints were not satisfied.');
        end
    end

    rayleighScale = positive(HM0/sqrt(2),'Rayleigh scale');
    statistics = repmat(struct('N',0,'threshold',0,'rayleigh_threshold',0, ...
        'mean_unlimited',0,'rayleigh_mean',0,'mean_final',0, ...
        'mean_normalized',0,'rayleigh_limit_applied',false),1,numel(ns));
    for i = 1:numel(ns)
        N = ns(i); lnN = log(N);
        s = statistics(i); s.N = N;
        s.rayleigh_threshold = positive(rayleighScale*sqrt(lnN),'Rayleigh threshold');
        s.rayleigh_mean = positive(rayleighFactors(i)*HM0,'Rayleigh mean');
        if results.uses_cwd
            if lnN < q
                hn = a1*lnN^(1/k1);
                moment = a1*gamma_interval(1+1/k1,lnN,q) + ...
                         a2*gamma_upper(1+1/k2,q);
            else
                hn = a2*lnN^(1/k2);
                moment = a2*gamma_upper(1+1/k2,lnN);
            end
            s.threshold = positive(hn*results.Hrms,'CWD threshold');
            s.mean_unlimited = positive(N*moment*results.Hrms,'CWD mean');
            s.mean_final = min(s.mean_unlimited,s.rayleigh_mean);
            s.rayleigh_limit_applied = s.mean_unlimited > s.rayleigh_mean*(1+8*eps);
        else
            s.threshold = s.rayleigh_threshold;
            s.mean_unlimited = s.rayleigh_mean;
            s.mean_final = s.rayleigh_mean;
        end
        s.mean_normalized = positive(s.mean_final/results.Hrms,'Normalized mean');
        if ~(s.mean_unlimited > s.threshold) || s.mean_final > s.rayleigh_mean*(1+1e-12)
            error('Inconsistent wave-height statistics.');
        end
        if i>1
            previous = statistics(i-1);
            if ~(s.threshold>previous.threshold && ...
                 s.mean_unlimited>previous.mean_unlimited && ...
                 s.mean_final>previous.mean_final)
                error('Wave-height ordering constraint not satisfied.');
            end
        end
        statistics(i) = s;
        suffix = sprintf('%d',N);
        results.(['H1_' suffix '_dim']) = s.mean_final;
        results.(['H1_' suffix '_Hrms']) = s.mean_normalized;
        results.(['H_' suffix '_threshold']) = s.threshold;
        results.(['H1_' suffix '_uncapped']) = s.mean_unlimited;
        results.(['H1_' suffix '_rayleigh']) = s.rayleigh_mean;
        results.(['H1_' suffix '_rayleigh_limit_applied']) = s.rayleigh_limit_applied;
        if i>1
            results.(['ratio_1_' suffix '_div_1_3']) = s.mean_final/statistics(1).mean_final;
        end
    end
    results.statistics = statistics;
    report = format_report(results);
    % FreeMat 4.2: fwrite returns the byte count for a binary write.
    % The report contains ASCII text; write its bytes without text-mode
    % translation or assumptions about fprintf's return value.
    bytes = uint8(report);
    fid = fopen('report.txt','w');
    if fid < 0, error('Could not open report.txt for writing.'); end
    nwrite = fwrite(fid, bytes, 'uint8');
    fclose(fid);
    if nwrite ~= length(bytes)
        error('Could not write report.txt completely.');
    end
    if nargout==0, fprintf('%s\n',report); end
end

% Finite-value predicate compatible with FreeMat 4.2.
% For finite real values x-x is zero; NaN and infinities produce NaN.
function ok = finite_value(x)
    ok = (x == x) && ((x - x) == 0);
end

function v = positive(v, name)
    if ~(finite_value(v) && v>0), error('%s must be finite and positive.',name); end
end

function validate_positive_scalar(x, name)
    if ~isnumeric(x) || ~isreal(x) || length(x) ~= 1 || ...
            ~finite_value(x) || x <= 0
        error('%s must be a finite positive real scalar.', name);
    end
end

function y = gamma_lower(a,x)
    [p,q] = gamma_regularized_pair(a,x);
    y = gamma(a)*p;
end

function y = gamma_upper(a,x)
    [p,q] = gamma_regularized_pair(a,x);
    y = gamma(a)*q;
end

function [p,q] = gamma_regularized_pair(a,x)
    if ~finite_value(a) || a<=0 || ~finite_value(x) || x<0
        error('Invalid incomplete gamma arguments.');
    end
    if x==0
        p=0; q=1; return;
    end
    tol=2.5e-14;
    maxIter=10000;
    factor=exp(a*log(x)-x-log(gamma(a)));
    if x<a+1
        ap=a; delta=1/a; total=delta; converged=0;
        for n=1:maxIter
            ap=ap+1;
            delta=delta*x/ap;
            total=total+delta;
            if abs(delta)<=tol*abs(total)
                converged=1; break;
            end
        end
        if ~converged, error('Incomplete gamma series did not converge.'); end
        p=total*factor; q=1-p;
    else
        fpmin=1e-300;
        b=x+1-a; c=1/fpmin; d=1/b; fraction=d; converged=0;
        for n=1:maxIter
            an=-n*(n-a);
            b=b+2;
            d=an*d+b;
            if abs(d)<fpmin, d=sign_nonzero(d)*fpmin; end
            c=b+an/c;
            if abs(c)<fpmin, c=sign_nonzero(c)*fpmin; end
            d=1/d;
            delta=c*d;
            fraction=fraction*delta;
            if abs(delta-1)<=tol
                converged=1; break;
            end
        end
        if ~converged, error('Incomplete gamma continued fraction did not converge.'); end
        q=fraction*factor; p=1-q;
    end
    if ~finite_value(p) || ~finite_value(q) || p < -2e-13 || q < -2e-13 || ...
            p > 1+2e-13 || q > 1+2e-13
        error('Incomplete gamma probability outside valid range.');
    end
    p=max(0,min(1,p)); q=max(0,min(1,q));
end

function y=sign_nonzero(x)
    if x<0, y=-1; else, y=1; end
end

function y = gamma_interval(a,lo,hi)
    if hi==lo, y=0; return; end
    if ~finite_value(lo) || ~finite_value(hi) || lo<0 || hi<lo
        error('Invalid incomplete gamma integration interval.');
    end
    % Sixteen-point Gauss-Legendre quadrature.
    nodes=[0.09501250983763744 0.28160355077925891 ...
           0.45801677765722739 0.61787624440264375 ...
           0.75540440835500303 0.86563120238783174 ...
           0.94457502307323258 0.98940093499164993];
    weights=[0.18945061045506850 0.18260341504492359 ...
             0.16915651939500254 0.14959598881657673 ...
             0.12462897125553387 0.09515851168249278 ...
             0.06225352393864789 0.02715245941175409];
    midpoint=(lo+hi)/2; half=(hi-lo)/2; total=0;
    for i=1:8
        u1=midpoint-half*nodes(i);
        u2=midpoint+half*nodes(i);
        total=total+weights(i)*(u1^(a-1)*exp(-u1)+u2^(a-1)*exp(-u2));
    end
    y=positive(total*half,'Gamma interval integral');
end

function [q, z, iterations] = solve_transition_exponent(t, k1, k2)
    positive(t, 'Normalized transition height');
    q0 = 1.002054709*t*t*tanh(t^1.652259933) * ...
         exp(-tanh(tanh(t))*exp(-2.146235748*t));
    positive(q0, 'Newton initial transition exponent');
    z = log(q0);
    logt = log(t);
    maxIterations = 40;
    rootTolerance = 5e-13;
    for iteration = 0:maxIterations
        [residual, derivative] = evaluate_log_moment(z, logt, k1, k2);
        if abs(residual) <= rootTolerance
            iterations = iteration;
            q = positive(exp(z), 'CWD transition exponent');
            return;
        end
        if iteration == maxIterations
            break;
        end
        z = z - residual/derivative;
        if ~finite_value(z)
            error('Non-finite Newton iteration.');
        end
    end
    error('CWD Newton solver did not converge.');
end

function [residual, derivative] = evaluate_log_moment(z, logt, k1, k2)
    q = positive(exp(z), 'Transition exponent q');
    lower = exp(-2*z/k1)*gamma_lower(1+2/k1,q);
    upper = exp(-2*z/k2)*gamma_upper(1+2/k2,q);
    kernel = positive(lower+upper, 'Second-moment kernel');
    residual = 2*logt + log(kernel);
    % Exact logarithmic derivative; the gamma boundary terms cancel.
    derivative = -((2/k1)*lower+(2/k2)*upper)/kernel;
    if ~finite_value(residual) || ~finite_value(derivative) || ~(derivative < 0)
        error('Invalid Newton residual or derivative.');
    end
end

function text = format_report(r)
    lines = {};
    lines{end+1}='SHALLOW-WATER INDIVIDUAL WAVE-HEIGHT STATISTICS';
    lines{end+1}=repmat('=',1,72);
    lines{end+1}=''; lines{end+1}='INPUT PARAMETERS';
    lines{end+1}=field('Spectral significant wave height Hm0 (m)',sprintf('%.4f',r.Hm0));
    lines{end+1}=field('Still-water depth d (m)',sprintf('%.4f',r.d));
    slopeText = sprintf('%.4f',r.slopeM);
    while length(slopeText)>0 && slopeText(length(slopeText))=='0'
        slopeText=slopeText(1:length(slopeText)-1);
    end
    if length(slopeText)>0 && slopeText(length(slopeText))=='.'
        slopeText=slopeText(1:length(slopeText)-1);
    end
    lines{end+1}=field('Foreshore slope',['1:' slopeText]);
    lines{end+1}=field('tan(alpha)',sprintf('%.6f',r.tanAlpha));
    lines{end+1}=''; lines{end+1}='LOCAL WAVE PARAMETERS';
    lines{end+1}=field('Spectral variance m0 (m^2)',sprintf('%.6f',r.m0));
    lines{end+1}=field('Individual-wave Hrms (m)',sprintf('%.6f',r.Hrms));
    lines{end+1}=field('Degree of saturation sqrt(m0)/d',sprintf('%.6f',r.saturation));
    lines{end+1}=field('Transitional wave height Htr (m)',sprintf('%.6f',r.Htr));
    lines{end+1}=field('Normalized transition Htr/Hrms',sprintf('%.6f',r.Htr_tilde));
    lines{end+1}=field('Distribution',r.distribution_type,false);
    lines{end+1}=''; lines{end+1}='DISTRIBUTION PARAMETERS';
    if r.uses_cwd
        lines{end+1}=field('Normalized first Weibull scale H1/Hrms',sprintf('%.6f',r.H1_Hrms));
        lines{end+1}=field('Normalized second Weibull scale H2/Hrms',sprintf('%.6f',r.H2_Hrms));
        lines{end+1}=field('Dimensional Weibull scale H1 (m)',sprintf('%.6f',r.H1_dim));
        lines{end+1}=field('Dimensional Weibull scale H2 (m)',sprintf('%.6f',r.H2_dim));
        lines{end+1}=field('Transition exponent q',sprintf('%.6f',r.transition_exponent));
        lines{end+1}=field('Probability P(H > Htr)',sprintf('%.9f',r.transition_exceedance_probability));
        lines{end+1}=field('Equivalent transition recurrence Ntr',sprintf('%.6f waves',r.transition_recurrence_count));
    else
        lines{end+1}=field('Composite Weibull scales H1 and H2','N/A',false);
        lines{end+1}=field('Composite Weibull transition exponent','N/A',false);
    end
    lines{end+1}=''; lines{end+1}='WAVE-HEIGHT STATISTICS (m)';
    lines{end+1}=statrow('HIGHEST','H_N','RAW MEAN','RAYLEIGH','FINAL MEAN','LIMIT');
    lines{end+1}=statrow('FRACTION','THRESHOLD','H_(1/N)','H_(1/N)','H_(1/N)','APPLIED');
    lines{end+1}=['  ' repmat('-',1,70)];
    for si=1:length(r.statistics)
        s=r.statistics(si);
        if s.rayleigh_limit_applied, limit='YES'; else, limit='NO'; end
        lines{end+1}=statrow(sprintf('1/%d',s.N),sprintf('%.4f',s.threshold), ...
            sprintf('%.4f',s.mean_unlimited),sprintf('%.4f',s.rayleigh_mean), ...
            sprintf('%.4f',s.mean_final),limit);
    end
    lines{end+1}=''; lines{end+1}='FINAL NORMALIZED MEAN-HIGH-WAVE HEIGHTS';
    for si=1:length(r.statistics)
        s=r.statistics(si);
        lines{end+1}=field(sprintf('H_(1/%d) / Hrms',s.N),sprintf('%.7f',s.mean_normalized));
    end
    lines{end+1}=''; lines{end+1}='RATIOS RELATIVE TO THE HIGHEST ONE THIRD';
    for si=2:length(r.statistics)
        s=r.statistics(si);
        lines{end+1}=field(sprintf('H_(1/%d) / H_(1/3)',s.N), ...
            sprintf('%.6f',s.mean_final/r.statistics(1).mean_final));
    end
    lines{end+1}=''; lines{end+1}='NUMERICAL DIAGNOSTICS';
    if r.uses_cwd
        lines{end+1}=field('Scalar root iterations',sprintf('%d',r.solver_iterations));
        lines{end+1}=field('Second-moment residual',sprintf('%.10e',r.moment_residual));
        lines{end+1}=field('Continuity exponent residual',sprintf('%.10e',r.continuity_residual));
    else
        lines{end+1}=field('CWD nonlinear solver','N/A',false);
    end
    lines{end+1}=''; lines{end+1}='INTERPRETATION';
    lines{end+1}='  H_N is the height exceeded with probability 1/N.';
    lines{end+1}='  H_(1/N) is the mean height of the highest 1/N fraction.';
    lines{end+1}='  Both are individual-wave statistics, not storm maxima.';
    lines{end+1}='  Mean-high-wave Rayleigh limits are applied separately.';
    lines{end+1}='  H_N always refers to the selected underlying distribution.';
    lines{end+1}='  The program uses local wave conditions; it does not propagate waves.';
    if r.slope_outside_calibration || (r.uses_cwd && r.saturation_outside_observed_range)
        lines{end+1}=''; lines{end+1}='APPLICABILITY NOTICE';
        if r.slope_outside_calibration
            lines{end+1}='  Foreshore slope is outside the laboratory slope range';
            lines{end+1}='  1:20 to 1:250 of the Composite Weibull formulation.';
            lines{end+1}='  Predictions for an extensive horizontal, constant-depth';
            lines{end+1}='  seabed require a separate applicability assessment.';
        end
        if r.uses_cwd && r.saturation_outside_observed_range
            lines{end+1}='  The degree of saturation is outside the approximate';
            lines{end+1}='  interval represented by the original foreshore tests.';
            lines{end+1}='  Interpretation requires engineering assessment.';
        end
    end
    text = '';
    for li=1:length(lines)
        if li>1, text=[text sprintf('\n')]; end
        text=[text lines{li}];
    end
end

function str = field(label,value,right)
    if nargin<3, right=true; end
    if right, str=sprintf('  %-42s : %16s',label,value);
    else, str=sprintf('  %-42s : %s',label,value); end
end

function str = statrow(fr,hn,raw,ray,fin,lim)
    str=sprintf('  %-10s %11s %12s %12s %12s %8s',fr,hn,raw,ray,fin,lim);
end
