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

    validateattributes(HM0, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'HM0', 1);
    validateattributes(DEPTH, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'DEPTH', 2);
    validateattributes(M, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'M', 3);
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
        logt = log(t);
        f = @(z) scalar_log_moment(z, logt, k1, k2);
        za = -700; zb = log(64);
        if ~(isfinite(f(za)) && isfinite(f(zb)) && f(za)>0 && f(zb)<0)
            error('CWD scalar solution is outside the representable range.');
        end
        % MATLAB fzero employs a safeguarded Brent-type scalar algorithm.
        options = optimset('TolX',2e-14,'MaxIter',200,'Display','off');
        [z, residual, exitflag, output] = fzero(f,[za zb],options);
        if exitflag <= 0 || ~isfinite(z) || abs(residual)>1e-10
            error('CWD scalar solver did not converge.');
        end
        q = positive(exp(z), 'Transition exponent');
        a1 = positive(t*exp(-z/k1), 'First Weibull scale');
        a2 = positive(t*exp(-z/k2), 'Second Weibull scale');
        results.H1_Hrms = a1; results.H2_Hrms = a2;
        results.H1_dim = a1*results.Hrms;
        results.H2_dim = a2*results.Hrms;
        results.transition_exponent = q;
        results.transition_exceedance_probability = exp(-q);
        results.transition_recurrence_count = exp(q);
        results.solver_iterations = output.iterations;
        m2 = a1^2*gamma_lower(1+2/k1,q) + a2^2*gamma_upper(1+2/k2,q);
        results.moment_residual = sqrt(m2)-1;
        results.continuity_residual = (t/a1)^k1-(t/a2)^k2;
        if ~isfinite(m2) || abs(results.moment_residual)>1e-10 || ...
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
    fid = fopen('report.txt','w','n','UTF-8');
    if fid < 0, error('Could not open report.txt for writing.'); end
    try
        nwrite = fprintf(fid,'%s',report);
        if nwrite < numel(report), error('Incomplete report write.'); end
        if fclose(fid)~=0, error('Could not close report.txt.'); end
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    if nargout==0, fprintf('%s\n',report); end
end

function v = positive(v, name)
    if ~(isfinite(v) && v>0), error('%s must be finite and positive.',name); end
end

function y = gamma_lower(a,x)
    y = gamma(a)*gammainc(x,a,'lower');
end

function y = gamma_upper(a,x)
    y = gamma(a)*gammainc(x,a,'upper');
end

function y = gamma_interval(a,lo,hi)
    if hi==lo, y=0; return; end
    if ~(isfinite(lo) && isfinite(hi) && lo>=0 && hi>=lo)
        error('Invalid incomplete gamma integration interval.');
    end
    % Direct integration avoids subtracting nearly equal tail functions.
    y = integral(@(u)u.^(a-1).*exp(-u),lo,hi, ...
        'RelTol',2e-13,'AbsTol',1e-14);
    y = positive(y,'Gamma interval integral');
end

function y = scalar_log_moment(z, logt, k1, k2)
    q = exp(z);
    kernel = exp(-2*z/k1)*gamma_lower(1+2/k1,q) + ...
             exp(-2*z/k2)*gamma_upper(1+2/k2,q);
    if ~(isfinite(kernel)&&kernel>0), error('Invalid CWD second-moment kernel.'); end
    y = 2*logt + log(kernel);
end

function text = format_report(r)
    lines = {};
    lines{end+1}='SHALLOW-WATER INDIVIDUAL WAVE-HEIGHT STATISTICS';
    lines{end+1}=repmat('=',1,72);
    lines{end+1}=''; lines{end+1}='INPUT PARAMETERS';
    lines{end+1}=field('Spectral significant wave height Hm0 (m)',sprintf('%.4f',r.Hm0));
    lines{end+1}=field('Still-water depth d (m)',sprintf('%.4f',r.d));
    slopeText = regexprep(sprintf('%.4f',r.slopeM),'0+$','');
    slopeText = regexprep(slopeText,'\.$','');
    lines{end+1}=field('Foreshore slope',['1:' slopeText]);
    lines{end+1}=field('tan(alpha)',sprintf('%.6f',r.tanAlpha));
    lines{end+1}=''; lines{end+1}='LOCAL WAVE PARAMETERS';
    lines{end+1}=field('Spectral variance m0 (m^2)',sprintf('%.9e',r.m0));
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
    for s=r.statistics
        if s.rayleigh_limit_applied, limit='YES'; else, limit='NO'; end
        lines{end+1}=statrow(sprintf('1/%d',s.N),sprintf('%.4f',s.threshold), ...
            sprintf('%.4f',s.mean_unlimited),sprintf('%.4f',s.rayleigh_mean), ...
            sprintf('%.4f',s.mean_final),limit);
    end
    lines{end+1}=''; lines{end+1}='FINAL NORMALIZED MEAN-HIGH-WAVE HEIGHTS';
    for s=r.statistics
        lines{end+1}=field(sprintf('H_(1/%d) / Hrms',s.N),sprintf('%.7f',s.mean_normalized));
    end
    lines{end+1}=''; lines{end+1}='RATIOS RELATIVE TO THE HIGHEST ONE THIRD';
    for s=r.statistics(2:end)
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
    text = strjoin(lines,sprintf('\n'));
end

function str = field(label,value,right)
    if nargin<3, right=true; end
    if right, str=sprintf('  %-42s : %16s',label,value);
    else, str=sprintf('  %-42s : %s',label,value); end
end

function str = statrow(fr,hn,raw,ray,fin,lim)
    str=sprintf('  %-10s %11s %12s %12s %12s %8s',fr,hn,raw,ray,fin,lim);
end
