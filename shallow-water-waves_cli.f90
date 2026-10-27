! Shallow-Water Individual Wave-Height Statistics
!
! Command-line calculation of local zero-crossing wave-height statistics with
! the Composite Weibull Distribution and Rayleigh limits. Input parameters:
! Hm0 (m), still-water depth d (m), foreshore slope denominator M in 1:M.
!
! Build:
!   gfortran -O2 -std=f2018 -Wall -Wextra -pedantic -ffree-line-length-none &
!     shallow-water-waves_cli.f90 -o shallow-water-waves_cli_f
!
! Use:
!   ./shallow-water-waves_cli_f 2.5 5 100
!   ./shallow-water-waves_cli_f 2.5 5 100 --output report.txt
!   ./shallow-water-waves_cli_f --no-file
!
module shallow_water_wave_statistics
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter :: dp = real64
    real(dp), parameter :: k1 = 2.0_dp, k2 = 3.6_dp, rayleigh_switch = 2.75_dp
    real(dp), parameter :: numerical_tol = 2.5e-14_dp
    integer, parameter :: max_gamma_iterations = 10000, max_root_iterations = 40
    integer, parameter :: n_values(6) = [3, 10, 50, 100, 250, 1000]
    real(dp), parameter :: rayleigh_mean_factors(6) = [ &
        1.001075736951740_dp, 1.272734273369137_dp, &
        1.560113379974762_dp, 1.668233372358517_dp, &
        1.801017222497626_dp, 1.984835590575388_dp]

    type :: statistic_t
        integer :: n = 0
        real(dp) :: threshold = 0.0_dp
        real(dp) :: rayleigh_threshold = 0.0_dp
        real(dp) :: mean_unlimited = 0.0_dp
        real(dp) :: rayleigh_mean = 0.0_dp
        real(dp) :: mean_final = 0.0_dp
        real(dp) :: mean_normalized = 0.0_dp
        logical :: rayleigh_limit_applied = .false.
    end type

    type, public :: analysis_t
        real(dp) :: hm0 = 0.0_dp, depth = 0.0_dp, slope_denominator = 0.0_dp
        real(dp) :: slope = 0.0_dp, m0 = 0.0_dp, hrms = 0.0_dp
        real(dp) :: saturation = 0.0_dp, htr = 0.0_dp, transition_normalized = 0.0_dp
        logical :: uses_cwd = .false., slope_outside_calibration = .false.
        logical :: saturation_outside_observed_range = .false.
        real(dp) :: scale1_normalized = 0.0_dp, scale2_normalized = 0.0_dp
        real(dp) :: transition_exponent = 0.0_dp, transition_exceedance_probability = 0.0_dp
        real(dp) :: transition_recurrence_count = 0.0_dp
        real(dp) :: moment_residual = 0.0_dp, continuity_residual = 0.0_dp
        integer :: solver_iterations = 0
        type(statistic_t) :: statistics(6)
    end type

    type :: report_builder_t
        character(len=:), allocatable :: content
        logical :: first = .true.
    end type

    public :: calculate, format_report

contains
    subroutine fail(message)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: message
        write(error_unit, '(a)') 'Error: '//message
        stop 1
    end subroutine

    subroutine require_positive(value, name)
        real(dp), intent(in) :: value
        character(len=*), intent(in) :: name
        if (.not. ieee_is_finite(value)) call fail(name//' must be finite and positive.')
        if (value <= 0.0_dp) call fail(name//' must be finite and positive.')
    end subroutine

    subroutine regularized_gamma(a, x, p, q)
        real(dp), intent(in) :: a, x
        real(dp), intent(out) :: p, q
        real(dp), parameter :: fpmin = 1.0e-300_dp, allowance = 2.0e-13_dp
        real(dp) :: prefactor, ap, delta, sum_value, b, c, d, fraction, an, nd
        integer :: i
        logical :: converged
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(x)) &
            call fail('Invalid argument in incomplete gamma function.')
        if (a <= 0.0_dp .or. x < 0.0_dp) call fail('Invalid argument in incomplete gamma function.')
        if (x <= 0.0_dp) then
            p = 0.0_dp
            q = 1.0_dp
            return
        end if
        prefactor = exp(a * log(x) - x - log_gamma(a))
        converged = .false.
        if (x < a + 1.0_dp) then
            ap = a
            delta = 1.0_dp / a
            sum_value = delta
            do i = 1, max_gamma_iterations
                ap = ap + 1.0_dp
                delta = delta * x / ap
                sum_value = sum_value + delta
                if (.not. ieee_is_finite(sum_value)) call fail('Gamma series became non-finite.')
                if (abs(delta) <= numerical_tol * abs(sum_value)) then
                    converged = .true.
                    exit
                end if
            end do
            if (.not. converged) call fail('Incomplete gamma series did not converge.')
            p = sum_value * prefactor
            q = 1.0_dp - p
        else
            b = x + 1.0_dp - a
            c = 1.0_dp / fpmin
            d = 1.0_dp / b
            fraction = d
            do i = 1, max_gamma_iterations
                nd = real(i, dp)
                an = -nd * (nd - a)
                b = b + 2.0_dp
                d = an * d + b
                if (abs(d) < fpmin) d = sign(fpmin, d)
                c = b + an / c
                if (abs(c) < fpmin) c = sign(fpmin, c)
                d = 1.0_dp / d
                delta = c * d
                fraction = fraction * delta
                if (.not. ieee_is_finite(fraction)) call fail('Gamma continued fraction became non-finite.')
                if (abs(delta - 1.0_dp) <= numerical_tol) then
                    converged = .true.
                    exit
                end if
            end do
            if (.not. converged) call fail('Incomplete gamma continued fraction did not converge.')
            q = fraction * prefactor
            p = 1.0_dp - q
        end if
        if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(q)) &
            call fail('Incomplete gamma result outside the probability range.')
        if (p < -allowance .or. p > 1.0_dp + allowance .or. &
            q < -allowance .or. q > 1.0_dp + allowance) &
            call fail('Incomplete gamma result outside the probability range.')
        p = max(0.0_dp, min(1.0_dp, p))
        q = max(0.0_dp, min(1.0_dp, q))
    end subroutine

    function gamma_lower(a, x) result(value)
        real(dp), intent(in) :: a, x
        real(dp) :: value, p, q
        call regularized_gamma(a, x, p, q)
        value = gamma(a) * p
    end function

    function gamma_upper(a, x) result(value)
        real(dp), intent(in) :: a, x
        real(dp) :: value, p, q
        call regularized_gamma(a, x, p, q)
        value = gamma(a) * q
    end function

    function gamma_interval(a, lower_bound, upper_bound) result(value)
        real(dp), intent(in) :: a, lower_bound, upper_bound
        real(dp) :: value, middle, half_width, displacement, u1, u2
        real(dp), parameter :: abscissae(8) = [ &
            0.09501250983763744_dp, 0.28160355077925891_dp, &
            0.45801677765722739_dp, 0.61787624440264375_dp, &
            0.75540440835500303_dp, 0.86563120238783174_dp, &
            0.94457502307323258_dp, 0.98940093499164993_dp ]
        real(dp), parameter :: weights(8) = [ &
            0.18945061045506850_dp, 0.18260341504492359_dp, &
            0.16915651939500254_dp, 0.14959598881657673_dp, &
            0.12462897125553387_dp, 0.09515851168249278_dp, &
            0.06225352393864789_dp, 0.02715245941175409_dp ]
        integer :: i
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(lower_bound) .or. &
            .not. ieee_is_finite(upper_bound)) call fail('Invalid incomplete gamma integration interval.')
        if (a <= 0.0_dp .or. lower_bound < 0.0_dp .or. upper_bound < lower_bound) &
            call fail('Invalid incomplete gamma integration interval.')
        if (upper_bound <= lower_bound) then
            value = 0.0_dp
            return
        end if
        middle = lower_bound + 0.5_dp * (upper_bound - lower_bound)
        half_width = 0.5_dp * (upper_bound - lower_bound)
        value = 0.0_dp
        do i = 1, 8
            displacement = half_width * abscissae(i)
            u1 = middle - displacement
            u2 = middle + displacement
            value = value + weights(i) * &
                (u1**(a - 1.0_dp) * exp(-u1) + u2**(a - 1.0_dp) * exp(-u2))
        end do
        value = value * half_width
        call require_positive(value, 'Gamma interval integral')
    end function

    ! Scalar Composite Weibull RMS closure in z = log(q).
    ! The derivative is analytic: the incomplete-gamma boundary terms cancel.
    subroutine evaluate_log_moment(log_q, log_t, residual, derivative)
        real(dp), intent(in) :: log_q, log_t
        real(dp), intent(out) :: residual, derivative
        real(dp) :: q, lower, upper, kernel
        q = exp(log_q)
        call require_positive(q, 'Transition exponent q')
        lower = exp(-(2.0_dp/k1)*log_q) * gamma_lower(1.0_dp + 2.0_dp/k1, q)
        upper = exp(-(2.0_dp/k2)*log_q) * gamma_upper(1.0_dp + 2.0_dp/k2, q)
        kernel = lower + upper
        call require_positive(kernel, 'Second-moment kernel')
        residual = 2.0_dp*log_t + log(kernel)
        derivative = -((2.0_dp/k1)*lower + (2.0_dp/k2)*upper)/kernel
        if (.not. ieee_is_finite(residual) .or. .not. ieee_is_finite(derivative)) &
            call fail('Invalid Newton residual or derivative.')
        if (derivative >= 0.0_dp) call fail('Invalid Newton residual or derivative.')
    end subroutine

    function initial_transition_exponent(t) result(q0)
        real(dp), intent(in) :: t
        real(dp) :: q0
        q0 = 1.002054709_dp*t*t*tanh(t**1.652259933_dp) * &
             exp(-tanh(tanh(t))*exp(-2.146235748_dp*t))
    end function

    function solve_transition_exponent(t, iterations) result(q)
        real(dp), intent(in) :: t
        integer, intent(out) :: iterations
        real(dp), parameter :: root_tol = 5.0e-13_dp
        real(dp) :: q, q0, z, log_t, residual, derivative, delta_z
        integer :: iter
        call require_positive(t, 'Normalized transition height')
        q0 = initial_transition_exponent(t)
        call require_positive(q0, 'Newton initial transition exponent')
        z = log(q0)
        log_t = log(t)
        do iter = 0, max_root_iterations
            call evaluate_log_moment(z, log_t, residual, derivative)
            if (abs(residual) <= root_tol) then
                iterations = iter
                q = exp(z)
                call require_positive(q, 'CWD transition exponent')
                return
            end if
            if (iter == max_root_iterations) exit
            delta_z = -residual/derivative
            z = z + delta_z
            if (.not. ieee_is_finite(z)) call fail('Non-finite Newton iteration.')
        end do
        call fail('CWD Newton solver did not converge.')
        q = 0.0_dp
        iterations = 0
    end function

    function cwd_second_moment(a1, a2, q) result(value)
        real(dp), intent(in) :: a1, a2, q
        real(dp) :: value
        value = a1*a1*gamma_lower(1.0_dp + 2.0_dp/k1, q) + &
                a2*a2*gamma_upper(1.0_dp + 2.0_dp/k2, q)
    end function

    function cwd_threshold_normalized(n, a1, a2, q) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: a1, a2, q
        real(dp) :: value, log_n
        log_n = log(real(n,dp))
        if (log_n < q) then
            value = a1 * log_n**(1.0_dp/k1)
        else
            value = a2 * log_n**(1.0_dp/k2)
        end if
    end function

    function cwd_mean_highest_normalized(n, a1, a2, q) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: a1, a2, q
        real(dp) :: value, log_n, tail_moment
        log_n = log(real(n,dp))
        if (log_n < q) then
            tail_moment = a1 * gamma_interval(1.0_dp + 1.0_dp/k1, log_n, q) + &
                          a2 * gamma_upper(1.0_dp + 1.0_dp/k2, q)
        else
            tail_moment = a2 * gamma_upper(1.0_dp + 1.0_dp/k2, log_n)
        end if
        value = real(n,dp) * tail_moment
        call require_positive(value, 'CWD mean-high-wave height')
    end function

    function calculate(hm0, depth, slope_denominator) result(r)
        real(dp), intent(in) :: hm0, depth, slope_denominator
        type(analysis_t) :: r
        real(dp) :: sigma, t, q, m2, e1, e2, rayleigh_scale, log_n
        integer :: i
        call require_positive(hm0, 'Hm0')
        call require_positive(depth, 'Water depth')
        call require_positive(slope_denominator, 'Slope denominator M')
        r%hm0 = hm0
        r%depth = depth
        r%slope_denominator = slope_denominator
        r%slope = 1.0_dp / slope_denominator
        call require_positive(r%slope, 'Foreshore slope')
        r%slope_outside_calibration = slope_denominator < 20.0_dp .or. slope_denominator > 250.0_dp
        sigma = hm0/4.0_dp
        call require_positive(sigma, 'Spectral standard deviation')
        r%m0 = sigma*sigma
        call require_positive(r%m0, 'Spectral variance m0')
        r%saturation = sigma/depth
        call require_positive(r%saturation, 'Degree of saturation')
        r%saturation_outside_observed_range = &
            r%saturation < 0.025_dp .or. r%saturation > 0.16_dp
        r%hrms = (2.69_dp + 3.24_dp * r%saturation) * sigma
        call require_positive(r%hrms, 'Root-mean-square wave height')
        r%htr = (0.35_dp + 5.8_dp*r%slope) * depth
        call require_positive(r%htr, 'Transitional wave height')
        r%transition_normalized = r%htr/r%hrms
        call require_positive(r%transition_normalized, 'Normalized transition height')
        r%uses_cwd = r%transition_normalized <= rayleigh_switch
        q = 0.0_dp
        if (r%uses_cwd) then
            t = r%transition_normalized
            q = solve_transition_exponent(t, r%solver_iterations)
            r%transition_exponent = q
            r%scale1_normalized = t * exp(-log(q)/k1)
            r%scale2_normalized = t * exp(-log(q)/k2)
            call require_positive(r%scale1_normalized, 'First CWD scale')
            call require_positive(r%scale2_normalized, 'Second CWD scale')
            r%transition_exceedance_probability = exp(-q)
            r%transition_recurrence_count = exp(q)
            m2 = cwd_second_moment(r%scale1_normalized, r%scale2_normalized, q)
            r%moment_residual = sqrt(m2) - 1.0_dp
            e1 = (t/r%scale1_normalized)**k1
            e2 = (t/r%scale2_normalized)**k2
            r%continuity_residual = e1-e2
            if (.not. ieee_is_finite(m2)) call fail('CWD probability constraints were not satisfied.')
            if (abs(r%moment_residual) > 1.0e-10_dp .or. &
                abs(r%continuity_residual) > 1.0e-10_dp) &
                call fail('CWD probability constraints were not satisfied.')
        end if
        rayleigh_scale = hm0/sqrt(2.0_dp)
        call require_positive(rayleigh_scale, 'Rayleigh height scale')
        do i = 1, size(n_values)
            r%statistics(i)%n = n_values(i)
            log_n = log(real(n_values(i),dp))
            r%statistics(i)%rayleigh_threshold = rayleigh_scale*sqrt(log_n)
            r%statistics(i)%rayleigh_mean = rayleigh_mean_factors(i)*hm0
            call require_positive(r%statistics(i)%rayleigh_threshold, 'Rayleigh exceedance threshold')
            call require_positive(r%statistics(i)%rayleigh_mean, 'Rayleigh high-wave mean')
            if (r%uses_cwd) then
                r%statistics(i)%threshold = cwd_threshold_normalized(n_values(i), &
                    r%scale1_normalized, r%scale2_normalized, q) * r%hrms
                r%statistics(i)%mean_unlimited = cwd_mean_highest_normalized(n_values(i), &
                    r%scale1_normalized, r%scale2_normalized, q) * r%hrms
                call require_positive(r%statistics(i)%threshold, 'CWD exceedance threshold')
                call require_positive(r%statistics(i)%mean_unlimited, 'CWD high-wave mean')
                r%statistics(i)%mean_final = min(r%statistics(i)%mean_unlimited, r%statistics(i)%rayleigh_mean)
                r%statistics(i)%rayleigh_limit_applied = &
                    r%statistics(i)%mean_unlimited > r%statistics(i)%rayleigh_mean * &
                    (1.0_dp + 8.0_dp*epsilon(1.0_dp))
            else
                r%statistics(i)%threshold = r%statistics(i)%rayleigh_threshold
                r%statistics(i)%mean_unlimited = r%statistics(i)%rayleigh_mean
                r%statistics(i)%mean_final = r%statistics(i)%rayleigh_mean
            end if
            r%statistics(i)%mean_normalized = r%statistics(i)%mean_final / r%hrms
            call require_positive(r%statistics(i)%mean_normalized, 'Normalized high-wave mean')
            if (.not. (r%statistics(i)%mean_unlimited > r%statistics(i)%threshold)) &
                call fail('Inconsistent wave-height statistics.')
            if (r%statistics(i)%mean_final > r%statistics(i)%rayleigh_mean*(1.0_dp+1.0e-12_dp)) &
                call fail('Inconsistent wave-height statistics.')
            if (i > 1) call verify_statistic_order(r%statistics(max(1,i-1)), r%statistics(i))
        end do
    end function

    subroutine verify_statistic_order(previous, current)
        type(statistic_t), intent(in) :: previous, current
        if (.not. (current%threshold > previous%threshold) .or. &
            .not. (current%mean_unlimited > previous%mean_unlimited) .or. &
            .not. (current%mean_final > previous%mean_final)) &
            call fail('Wave-height ordering constraint not satisfied.')
    end subroutine

    function int_string(value) result(s)
        integer, intent(in) :: value
        character(len=:), allocatable :: s
        character(len=24) :: tmp
        write(tmp,'(i0)') value
        s = trim(tmp)
    end function

    function report_number(value, decimals, scientific) result(s)
        real(dp), intent(in) :: value
        integer, intent(in), optional :: decimals
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: s
        character(len=160) :: fmt, buffer
        integer :: nd, i
        logical :: as_scientific
        nd = 6
        if (present(decimals)) nd = decimals
        as_scientific = .false.
        if (present(scientific)) as_scientific = scientific
        if (as_scientific) then
            write(fmt,'("(ES160.",i0,"E2)")') nd
        else
            write(fmt,'("(F160.",i0,")")') nd
        end if
        write(buffer,trim(fmt)) value
        s = trim(adjustl(buffer))
        if (as_scientific) then
            do i=1,len(s)
                if (s(i:i) == 'E') s(i:i) = 'e'
            end do
        end if
    end function

    function report_slope(value) result(s)
        real(dp), intent(in) :: value
        character(len=:), allocatable :: s
        character(len=:), allocatable :: number
        integer :: n
        number = report_number(value,4)
        if (index(number,'.') > 0) then
            n = len(number)
            do while (n > 0)
                if (number(n:n) /= '0') exit
                n = n-1
            end do
            if (number(n:n) == '.') n = n-1
            number = number(:n)
        end if
        s = '1:'//number
    end function

    function pad_right(s, width) result(out)
        character(len=*), intent(in) :: s
        integer, intent(in) :: width
        character(len=:), allocatable :: out
        out = s // repeat(' ', max(0,width-len(s)))
    end function

    function pad_left(s, width) result(out)
        character(len=*), intent(in) :: s
        integer, intent(in) :: width
        character(len=:), allocatable :: out
        out = repeat(' ', max(0,width-len(s))) // s
    end function

    subroutine line_add(builder, line)
        type(report_builder_t), intent(inout) :: builder
        character(len=*), intent(in) :: line
        if (builder%first) then
            builder%content = line
            builder%first = .false.
        else
            builder%content = builder%content // achar(10) // line
        end if
    end subroutine

    subroutine field_add(builder, label, value, right_align)
        type(report_builder_t), intent(inout) :: builder
        character(len=*), intent(in) :: label, value
        logical, intent(in), optional :: right_align
        logical :: right
        right = .true.
        if (present(right_align)) right = right_align
        if (right) then
            call line_add(builder, '  '//pad_right(label,42)//' : '//pad_left(value,16))
        else
            call line_add(builder, '  '//pad_right(label,42)//' : '//value)
        end if
    end subroutine

    subroutine stat_row_add(builder, fraction, threshold, raw_mean, rayleigh_mean, final_mean, lim)
        type(report_builder_t), intent(inout) :: builder
        character(len=*), intent(in) :: fraction, threshold, raw_mean, rayleigh_mean, final_mean, lim
        call line_add(builder, '  '//pad_right(fraction,10)//' '//pad_left(threshold,11)// &
            ' '//pad_left(raw_mean,12)//' '//pad_left(rayleigh_mean,12)// &
            ' '//pad_left(final_mean,12)//' '//pad_left(lim,8))
    end subroutine

    function format_report(r) result(report)
        type(analysis_t), intent(in) :: r
        character(len=:), allocatable :: report, ns, lim
        type(report_builder_t) :: b
        integer :: i
        call line_add(b, 'SHALLOW-WATER INDIVIDUAL WAVE-HEIGHT STATISTICS')
        call line_add(b, repeat('=',72))
        call line_add(b, '')
        call line_add(b, 'INPUT PARAMETERS')
        call field_add(b, 'Spectral significant wave height Hm0 (m)', report_number(r%hm0,4))
        call field_add(b, 'Still-water depth d (m)', report_number(r%depth,4))
        call field_add(b, 'Foreshore slope', report_slope(r%slope_denominator))
        call field_add(b, 'tan(alpha)', report_number(r%slope))
        call line_add(b, '')
        call line_add(b, 'LOCAL WAVE PARAMETERS')
        call field_add(b, 'Spectral variance m0 (m^2)', report_number(r%m0))
        call field_add(b, 'Individual-wave Hrms (m)', report_number(r%hrms))
        call field_add(b, 'Degree of saturation sqrt(m0)/d', report_number(r%saturation))
        call field_add(b, 'Transitional wave height Htr (m)', report_number(r%htr))
        call field_add(b, 'Normalized transition Htr/Hrms', report_number(r%transition_normalized))
        if (r%uses_cwd) then
            call field_add(b, 'Distribution', 'Composite Weibull (B&G)', .false.)
        else
            call field_add(b, 'Distribution', 'Rayleigh', .false.)
        end if
        call line_add(b, '')
        call line_add(b, 'DISTRIBUTION PARAMETERS')
        if (r%uses_cwd) then
            call field_add(b, 'Normalized first Weibull scale H1/Hrms', report_number(r%scale1_normalized))
            call field_add(b, 'Normalized second Weibull scale H2/Hrms', report_number(r%scale2_normalized))
            call field_add(b, 'Dimensional Weibull scale H1 (m)', report_number(r%scale1_normalized*r%hrms))
            call field_add(b, 'Dimensional Weibull scale H2 (m)', report_number(r%scale2_normalized*r%hrms))
            call field_add(b, 'Transition exponent q', report_number(r%transition_exponent))
            call field_add(b, 'Probability P(H > Htr)', report_number(r%transition_exceedance_probability,9))
            call field_add(b, 'Equivalent transition recurrence Ntr', &
                report_number(r%transition_recurrence_count)//' waves')
        else
            call field_add(b, 'Composite Weibull scales H1 and H2', 'N/A', .false.)
            call field_add(b, 'Composite Weibull transition exponent', 'N/A', .false.)
        end if
        call line_add(b, '')
        call line_add(b, 'WAVE-HEIGHT STATISTICS (m)')
        call stat_row_add(b, 'HIGHEST', 'H_N', 'RAW MEAN', 'RAYLEIGH', 'FINAL MEAN', 'LIMIT')
        call stat_row_add(b, 'FRACTION', 'THRESHOLD', 'H_(1/N)', 'H_(1/N)', 'H_(1/N)', 'APPLIED')
        call line_add(b, '  '//repeat('-',70))
        do i=1,6
            ns = int_string(r%statistics(i)%n)
            lim = 'NO'
            if (r%statistics(i)%rayleigh_limit_applied) lim='YES'
            call stat_row_add(b, '1/'//ns, report_number(r%statistics(i)%threshold,4), &
                report_number(r%statistics(i)%mean_unlimited,4), &
                report_number(r%statistics(i)%rayleigh_mean,4), &
                report_number(r%statistics(i)%mean_final,4), lim)
        end do
        call line_add(b, '')
        call line_add(b, 'FINAL NORMALIZED MEAN-HIGH-WAVE HEIGHTS')
        do i=1,6
            ns=int_string(r%statistics(i)%n)
            call field_add(b, 'H_(1/'//ns//') / Hrms', report_number(r%statistics(i)%mean_normalized,7))
        end do
        call line_add(b, '')
        call line_add(b, 'RATIOS RELATIVE TO THE HIGHEST ONE THIRD')
        do i=2,6
            ns=int_string(r%statistics(i)%n)
            call field_add(b, 'H_(1/'//ns//') / H_(1/3)', &
                report_number(r%statistics(i)%mean_final/r%statistics(1)%mean_final))
        end do
        call line_add(b, '')
        call line_add(b, 'NUMERICAL DIAGNOSTICS')
        if (r%uses_cwd) then
            call field_add(b, 'Scalar root iterations', int_string(r%solver_iterations))
            call field_add(b, 'Second-moment residual', report_number(r%moment_residual,10,.true.))
            call field_add(b, 'Continuity exponent residual', report_number(r%continuity_residual,10,.true.))
        else
            call field_add(b, 'CWD nonlinear solver', 'N/A', .false.)
        end if
        call line_add(b, '')
        call line_add(b, 'INTERPRETATION')
        call line_add(b, '  H_N is the height exceeded with probability 1/N.')
        call line_add(b, '  H_(1/N) is the mean height of the highest 1/N fraction.')
        call line_add(b, '  Both are individual-wave statistics, not storm maxima.')
        call line_add(b, '  Mean-high-wave Rayleigh limits are applied separately.')
        call line_add(b, '  H_N always refers to the selected underlying distribution.')
        call line_add(b, '  The program uses local wave conditions; it does not propagate waves.')
        if (r%slope_outside_calibration .or. &
            (r%uses_cwd .and. r%saturation_outside_observed_range)) then
            call line_add(b, '')
            call line_add(b, 'APPLICABILITY NOTICE')
            if (r%slope_outside_calibration) then
                call line_add(b, '  Foreshore slope is outside the laboratory slope range')
                call line_add(b, '  1:20 to 1:250 of the Composite Weibull formulation.')
                call line_add(b, '  Predictions for an extensive horizontal, constant-depth')
                call line_add(b, '  seabed require a separate applicability assessment.')
            end if
            if (r%uses_cwd .and. r%saturation_outside_observed_range) then
                call line_add(b, '  The degree of saturation is outside the approximate')
                call line_add(b, '  interval represented by the original foreshore tests.')
                call line_add(b, '  Interpretation requires engineering assessment.')
            end if
        end if
        report = b%content
    end function
end module shallow_water_wave_statistics

program shallow_water_waves_cli
    use iso_fortran_env, only: real64, output_unit, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use shallow_water_wave_statistics, only: analysis_t, calculate, format_report
    implicit none
    integer, parameter :: dp = real64
    type(analysis_t) :: result
    character(len=4096) :: arg, pos(3), output_path
    character(len=:), allocatable :: report
    real(dp) :: hm0, depth, slope_m
    integer :: argc, i, positional_count, ios, unit
    logical :: write_file, show_help
    character(len=512) :: iomsg

    argc = command_argument_count()
    write_file = .true.
    show_help = .false.
    output_path = 'report.txt'
    positional_count = 0
    i = 1
    do while (i <= argc)
        arg = ''
        call get_command_argument(i, arg, status=ios)
        if (ios /= 0) call die('Unable to read command-line argument.')
        select case(trim(arg))
        case('--help', '-h')
            show_help = .true.
        case('--no-file')
            write_file = .false.
        case('--output', '-o')
            i = i+1
            if (i > argc) call die('Missing path after --output.')
            call get_command_argument(i, output_path, status=ios)
            if (ios /= 0) call die('Output path is too long.')
            if (len_trim(output_path) == 0) call die('Output path must not be empty.')
        case default
            if (index(trim(arg), '--') == 1) call die('Unrecognized option: '//trim(arg))
            positional_count = positional_count+1
            if (positional_count > 3) call die('Provide exactly three numerical inputs or none.')
            pos(positional_count) = arg
        end select
        i = i+1
    end do

    if (show_help) then
        call usage()
        stop
    end if
    if (positional_count /= 0 .and. positional_count /= 3) &
        call die('Provide exactly three numerical inputs or none.')

    if (positional_count == 3) then
        hm0 = parse_positive_number(pos(1), 'Hm0')
        depth = parse_positive_number(pos(2), 'water depth')
        slope_m = parse_positive_number(pos(3), 'slope denominator')
    else
        write(output_unit,'(a)',advance='no') 'Hm0 (m): '
        hm0 = prompt_value('Hm0')
        write(output_unit,'(a)',advance='no') 'Water depth d (m): '
        depth = prompt_value('water depth')
        write(output_unit,'(a)',advance='no') 'Foreshore slope 1:M, enter M: '
        slope_m = prompt_value('slope denominator')
        write(output_unit,'(a)') ''
    end if

    result = calculate(hm0, depth, slope_m)
    report = format_report(result)
    if (write_file) then
        open(newunit=unit, file=trim(output_path), status='replace', &
            action='write', access='stream', form='unformatted', iostat=ios, iomsg=iomsg)
        if (ios /= 0) call die('Cannot open output file: '//trim(output_path))
        write(unit,iostat=ios,iomsg=iomsg) report
        if (ios /= 0) call die('Cannot write report completely: '//trim(output_path))
        close(unit,iostat=ios,iomsg=iomsg)
        if (ios /= 0) call die('Cannot close output file: '//trim(output_path))
    end if
    write(output_unit,'(a)') report

contains
    subroutine die(message)
        character(len=*), intent(in) :: message
        write(error_unit,'(a)') 'Error: '//message
        stop 1
    end subroutine

    subroutine usage()
        write(output_unit,'(a)') &
            'Usage:', &
            '  shallow-water-waves_cli_f Hm0 depth slope_M [--output path] [--no-file]', &
            '  shallow-water-waves_cli_f [--output path] [--no-file]', &
            '  shallow-water-waves_cli_f --help', &
            '', &
            'Hm0       Local spectral significant wave height (m).', &
            'depth     Local still-water depth (m).', &
            'slope_M   Positive denominator M of a 1:M foreshore slope.', &
            '', &
            'With no numeric arguments, the program prompts for inputs.', &
            'The report is displayed on the console and saved as report.txt', &
            'unless --no-file is specified. --output selects another file.'
    end subroutine

    function prompt_value(name) result(value)
        character(len=*), intent(in) :: name
        real(dp) :: value
        character(len=4096) :: text
        integer :: stat
        read(*,'(a)',iostat=stat) text
        if (stat /= 0) call die('No input available for '//name//'.')
        value = parse_positive_number(text,name)
    end function

    function parse_positive_number(raw, name) result(value)
        character(len=*), intent(in) :: raw, name
        real(dp) :: value
        character(len=:), allocatable :: s
        integer :: j, code, dot_count, digit_count, exponent_digits, exp_pos, local_status
        s = trim(adjustl(raw))
        if (len(s) == 0) call die('Invalid numeric value for '//name//'.')
        do j=1,len(s)
            if (s(j:j) == ',') s(j:j) = '.'
        end do
        exp_pos=0
        do j=1,len(s)
            if (s(j:j) == 'e' .or. s(j:j)=='E') then
                if (exp_pos /= 0) call die('Invalid numeric value for '//name//'.')
                exp_pos=j
            end if
        end do
        dot_count=0; digit_count=0
        do j=1,merge(exp_pos-1,len(s),exp_pos>0)
            if (s(j:j) == '+' .or. s(j:j) == '-') then
                if (j/=1) call die('Invalid numeric value for '//name//'.')
            else if (s(j:j) == '.') then
                dot_count=dot_count+1
                if (dot_count>1) call die('Invalid numeric value for '//name//'.')
            else
                code=iachar(s(j:j))
                if (code < iachar('0') .or. code>iachar('9')) &
                    call die('Invalid numeric value for '//name//'.')
                digit_count=digit_count+1
            end if
        end do
        if (digit_count==0) call die('Invalid numeric value for '//name//'.')
        if (exp_pos>0) then
            exponent_digits=0
            do j=exp_pos+1,len(s)
                if (s(j:j)=='+' .or. s(j:j)=='-') then
                    if (j/=exp_pos+1) call die('Invalid numeric value for '//name//'.')
                else
                    code=iachar(s(j:j))
                    if (code<iachar('0') .or. code>iachar('9')) &
                        call die('Invalid numeric value for '//name//'.')
                    exponent_digits=exponent_digits+1
                end if
            end do
            if (exponent_digits==0) call die('Invalid numeric value for '//name//'.')
        end if
        read(s,*,iostat=local_status) value
        if (local_status/=0) call die('Invalid numeric value for '//name//'.')
        if (.not. ieee_is_finite(value)) call die('Invalid numeric value for '//name//'.')
        if (value <= 0.0_dp) call die('Invalid numeric value for '//name//'.')
    end function
end program shallow_water_waves_cli
