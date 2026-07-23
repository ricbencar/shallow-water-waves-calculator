!***********************************************************************
! Program: shallow-water-waves_cli.f90
!
! Detailed Description:
! This program computes local shallow-foreshore wave-height distribution
! parameters using a model based on the Composed Weibull distribution.
!
! The command-line application performs the following:
! 1. If three command-line arguments are provided, they are used as
!    Hm0 (local significant spectral wave height), d (local water depth),
!    and slopeM (beach slope 1:m). Otherwise, the program
!    prompts the user for these values.
! 2. Computes the following intermediate values:
!    - Free-surface variance (m0) is calculated directly using m0 = (Hm0/4)^2.
!    - Root-mean-square wave height (Hrms) is then calculated from m0
!      using the equation: Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0).
!    - A dimensional transitional wave height: Htr = (0.35 + 5.8*(1/m)) * d.
!    - The dimensionless transitional parameter: H̃_tr = Htr / Hrms.
!
! 3. Calculates the dimensionless wave-height ratios (Hᵢ/Hrms)
!    by solving a system of non-linear equations derived from the Composite Weibull
!    distribution, ensuring the normalized Hrms of the distribution equals one.
!    This involves using a Newton-Raphson matrix method for simultaneous root-finding
!    and functions for unnormalized incomplete gamma calculations. These ratios are then
!    converted to dimensional quantities (in meters) by multiplying with Hrms.
!
! 4. A detailed report is then generated (and written to "report.txt") with
!    the input parameters, intermediate values, calculated ratios and computed
!    dimensional wave heights, as well as diagnostic ratios.
!
! Overshoot-prevention logic:
! To prevent the composite Weibull distribution from predicting wave heights
! greater than the theoretical Rayleigh distribution (which is valid for deep
! water), the following logic is applied:
!
! a) H_tr Threshold Switch: If the dimensionless transitional height (Htr/Hrms)
!    is greater than 2.75, the conditions are considered deep-water dominant.
!    The program bypasses the B&G calculation and uses established Rayleigh
!    distribution values for all H1/N statistics.
!
! b) Capping Statistical Parameters: If the B&G calculation is performed, the
!    resulting dimensional values for all H1/N are capped at their
!    theoretical Rayleigh limits (e.g., H1/3 <= 1.001075*Hm0, H1/10 <= 1.272734*Hm0, etc.).
!
! Compilation Instructions (example using gfortran):
!
! gfortran -O3 -march=native -std=f2018 -Wall -Wextra -pedantic -Wconversion -static -fno-underscoring -o shallow-water-waves_cli shallow-water-waves_cli.f90
!
! To run with command-line arguments (e.g., Hm0=2.5, d=5, slopeM=100):
! ./shallow-water-waves_cli 2.5 5 100
!***********************************************************************
PROGRAM ShallowWaterWaves
    IMPLICIT NONE

    ! Define a portable double precision kind
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

    ! Global parameters for the Composite Weibull distribution.
    REAL(KIND=dp), PARAMETER :: K1 = 2.0_dp ! Exponent for the first part of the Composite Weibull distribution (Rayleigh-shaped)
    REAL(KIND=dp), PARAMETER :: K2 = 3.6_dp ! Exponent for the second part of the Composite Weibull distribution

    ! Precision for Numerical Methods.
    REAL(KIND=dp), PARAMETER :: EPSILON = 1.0E-12_dp     ! Max allowable error for solver convergence.
    REAL(KIND=dp), PARAMETER :: JACOBIAN_DX = 1.0E-8_dp  ! Step size for finite difference Jacobian.
    REAL(KIND=dp), PARAMETER :: LOCAL_EPS = 1.0E-16_dp   ! Tolerance for gamma function convergence.

    ! File unit for the report
    INTEGER, PARAMETER :: REPORT_UNIT = 10

    ! Derived type to hold all results for cleaner data management.
    TYPE :: wave_results_t
        REAL(KIND=dp) :: Hm0 = 0.0_dp
        REAL(KIND=dp) :: d = 0.0_dp
        REAL(KIND=dp) :: slopeM = 0.0_dp
        CHARACTER(LEN=10) :: distribution_type = ''
        REAL(KIND=dp) :: Hrms = 0.0_dp
        REAL(KIND=dp) :: m0 = 0.0_dp
        REAL(KIND=dp) :: tanAlpha = 0.0_dp
        REAL(KIND=dp) :: Htr_dim = 0.0_dp
        REAL(KIND=dp) :: Htr_tilde = 0.0_dp
        REAL(KIND=dp) :: H1_Hrms = 0.0_dp
        REAL(KIND=dp) :: H2_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_3_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_10_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_50_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_100_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_250_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_1000_Hrms = 0.0_dp
        REAL(KIND=dp) :: H1_dim = 0.0_dp
        REAL(KIND=dp) :: H2_dim = 0.0_dp
        REAL(KIND=dp) :: H1_3_dim = 0.0_dp
        REAL(KIND=dp) :: H1_10_dim = 0.0_dp
        REAL(KIND=dp) :: H1_50_dim = 0.0_dp
        REAL(KIND=dp) :: H1_100_dim = 0.0_dp
        REAL(KIND=dp) :: H1_250_dim = 0.0_dp
        REAL(KIND=dp) :: H1_1000_dim = 0.0_dp
        REAL(KIND=dp) :: ratio_1_10_div_1_3 = 0.0_dp
        REAL(KIND=dp) :: ratio_1_50_div_1_3 = 0.0_dp
        REAL(KIND=dp) :: ratio_1_100_div_1_3 = 0.0_dp
        REAL(KIND=dp) :: ratio_1_250_div_1_3 = 0.0_dp
        REAL(KIND=dp) :: ratio_1_1000_div_1_3 = 0.0_dp
    END TYPE wave_results_t

    REAL(KIND=dp) :: Hm0_in, d_in, slopeM_in
    TYPE(wave_results_t) :: results
    LOGICAL :: success

    CALL get_user_input(Hm0_in, d_in, slopeM_in)

    results%Hm0 = Hm0_in
    results%d = d_in
    results%slopeM = slopeM_in

    ! Perform calculations
    success = perform_wave_analysis(results)

    ! Check status and generate report
    IF (success) THEN
        CALL generate_and_write_report(results)
    ELSE
        WRITE(*,*) "ERROR: A calculation failed."
        OPEN(UNIT=REPORT_UNIT, FILE="report.txt", STATUS="REPLACE", ACTION="WRITE")
        WRITE(REPORT_UNIT, '(A)') "ERROR: A calculation failed."
        CLOSE(REPORT_UNIT)
        ERROR STOP 1
    END IF

CONTAINS

    ! =======================================================================
    ! Gamma Functions
    ! =======================================================================

    FUNCTION incomplete_gamma_lower(a, x) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: a, x
        REAL(KIND=dp) :: res
        INTEGER, PARAMETER :: MAXIT = 500
        REAL(KIND=dp) :: gln, ap, sum_val, del_val, b, c, d_val, h, an, del, exp_term
        INTEGER :: n_iter

        gln = LOG_GAMMA(a)

        IF (x < a + 1.0_dp) THEN
            ap = a
            sum_val = 1.0_dp / a
            del_val = sum_val
            DO n_iter = 1, MAXIT
                ap = ap + 1.0_dp
                del_val = del_val * x / ap
                sum_val = sum_val + del_val
                IF (ABS(del_val) < ABS(sum_val) * LOCAL_EPS) THEN
                    exp_term = EXP(-x + a * LOG(x) - gln)
                    res = sum_val * exp_term * GAMMA(a)
                    RETURN
                END IF
            END DO
            res = HUGE(1.0_dp)
        ELSE
            b = x + 1.0_dp - a
            c = 1.0_dp / TINY(1.0_dp)
            d_val = 1.0_dp / b
            h = d_val
            DO n_iter = 1, MAXIT
                an = -1.0_dp * REAL(n_iter, KIND=dp) * (REAL(n_iter, KIND=dp) - a)
                b = b + 2.0_dp
                d_val = an * d_val + b
                IF (ABS(d_val) < TINY(1.0_dp)) d_val = TINY(1.0_dp)
                c = b + an / c
                IF (ABS(c) < TINY(1.0_dp)) c = TINY(1.0_dp)
                d_val = 1.0_dp / d_val
                del = d_val * c
                h = h * del
                IF (ABS(del - 1.0_dp) < LOCAL_EPS) EXIT
            END DO
            exp_term = EXP(-x + a * LOG(x) - gln)
            res = (1.0_dp - exp_term * h) * GAMMA(a)
        END IF
    END FUNCTION incomplete_gamma_lower

    FUNCTION incomplete_gamma_upper(a, x) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: a, x
        REAL(KIND=dp) :: res
        res = GAMMA(a) - incomplete_gamma_lower(a, x)
    END FUNCTION incomplete_gamma_upper

    ! =======================================================================
    ! Core Wave Model Functions
    ! =======================================================================

    FUNCTION F1(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val
        REAL(KIND=dp) :: res, u, v, term1, term2, sum_terms
        IF (H1_Hrms_val <= 0.0_dp .OR. H2_Hrms_val <= 0.0_dp) THEN
            res = HUGE(1.0_dp); RETURN
        END IF
        u = (Htr_Hrms_val / H1_Hrms_val)**K1
        v = (Htr_Hrms_val / H2_Hrms_val)**K2
        term1 = H1_Hrms_val**2 * incomplete_gamma_lower(2.0_dp / K1 + 1.0_dp, u)
        term2 = H2_Hrms_val**2 * incomplete_gamma_upper(2.0_dp / K2 + 1.0_dp, v)
        sum_terms = term1 + term2
        res = SQRT(MAX(0.0_dp, sum_terms)) - 1.0_dp
    END FUNCTION F1

    FUNCTION F2(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val
        REAL(KIND=dp) :: res
        IF (H1_Hrms_val <= 0.0_dp .OR. H2_Hrms_val <= 0.0_dp) THEN
            res = HUGE(1.0_dp); RETURN
        END IF
        res = (Htr_Hrms_val / H1_Hrms_val)**K1 - (Htr_Hrms_val / H2_Hrms_val)**K2
    END FUNCTION F2

    SUBROUTINE solve_linear_system_2x2(J11, J12, J21, J22, b1, b2, dx1, dx2)
        REAL(KIND=dp), INTENT(IN) :: J11, J12, J21, J22, b1, b2
        REAL(KIND=dp), INTENT(OUT) :: dx1, dx2
        REAL(KIND=dp) :: determinant
        determinant = J11 * J22 - J12 * J21
        IF (ABS(determinant) < 1.0E-20_dp) THEN
            dx1 = 0.0_dp; dx2 = 0.0_dp; RETURN
        END IF
        dx1 = (b1 * J22 - b2 * J12) / determinant
        dx2 = (J11 * b2 - J21 * b1) / determinant
    END SUBROUTINE solve_linear_system_2x2

    SUBROUTINE get_initial_guesses(Htr_Hrms_val, H1_initial, H2_initial)
        REAL(KIND=dp), INTENT(IN) :: Htr_Hrms_val
        REAL(KIND=dp), INTENT(OUT) :: H1_initial, H2_initial
        H1_initial = 2.244660800090239e-03_dp + &
             tanh(1.918610494219390e+00_dp * Htr_Hrms_val)**1.780892753373355e-01_dp / &
             tanh(sinh(1.009497360864962e+00_dp * Htr_Hrms_val))**9.777939607559606e-01_dp
        H2_initial = 1.059259665431797_dp + (0.2059286860468916_dp * Htr_Hrms_val) / &
                     (1.0_dp + 3.865701948059343_dp * Htr_Hrms_val**(-3.479682433107255_dp))
        IF (H1_initial <= 0.0_dp) H1_initial = TINY(1.0_dp)
        IF (H2_initial <= 0.0_dp) H2_initial = TINY(1.0_dp)
    END SUBROUTINE get_initial_guesses

    FUNCTION newtonRaphsonSystemSolver(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out, maxit) RESULT(converged)
        REAL(KIND=dp), INTENT(IN) :: Htr_Hrms_val
        REAL(KIND=dp), INTENT(INOUT) :: H1_Hrms_out, H2_Hrms_out
        INTEGER, INTENT(IN) :: maxit
        LOGICAL :: converged
        REAL(KIND=dp) :: f1_val, f2_val, J11, J12, J21, J22, dH1, dH2
        INTEGER :: iter

        CALL get_initial_guesses(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out)

        DO iter = 0, maxit - 1
            f1_val = F1(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)
            f2_val = F2(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)
            IF (ABS(f1_val) < EPSILON .AND. ABS(f2_val) < EPSILON) THEN
                converged = .TRUE.; RETURN
            END IF
            J11 = (F1(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F1(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)
            J12 = (F1(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F1(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)
            J21 = (F2(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F2(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)
            J22 = (F2(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F2(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)
            CALL solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2)
            H1_Hrms_out = H1_Hrms_out + dH1
            H2_Hrms_out = H2_Hrms_out + dH2
            IF (H1_Hrms_out <= 0.0_dp) H1_Hrms_out = TINY(1.0_dp)
            IF (H2_Hrms_out <= 0.0_dp) H2_Hrms_out = TINY(1.0_dp)
        END DO
        WRITE(*,*) "Newton-Raphson solver failed to converge for Htr_Hrms =", Htr_Hrms_val
        converged = .FALSE.
    END FUNCTION newtonRaphsonSystemSolver

    FUNCTION calculate_HN(N, H1, H2, Htr) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: N, H1, H2, Htr
        REAL(KIND=dp) :: res, HN_candidate1
        HN_candidate1 = H1 * (LOG(N))**(1.0_dp / K1)
        IF (HN_candidate1 < Htr - EPSILON) THEN
            res = HN_candidate1
        ELSE
            res = H2 * (LOG(N))**(1.0_dp / K2)
        END IF
    END FUNCTION calculate_HN

    FUNCTION calculate_H1N(N_val, H1, H2, Htr) RESULT(res)
        REAL(KIND=dp), INTENT(IN) :: N_val, H1, H2, Htr
        REAL(KIND=dp) :: res, H_N_val, term1_a, term1_x_ln_Nval, term1_x_HtrH1, &
                         gamma1, term2_a, term2_x_HtrH2, gamma2
        H_N_val = calculate_HN(N_val, H1, H2, Htr)
        term1_x_ln_Nval = LOG(N_val)
        term2_a = 1.0_dp / K2 + 1.0_dp
        IF (H_N_val < Htr - EPSILON) THEN
            term1_a = 1.0_dp / K1 + 1.0_dp
            term1_x_HtrH1 = (Htr / H1)**K1
            term2_x_HtrH2 = (Htr / H2)**K2
            gamma1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval) - incomplete_gamma_upper(term1_a, term1_x_HtrH1)
            gamma2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2)
            res = N_val * (H1 * gamma1 + H2 * gamma2)
        ELSE
            res = N_val * H2 * incomplete_gamma_upper(term2_a, term1_x_ln_Nval)
        END IF
    END FUNCTION calculate_H1N

    FUNCTION perform_wave_analysis(results) RESULT(success)
        TYPE(wave_results_t), INTENT(INOUT) :: results
        LOGICAL :: success
        REAL(KIND=dp) :: sqrt_m0
        REAL(KIND=dp) :: H1_3_Hrms_capped, H1_10_Hrms_capped, H1_50_Hrms_capped, &
                         H1_100_Hrms_capped, H1_250_Hrms_capped, H1_1000_Hrms_capped

        success = .TRUE.

        ! Step 1: Calculate primary parameters
        results%m0 = (results%Hm0 / 4.0_dp)**2
        sqrt_m0 = SQRT(results%m0)
        results%Hrms = (2.69_dp + 3.24_dp * sqrt_m0 / results%d) * sqrt_m0
        results%tanAlpha = 1.0_dp / results%slopeM
        results%Htr_dim = (0.35_dp + 5.8_dp * results%tanAlpha) * results%d
        IF (results%Hrms > 0.0_dp) THEN
            results%Htr_tilde = results%Htr_dim / results%Hrms
        ELSE
            results%Htr_tilde = 0.0_dp
        END IF

        ! OVERSHOOT-PREVENTION: Method 1 - H_tr Threshold Switch
        IF (results%Htr_tilde > 2.75_dp) THEN
            results%distribution_type = "Rayleigh"
            ! Use theoretically exact H(1/N)/Hm0 ratios for a pure Rayleigh distribution.
            results%H1_3_dim    = 1.001075736951740_dp * results%Hm0
            results%H1_10_dim   = 1.272734273369137_dp * results%Hm0
            results%H1_50_dim   = 1.560113379974762_dp * results%Hm0
            results%H1_100_dim  = 1.668233372358517_dp * results%Hm0
            results%H1_250_dim  = 1.801017222497626_dp * results%Hm0
            results%H1_1000_dim = 1.984835590575388_dp * results%Hm0

            IF (results%Hrms > 0.0_dp) THEN
                results%H1_3_Hrms = results%H1_3_dim / results%Hrms
                results%H1_10_Hrms = results%H1_10_dim / results%Hrms
                results%H1_50_Hrms = results%H1_50_dim / results%Hrms
                results%H1_100_Hrms = results%H1_100_dim / results%Hrms
                results%H1_250_Hrms = results%H1_250_dim / results%Hrms
                results%H1_1000_Hrms = results%H1_1000_dim / results%Hrms
            END IF
        ELSE
            results%distribution_type = "B&G"
            ! Step 2: Solve for H1/Hrms and H2/Hrms
            IF (.NOT. newtonRaphsonSystemSolver(results%Htr_tilde, results%H1_Hrms, results%H2_Hrms, 100)) THEN
                success = .FALSE.; RETURN
            END IF

            ! Step 3: Calculate H1/N quantiles
            results%H1_3_Hrms = calculate_H1N(3.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)
            results%H1_10_Hrms = calculate_H1N(10.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)
            results%H1_50_Hrms = calculate_H1N(50.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)
            results%H1_100_Hrms = calculate_H1N(100.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)
            results%H1_250_Hrms = calculate_H1N(250.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)
            results%H1_1000_Hrms = calculate_H1N(1000.0_dp, results%H1_Hrms, results%H2_Hrms, results%Htr_tilde)

            ! Step 4: Convert to dimensional values
            results%H1_dim = results%H1_Hrms * results%Hrms
            results%H2_dim = results%H2_Hrms * results%Hrms
            results%H1_3_dim = results%H1_3_Hrms * results%Hrms
            results%H1_10_dim = results%H1_10_Hrms * results%Hrms
            results%H1_50_dim = results%H1_50_Hrms * results%Hrms
            results%H1_100_dim = results%H1_100_Hrms * results%Hrms
            results%H1_250_dim = results%H1_250_Hrms * results%Hrms
            results%H1_1000_dim = results%H1_1000_Hrms * results%Hrms

            ! OVERSHOOT-PREVENTION: Method 2 - Capping
            ! Use theoretically exact H(1/N)/Hm0 ratios for a pure Rayleigh distribution.
            results%H1_3_dim    = MIN(results%H1_3_dim,    1.001075736951740_dp * results%Hm0)
            results%H1_10_dim   = MIN(results%H1_10_dim,   1.272734273369137_dp * results%Hm0)
            results%H1_50_dim   = MIN(results%H1_50_dim,   1.560113379974762_dp * results%Hm0)
            results%H1_100_dim  = MIN(results%H1_100_dim,  1.668233372358517_dp * results%Hm0)
            results%H1_250_dim  = MIN(results%H1_250_dim,  1.801017222497626_dp * results%Hm0)
            results%H1_1000_dim = MIN(results%H1_1000_dim, 1.984835590575388_dp * results%Hm0)
        END IF

        ! Step 5: Calculate diagnostic ratios (always do this last)
        IF (results%H1_3_dim > 0.0_dp) THEN
            H1_3_Hrms_capped = results%H1_3_dim / results%Hrms
            H1_10_Hrms_capped = results%H1_10_dim / results%Hrms
            H1_50_Hrms_capped = results%H1_50_dim / results%Hrms
            H1_100_Hrms_capped = results%H1_100_dim / results%Hrms
            H1_250_Hrms_capped = results%H1_250_dim / results%Hrms
            H1_1000_Hrms_capped = results%H1_1000_dim / results%Hrms

            results%ratio_1_10_div_1_3 = H1_10_Hrms_capped / H1_3_Hrms_capped
            results%ratio_1_50_div_1_3 = H1_50_Hrms_capped / H1_3_Hrms_capped
            results%ratio_1_100_div_1_3 = H1_100_Hrms_capped / H1_3_Hrms_capped
            results%ratio_1_250_div_1_3 = H1_250_Hrms_capped / H1_3_Hrms_capped
            results%ratio_1_1000_div_1_3 = H1_1000_Hrms_capped / H1_3_Hrms_capped
        ELSE
            results%ratio_1_10_div_1_3 = 0.0_dp
            results%ratio_1_50_div_1_3 = 0.0_dp
            results%ratio_1_100_div_1_3 = 0.0_dp
            results%ratio_1_250_div_1_3 = 0.0_dp
            results%ratio_1_1000_div_1_3 = 0.0_dp
        END IF
    END FUNCTION perform_wave_analysis

    ! =======================================================================
    ! I/O and Report Generation
    ! =======================================================================
    SUBROUTINE get_user_input(Hm0_out, d_out, slopeM_out)
        REAL(KIND=dp), INTENT(OUT) :: Hm0_out, d_out, slopeM_out
        CHARACTER(LEN=256) :: arg_str
        INTEGER :: arg_count, iostat_val

        arg_count = COMMAND_ARGUMENT_COUNT()
        IF (arg_count >= 3) THEN
            CALL GET_COMMAND_ARGUMENT(1, arg_str)
            READ(arg_str, *, IOSTAT=iostat_val) Hm0_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Invalid Hm0 argument."
                ERROR STOP 2
            END IF
            CALL GET_COMMAND_ARGUMENT(2, arg_str)
            READ(arg_str, *, IOSTAT=iostat_val) d_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Invalid d argument."
                ERROR STOP 3
            END IF
            CALL GET_COMMAND_ARGUMENT(3, arg_str)
            READ(arg_str, *, IOSTAT=iostat_val) slopeM_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Invalid slopeM argument."
                ERROR STOP 4
            END IF
        ELSE
            WRITE(*, '(A)', ADVANCE='NO') "Enter Hm0 (m): "
            READ(*,*,IOSTAT=iostat_val) Hm0_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Input error for Hm0."
                ERROR STOP 5
            END IF
            WRITE(*, '(A)', ADVANCE='NO') "Enter water depth d (m): "
            READ(*,*,IOSTAT=iostat_val) d_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Input error for d."
                ERROR STOP 6
            END IF
            WRITE(*, '(A)', ADVANCE='NO') "Enter beach slope (m): "
            READ(*,*,IOSTAT=iostat_val) slopeM_out
            IF (iostat_val /= 0) THEN
                WRITE(*,*) "Input error for slopeM."
                ERROR STOP 7
            END IF
        END IF

        IF (Hm0_out <= 0.0_dp .OR. d_out <= 0.0_dp .OR. slopeM_out <= 0.0_dp) THEN
            WRITE(*,*) "ERROR: All inputs must be positive."
            ERROR STOP 8
        END IF
    END SUBROUTINE get_user_input

    SUBROUTINE generate_and_write_report(results)
        !---------------------------------------------------------------------
        ! Purpose:
        !   Builds a detailed report from the results object and writes it
        !   to "report.txt" and standard output.
        !---------------------------------------------------------------------
        TYPE(wave_results_t), INTENT(IN) :: results
        CHARACTER(LEN=4000) :: report_str ! A large enough string to hold the report
        CHARACTER(LEN=120) :: temp_line ! Temporary string for building each line
        CHARACTER(LEN=20)  :: num_str, slopeM_str, tanAlpha_str ! Temporary string for number formatting
        INTEGER :: iostat_val

        report_str = ''

        ! Build the report string with formatted output using internal write
        WRITE(temp_line, '(A)') "======================"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   INPUT PARAMETERS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "======================"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%Hm0
        WRITE(temp_line, '(A, A)') "Hm0 (m)         : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%d
        WRITE(temp_line, '(A, A)') "d (m)           : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(slopeM_str, '(F20.4)') results%slopeM
        WRITE(tanAlpha_str, '(F20.4)') results%tanAlpha
        temp_line = "Beach slope (m) : " // TRIM(ADJUSTL(slopeM_str)) &
                  // "   (tan(alpha) = " // TRIM(ADJUSTL(tanAlpha_str)) // ")"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A') ! Blank line

        WRITE(temp_line, '(A, A)') "Distribution Used : ", TRIM(results%distribution_type)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A') ! Blank line

        ! Section: CALCULATED PARAMETERS
        WRITE(temp_line, '(A)') "==========================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   CALCULATED PARAMETERS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "==========================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%m0
        WRITE(temp_line, '(A, A)') "Free-surface variance m0 (m^2)   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%Hrms
        WRITE(temp_line, '(A, A)') "Mean square wave height Hrms (m) : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%Htr_dim
        WRITE(temp_line, '(A, A)') "Transitional wave height Htr (m) : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%Htr_tilde
        WRITE(temp_line, '(A, A)') "Dimensionless H~_tr (Htr/Hrms)   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        report_str = TRIM(report_str) // NEW_LINE('A') ! Blank line

        ! Section: DIMENSIONLESS WAVE HEIGHTS (H/Hrms)
        WRITE(temp_line, '(A)') "========================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIMENSIONLESS WAVE HEIGHTS (H/Hrms)"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "========================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_Hrms
        WRITE(temp_line, '(A, A)') "H1/Hrms       : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H2_Hrms
        WRITE(temp_line, '(A, A)') "H2/Hrms       : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_3_Hrms
        WRITE(temp_line, '(A, A)') "H1/3 / Hrms   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_10_Hrms
        WRITE(temp_line, '(A, A)') "H1/10 / Hrms  : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_50_Hrms
        WRITE(temp_line, '(A, A)') "H1/50 / Hrms  : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_100_Hrms
        WRITE(temp_line, '(A, A)') "H1/100 / Hrms : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_250_Hrms
        WRITE(temp_line, '(A, A)') "H1/250 / Hrms : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_1000_Hrms
        WRITE(temp_line, '(A, A)') "H1/1000 /Hrms : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A') ! Blank line

        ! Section: DIMENSIONAL WAVE HEIGHTS (m)
        WRITE(temp_line, '(A)') "=================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIMENSIONAL WAVE HEIGHTS (m)"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "=================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_dim
        WRITE(temp_line, '(A, A)') "H1 (m)        : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H2_dim
        WRITE(temp_line, '(A, A)') "H2 (m)        : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_3_dim
        WRITE(temp_line, '(A, A)') "H1/3 (m)      : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_10_dim
        WRITE(temp_line, '(A, A)') "H1/10 (m)     : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_50_dim
        WRITE(temp_line, '(A, A)') "H1/50 (m)     : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_100_dim
        WRITE(temp_line, '(A, A)') "H1/100 (m)    : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_250_dim
        WRITE(temp_line, '(A, A)') "H1/250 (m)    : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%H1_1000_dim
        WRITE(temp_line, '(A, A)') "H1/1000 (m)   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A') ! Blank line

        ! Section: DIAGNOSTIC RATIOS
        WRITE(temp_line, '(A)') "======================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIAGNOSTIC RATIOS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "======================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%ratio_1_10_div_1_3
        WRITE(temp_line, '(A, A)') "(H1/10)/(H1/3)   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%ratio_1_50_div_1_3
        WRITE(temp_line, '(A, A)') "(H1/50)/(H1/3)   : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%ratio_1_100_div_1_3
        WRITE(temp_line, '(A, A)') "(H1/100)/(H1/3)  : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%ratio_1_250_div_1_3
        WRITE(temp_line, '(A, A)') "(H1/250)/(H1/3)  : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(num_str, '(F20.4)') results%ratio_1_1000_div_1_3
        WRITE(temp_line, '(A, A)') "(H1/1000)/(H1/3) : ", TRIM(ADJUSTL(num_str))
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A') // NEW_LINE('A')

        WRITE(temp_line, '(A)') "End of Report"
        report_str = TRIM(report_str) // TRIM(temp_line)

        ! Write to file
        OPEN(UNIT=REPORT_UNIT, FILE="report.txt", STATUS="REPLACE", ACTION="WRITE", IOSTAT=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Error: Could not open report.txt for writing."
        ELSE
            WRITE(REPORT_UNIT, '(A)') TRIM(report_str)
            CLOSE(REPORT_UNIT)
        END IF

        ! Write to standard output
        WRITE(*, '(A)') TRIM(report_str)

    END SUBROUTINE generate_and_write_report

END PROGRAM ShallowWaterWaves
