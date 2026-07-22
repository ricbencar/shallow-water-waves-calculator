!***********************************************************************
! Program: carvalho2025_table.f90
!
! Detailed Description:
! This program is designed to compute and tabulate normalized wave height parameters
! for shallow-water environments, utilizing the Composite Weibull distribution model.
! All wave heights are normalized by the root-mean-square wave height (Hrms) to provide
! dimensionless ratios, which are fundamental for understanding wave statistics in complex
! shallow foreshore conditions where wave breaking significantly alters traditional
! wave height distributions (e.g., Rayleigh distribution).
!
! The core functionality of this program involves a robust numerical solution using
! the Newton-Raphson method for a system of two non-linear equations:
!
! 1.  **Iterating through Normalized Transitional Wave Height (Htr/Hrms)**:
! The program systematically varies the normalized transitional wave height (Htr_Hrms)
! from 0.01 to 3.50 in increments of 0.01. The transitional wave height (Htr) is a
! critical parameter in the Composite Weibull distribution, marking the point where
! the wave height distribution transitions from a Rayleigh-like behavior (for smaller waves)
! to a different, more complex behavior influenced by wave breaking (for larger waves).
!
! 2.  **Solving for Normalized Scale Parameters (H1/Hrms and H2/Hrms) Simultaneously**:
! For each Htr_Hrms value, the program solves a system of two non-linear equations
! to determine H1_Hrms and H2_Hrms simultaneously. This is a significant improvement
! over sequential solving, ensuring consistency and accuracy in the relationship
! between the two Weibull components.
! The two equations are:
! a) The normalized Hrms equation, which states that the overall normalized Hrms of the Composite Weibull distribution
! must precisely equal one.
! b) The continuity condition, which ensures
! a smooth transition between the two Weibull components at Htr.
! The Newton-Raphson method for systems of equations is employed, which involves
! calculating the Jacobian matrix (matrix of partial derivatives) at each iteration.
!
! 3.  **Computing Normalized Quantile Wave Heights (H(1/N)/Hrms)**:
! Once H1_Hrms and H2_Hrms are successfully determined, the program computes the
! corresponding normalized quantile wave heights for a predefined set of exceedance
! probabilities (represented by N values, e.g., N=3 for H1/3, N=10 for H1/10, etc.).
! These calculations utilize specific formulas,
! depending on the relationship between Htr_Hrms and the specific quantile being calculated.
! These quantiles represent the mean of the highest 1/N-part of the wave heights,
! providing crucial insights into the statistical properties of extreme waves.
!
! Output:
! The program generates a formatted text file named "carvalho2025_table.txt". This file
! contains a comprehensive table of Htr_Hrms, H1_Hrms, H2_Hrms, and the various
! H(1/N)/Hrms values, providing a valuable reference for wave height analysis in
! shallow-water coastal engineering applications.
!
! All gamma functions applied, either complete or incomplete, are non-normalized.
!
! Compilation Instructions:
! To compile this program with high optimization levels, you can use a gfortran command similar to this:
!
! gfortran -O3 -march=native -std=f2008 -Wall -Wextra -pedantic \
! -fno-underscoring -o carvalho2025_table.exe carvalho2025_table.f90
!
! Note: Fortran compilers often link to system math libraries automatically.
! For special functions like incomplete gamma, it might be necessary to use a library,
! or implement them. This code provides an implementation.
!**********************************************************************

PROGRAM carvalho2025_table
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER :: k1 = 2.0_8    ! Exponent for the first part of the Composite Weibull distribution.
                                           ! As per the theoretical foundation (Groenendijk, 1998, Section 2.1 and 3.3.2),
                                           ! the initial part of the wave height distribution in shallow water is assumed
                                           ! to be Rayleigh-shaped. A Rayleigh distribution is a special case of the Weibull
                                           ! distribution with an exponent (shape parameter) of 2.0.
    REAL(KIND=8), PARAMETER :: k2 = 3.6_8    ! Exponent for the second part of the Composite Weibull distribution.
                                           ! This value was empirically determined through calibration and optimization
                                           ! processes described in Groenendijk (1998, Section 2.1) and Groenendijk & Van Gent (1998).
                                           ! It reflects the altered shape of the wave height distribution for larger waves
                                           ! due to depth-induced breaking.

    REAL(KIND=8), PARAMETER :: EPSILON = 1.0E-12_8 ! A small value (10^-12) indicating the maximum allowable error or difference
                                               ! between successive iterations for a solution to be considered converged.
    REAL(KIND=8), PARAMETER :: JACOBIAN_DX = 1.0E-8_8 ! Small step size for finite difference approximation of Jacobian derivatives.

    ! No INTERFACE block needed here because all procedures are internal procedures
    ! defined within the CONTAINS section of the main program.

    INTEGER, PARAMETER :: OUTPUT_UNIT = 10
    REAL(KIND=8) :: Htr_Hrms_start, Htr_Hrms_end, Htr_Hrms_step
    REAL(KIND=8) :: Htr_Hrms, H1_normalized, H2_normalized
    INTEGER :: N_values(6)
    INTEGER :: i, N_val_idx
    REAL(KIND=8) :: h1n_values(6)
    LOGICAL :: converged

    ! Define the range and step size for the normalized transitional wave height (Htr/Hrms).
    Htr_Hrms_start = 0.01_8
    Htr_Hrms_end = 3.50_8
    Htr_Hrms_step = 0.01_8

    ! Define an array of integer N values for which the H(1/N)/Hrms quantiles will be calculated.
    N_values = (/3, 10, 50, 100, 250, 1000/)

    ! Open an output file named "carvalho2025_table.txt" for writing.
    OPEN(UNIT=OUTPUT_UNIT, FILE="carvalho2025_table.txt", STATUS="REPLACE", ACTION="WRITE", IOSTAT=i)
    IF (i /= 0) THEN
        WRITE(*,*) "Error opening output file: carvalho2025_table.txt"
        STOP 1
    END IF

    ! Write the header row to the output file.
    WRITE(OUTPUT_UNIT, '(A9, A9, A9, A11, A11, A11, A12, A12, A13)') &
         "Htr/Hrms", "H1/Hrms", "H2/Hrms", "H1/3/Hrms", "H1/10/Hrms", &
         "H1/50/Hrms", "H1/100/Hrms", "H1/250/Hrms", "H1/1000/Hrms"

    ! Main loop: Iterate through the defined range of Htr_Hrms values.
    Htr_Hrms = Htr_Hrms_start
    DO WHILE (Htr_Hrms <= Htr_Hrms_end + EPSILON)
        ! Initialize h1n_values to avoid carrying over values from previous iterations if an error occurs
        h1n_values = -999.0_8 ! Sentinel value for error indication

        converged = .FALSE.
        H1_normalized = -1.0_8
        H2_normalized = -1.0_8

        ! Call the newtonRaphsonSystemSolver function to find H1_normalized and H2_normalized simultaneously.
        converged = newtonRaphsonSystemSolver(Htr_Hrms, H1_normalized, H2_normalized, EPSILON, 100)

        IF (.NOT. converged) THEN
            WRITE(*,*) "Calculation error for Htr_Hrms =", Htr_Hrms, &
                 ": Newton-Raphson system solver failed to find H1_Hrms and H2_Hrms."
            WRITE(OUTPUT_UNIT, '(F9.5, 8(A9))') Htr_Hrms, "ERROR", "ERROR", "ERROR", &
                 "ERROR", "ERROR", "ERROR", "ERROR", "ERROR"
        ELSE IF (H1_normalized <= 0.0_8 .OR. H2_normalized <= 0.0_8) THEN
            ! This case should ideally be handled by the solver, but as a safeguard.
            WRITE(*,*) "Calculation error for Htr_Hrms =", Htr_Hrms, &
                 ": H1_normalized or H2_normalized became non-positive."
            WRITE(OUTPUT_UNIT, '(F9.5, 8(A9))') Htr_Hrms, "ERROR", "ERROR", "ERROR", &
                 "ERROR", "ERROR", "ERROR", "ERROR", "ERROR"
        ELSE
            ! Calculate the H1/N values
            DO N_val_idx = 1, SIZE(N_values)
                h1n_values(N_val_idx) = calculate_H1N(REAL(N_values(N_val_idx), KIND=8), H1_normalized, &
                                                      H2_normalized, k1, k2, Htr_Hrms)
            END DO

            ! Write the primary normalized parameters (Htr/Hrms, H1/Hrms, H2/Hrms) to the file,
            ! formatted with specified widths.
            WRITE(OUTPUT_UNIT, '(F9.5, F9.5, F9.5, F11.5, F11.5, F11.5, F12.5, F12.5, F13.5)') &
                 Htr_Hrms, H1_normalized, H2_normalized, &
                 h1n_values(1), h1n_values(2), h1n_values(3), &
                 h1n_values(4), h1n_values(5), h1n_values(6)
        END IF

        Htr_Hrms = Htr_Hrms + Htr_Hrms_step
    END DO

    CLOSE(OUTPUT_UNIT)
    WRITE(*,*) "Calculation complete. Results saved to carvalho2025_table.txt"

CONTAINS

    ! All gamma functions applied, either complete or incomplete, are non-normalized.

    REAL(KIND=8) FUNCTION incomplete_gamma_lower(a, x)
        ! @brief Computes the unnormalized lower incomplete gamma function γ(a, x) = ∫[0,x] t^(a-1)e^(-t) dt.
        !
        ! This function is a critical component for calculating moments of Weibull distributions,
        ! which are fundamental to the Composite Weibull model. It employs a hybrid approach
        ! for numerical stability and accuracy:
        ! - For small values of 'x' (specifically, x < a + 1.0), it uses a series expansion.
        ! - For larger values of 'x', it utilizes a continued fraction expansion.
        ! This adaptive strategy ensures robust and precise computation across different input ranges.
        !
        ! @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
        ! @param x The 'x' parameter (upper integration limit) of the incomplete gamma function. Must be non-negative.
        ! @return The computed value of γ(a, x). Returns `Huge(1.0_8)` (representing an, error) if input parameters are invalid
        ! or if the series/continued fraction fails to converge within the maximum iterations.
        REAL(KIND=8), INTENT(IN) :: a, x
        INTEGER, PARAMETER :: MAXIT = 500      ! Maximum number of iterations allowed for either the series expansion or the continued fraction.
        REAL(KIND=8), PARAMETER :: LOCAL_EPS = 1.0E-16_8 ! A very small tolerance value (10^-16) used for checking convergence
                                                       ! within the gamma function calculations, ensuring high numerical precision.
        REAL(KIND=8) :: gln, ap, sum_val, del_val, b, c, d, h, an, del
        INTEGER :: n_iter, i

        ! Input validation: Essential for preventing mathematical errors and ensuring function robustness.
        IF (a <= 0.0_8 .OR. x < 0.0_8) THEN
            incomplete_gamma_lower = HUGE(1.0_8) ! Use HUGE for NaN equivalent
            RETURN
        END IF
        ! Corrected: Check if x is very close to 0.0_8 instead of exactly equal.
        IF (ABS(x - 0.0_8) < TINY(1.0_8)) THEN
            incomplete_gamma_lower = 0.0_8
            RETURN
        END IF

        ! Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
        ! `LOG_GAMMA` is used for improved numerical stability.
        gln = LOG_GAMMA(a)

        ! Conditional logic to select the appropriate numerical method (series or continued fraction).
        IF (x < a + 1.0_8) THEN  ! Use series expansion for x < a+1.
            ap = a
            sum_val = 1.0_8 / a
            del_val = sum_val
            DO n_iter = 1, MAXIT
                ap = ap + 1.0_8
                del_val = del_val * x / ap
                sum_val = sum_val + del_val
                IF (ABS(del_val) < ABS(sum_val) * LOCAL_EPS) THEN
                    incomplete_gamma_lower = sum_val * EXP(-x + a * LOG(x) - gln) * GAMMA(a)
                    RETURN
                END IF
            END DO
            incomplete_gamma_lower = HUGE(1.0_8) ! Convergence failed
        ELSE  ! Use continued fraction for x >= a+1.
            b = x + 1.0_8 - a
            c = 1.0_8 / TINY(1.0_8) ! Equivalent to C++ numeric_limits<long double>::min()
            d = 1.0_8 / b
            h = d
            DO i = 1, MAXIT
                an = -1.0_8 * REAL(i, KIND=8) * (REAL(i, KIND=8) - a)
                b = b + 2.0_8
                d = an * d + b
                IF (ABS(d) < TINY(1.0_8)) d = TINY(1.0_8)
                c = b + an / c
                IF (ABS(c) < TINY(1.0_8)) c = TINY(1.0_8)
                d = 1.0_8 / d
                del = d * c
                h = h * del
                IF (ABS(del - 1.0_8) < LOCAL_EPS) EXIT
            END DO
            incomplete_gamma_lower = (1.0_8 - EXP(-x + a * LOG(x) - gln) * h) * GAMMA(a)
            IF (incomplete_gamma_lower < 0.0_8) incomplete_gamma_lower = 0.0_8 ! Clamp to 0.0
        END IF
    END FUNCTION incomplete_gamma_lower


    REAL(KIND=8) FUNCTION incomplete_gamma_upper(a, x)
        ! @brief Computes the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
        !
        ! This function calculates the unnormalized form of the upper incomplete gamma function.
        ! It is derived from the complete gamma function Γ(a) and the unnormalized lower incomplete gamma function γ(a, x),
        ! using the identity Γ(a,x) = Γ(a) - γ(a,x).
        !
        ! @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
        ! @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
        ! @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
        ! Returns `Huge(1.0_8)` if input parameters are invalid.
        REAL(KIND=8), INTENT(IN) :: a, x

        IF (a <= 0.0_8) THEN
            WRITE(*,*) "Error: incomplete_gamma_upper: 'a' must be positive."
            incomplete_gamma_upper = HUGE(1.0_8)
            RETURN
        END IF
        IF (x < 0.0_8) THEN
            WRITE(*,*) "Error: incomplete_gamma_upper: 'x' must be non-negative."
            incomplete_gamma_upper = HUGE(1.0_8)
            RETURN
        END IF
        incomplete_gamma_upper = GAMMA(a) - incomplete_gamma_lower(a, x)
    END FUNCTION incomplete_gamma_upper

    REAL(KIND=8) FUNCTION calculate_HN(N, H1, H2, k1_val, k2_val, Htr)
        ! @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
        !
        ! This function determines the specific wave height (H) such that the probability of a wave
        ! exceeding this height is 1/N. This is a key statistical measure. The calculation depends
        ! on whether the target wave height falls into the first or second part of the Composite Weibull
        ! distribution, which is separated by the transitional wave height (Htr).
        ! The logic directly implements the formulas from Appendix A.2.1 of Groenendijk (1998).
        !
        ! @param N The N value (e.g., 3 for H1/3, 100 for H1%). Must be strictly greater than 1.0
        ! because `log(N)` is used, and `N=1` would result in `log(1)=0`.
        ! @param H1 The scale parameter of the first Weibull distribution. Must be positive.
        ! @param H2 The scale parameter of the second Weibull distribution. Must be positive.
        ! @param k1_val The exponent (shape parameter) of the first Weibull distribution. Must be positive.
        ! @param k2_val The exponent (shape parameter) of the second Weibull distribution. Must be positive.
        ! @param Htr The transitional wave height, which defines the boundary between the two Weibull parts.
        ! @return The calculated value of HN.
        ! Returns `Huge(1.0_8)` if input parameters are invalid.
        REAL(KIND=8), INTENT(IN) :: N, H1, H2, k1_val, k2_val, Htr
        REAL(KIND=8) :: HN_candidate1

        IF (N <= 1.0_8) THEN
            WRITE(*,*) "Error: calculate_HN: N must be greater than 1 for LOG(N)."
            calculate_HN = HUGE(1.0_8)
            RETURN
        END IF
        IF (H1 <= 0.0_8 .OR. H2 <= 0.0_8) THEN
            WRITE(*,*) "Error: calculate_HN: H1 and H2 must be positive."
            calculate_HN = HUGE(1.0_8)
            RETURN
        END IF
        IF (k1_val <= 0.0_8 .OR. k2_val <= 0.0_8) THEN
            WRITE(*,*) "Error: calculate_HN: k1 and k2 must be positive."
            calculate_HN = HUGE(1.0_8)
            RETURN
        END IF

        HN_candidate1 = H1 * (LOG(N))**(1.0_8 / k1_val)

        IF (HN_candidate1 < Htr - EPSILON) THEN
            calculate_HN = HN_candidate1
        ELSE
            calculate_HN = H2 * (LOG(N))**(1.0_8 / k2_val)
        END IF
    END FUNCTION calculate_HN


    REAL(KIND=8) FUNCTION calculate_H1N(N_val, H1, H2, k1_val, k2_val, Htr)
        ! @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N) for the Composite Weibull distribution.
        !
        ! This function computes a characteristic wave height that represents the average height of the
        ! highest N-th fraction of waves in a given wave field. This is a commonly used metric in
        ! oceanography and coastal engineering (e.g., H1/3 for significant wave height).
        ! The calculation involves integrals of the probability density function and depends on
        ! whether the relevant wave heights fall into the first or second part of the Composite Weibull distribution.
        ! The implementation follows the detailed derivations in Appendix A.2.2 of Groenendijk (1998).
        !
        ! @param N_val The N parameter for H1/N (e.g., 3 for H1/3, 10 for H1/10). Must be strictly greater than 1.
        ! @param H1 The scale parameter of the first Weibull distribution. Must be positive.
        ! @param H2 The scale parameter of the second Weibull distribution. Must be positive.
        ! @param k1_val The exponent (shape parameter) of the first Weibull distribution. Must be positive.
        ! @param k2_val The exponent (shape parameter) of the second Weibull distribution. Must be positive.
        ! @param Htr The transitional wave height.
        ! @return The calculated value of H1/N.
        ! Returns `Huge(1.0_8)` if input parameters are invalid.
        REAL(KIND=8), INTENT(IN) :: N_val, H1, H2, k1_val, k2_val, Htr
        REAL(KIND=8) :: H_N_val, term1_a, term1_x_ln_Nval, term1_x_HtrH1, &
                         gamma_term1_part1, gamma_term1_part2, gamma_term1, &
                         term2_a, term2_x_HtrH2, gamma_term2, &
                         term_a, term_x

        IF (H1 <= 0.0_8 .OR. H2 <= 0.0_8) THEN
            WRITE(*,*) "Error: calculate_H1N: H1 and H2 must be positive."
            calculate_H1N = HUGE(1.0_8)
            RETURN
        END IF
        IF (k1_val <= 0.0_8 .OR. k2_val <= 0.0_8) THEN
            WRITE(*,*) "Error: calculate_H1N: k1 and k2 must be positive."
            calculate_H1N = HUGE(1.0_8)
            RETURN
        END IF
        IF (N_val <= 1.0_8) THEN
            WRITE(*,*) "Error: calculate_H1N: N_val must be greater than 1."
            calculate_H1N = HUGE(1.0_8)
            RETURN
        END IF

        H_N_val = calculate_HN(N_val, H1, H2, k1_val, k2_val, Htr)

        IF (H_N_val < Htr - EPSILON) THEN
            term1_a = 1.0_8 / k1_val + 1.0_8
            term1_x_ln_Nval = LOG(N_val)
            term1_x_HtrH1 = (Htr / H1)**k1_val

            gamma_term1_part1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval)
            gamma_term1_part2 = incomplete_gamma_upper(term1_a, term1_x_HtrH1)
            gamma_term1 = gamma_term1_part1 - gamma_term1_part2

            term2_a = 1.0_8 / k2_val + 1.0_8
            term2_x_HtrH2 = (Htr / H2)**k2_val

            gamma_term2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2)

            calculate_H1N = N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2
        ELSE
            term_a = 1.0_8 / k2_val + 1.0_8
            term_x = LOG(N_val)
            calculate_H1N = N_val * H2 * incomplete_gamma_upper(term_a, term_x)
        END IF
    END FUNCTION calculate_H1N


    REAL(KIND=8) FUNCTION F1(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val)
        ! @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
        !
        ! This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
        ! It is derived from Equation 7.11 from Groenendijk (1998) and Equation 9 in Caires & Van Gent (2012),
        ! and states that the overall normalized Hrms of the Composite Weibull distribution must precisely equal 1.
        ! The equation directly uses the unnormalized lower incomplete gamma function, γ(a,x),
        ! and the unnormalized upper incomplete gamma function, Γ(a,x).
        !
        ! @param H1_Hrms_val The normalized scale parameter of the first Weibull distribution.
        ! @param H2_Hrms_val The normalized scale parameter of the second Weibull distribution.
        ! @param Htr_Hrms_val The normalized transitional wave height (constant for a given solve).
        ! @return The value of the first function, which should be driven to zero.
        REAL(KIND=8), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val
        REAL(KIND=8) :: arg1, arg2, term1, term2, sum_terms

        IF (H1_Hrms_val <= 0.0_8 .OR. H2_Hrms_val <= 0.0_8) THEN
            F1 = HUGE(1.0_8)
            RETURN
        END IF

        arg1 = (Htr_Hrms_val / H1_Hrms_val)**k1
        arg2 = (Htr_Hrms_val / H2_Hrms_val)**k2

        term1 = H1_Hrms_val * H1_Hrms_val * incomplete_gamma_lower(2.0_8 / k1 + 1.0_8, arg1)
        term2 = H2_Hrms_val * H2_Hrms_val * incomplete_gamma_upper(2.0_8 / k2 + 1.0_8, arg2)

        sum_terms = term1 + term2
        IF (sum_terms < 0.0_8) sum_terms = 0.0_8

        F1 = SQRT(sum_terms) - 1.0_8
    END FUNCTION F1


    REAL(KIND=8) FUNCTION F2(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val)
        ! @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
        !
        ! This equation represents the continuity condition between the two Weibull distributions
        ! at the transitional wave height Htr. It is derived from Equation 3.4 in Groenendijk (1998)
        ! and Equation 8 in Caires & Van Gent (2012):
        ! `(Htr/H1)^k1 = (Htr/H2)^k2`. Rearranging this gives `(Htr/H1)^k1 - (Htr/H2)^k2 = 0`.
        !
        ! @param H1_Hrms_val The normalized scale parameter of the first Weibull distribution.
        ! @param H2_Hrms_val The normalized scale parameter of the second Weibull distribution.
        ! @param Htr_Hrms_val The normalized transitional wave height (constant for a given solve).
        ! @return The value of the second function, which should be driven to zero.
        REAL(KIND=8), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val

        IF (H1_Hrms_val <= 0.0_8 .OR. H2_Hrms_val <= 0.0_8) THEN
            F2 = HUGE(1.0_8)
            RETURN
        END IF
        F2 = (Htr_Hrms_val / H1_Hrms_val)**k1 - (Htr_Hrms_val / H2_Hrms_val)**k2
    END FUNCTION F2


    SUBROUTINE solve_linear_system_2x2(J11, J12, J21, J22, b1, b2, dx1, dx2)
        ! @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
        !
        ! This function is a helper for the Newton-Raphson method for systems.
        ! It takes the Jacobian matrix elements (J11, J12, J21, J22) and the negative
        ! function values (-F1, -F2) as the right-hand side, and computes the updates
        ! (dx1, dx2) for H1_Hrms and H2_Hrms.
        !
        ! @param J11 Element (1,1) of the Jacobian matrix (dF1/dH1).
        ! @param J12 Element (1,2) of the Jacobian matrix (dF1/dH2).
        ! @param J21 Element (2,1) of the Jacobian matrix (dF2/dH1).
        ! @param J22 Element (2,2) of the Jacobian matrix (dF2/dH2).
        ! @param b1 Right-hand side for the first equation (-F1).
        ! @param b2 Right-hand side for the second equation (-F2).
        ! @param dx1 Output: The calculated change for H1_Hrms.
        ! @param dx2 Output: The calculated change for H2_Hrms.
        REAL(KIND=8), INTENT(IN) :: J11, J12, J21, J22, b1, b2
        REAL(KIND=8), INTENT(OUT) :: dx1, dx2
        REAL(KIND=8) :: determinant

        determinant = J11 * J22 - J12 * J21

        IF (ABS(determinant) < EPSILON * 100.0_8) THEN ! Using EPSILON for comparison with tiny
            dx1 = 0.0_8
            dx2 = 0.0_8
            RETURN
        END IF

        dx1 = (b1 * J22 - b2 * J12) / determinant
        dx2 = (J11 * b2 - J21 * b1) / determinant
    END SUBROUTINE solve_linear_system_2x2


    SUBROUTINE get_initial_guesses(Htr_Hrms_val, H1_initial, H2_initial)
        ! @brief Provides initial guesses for H1_Hrms and H2_Hrms based on Htr_Hrms.
        !
        ! A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
        ! This function uses an empirical regression for H1_Hrms, and then derives H2_Hrms from the
        ! continuity condition, assuming an initial relationship.
        ! The regression for H1_Hrms is based on Groenendijk's findings.
        !
        ! @param Htr_Hrms_val The normalized transitional wave height.
        ! @param H1_initial Output: The initial guess for H1_Hrms.
        ! @param H2_initial Output: The initial guess for H2_Hrms.
        REAL(KIND=8), INTENT(IN) :: Htr_Hrms_val
        REAL(KIND=8), INTENT(OUT) :: H1_initial, H2_initial
        integer, parameter :: dp = kind(1.0d0)

        H1_initial = 0.9718670705250743_dp + 1.115952604282648_dp * &
                     Htr_Hrms_val**(-0.7970446117540275_dp) * EXP(-1.449005086812895_dp * Htr_Hrms_val)
        H2_initial = 1.059259665431797_dp + (0.2059286860468916_dp * Htr_Hrms_val) / &
                     (1.0_dp + 3.865701948059343_dp * Htr_Hrms_val**(-3.479682433107255_dp))

        IF (H1_initial <= 0.0_8) H1_initial = TINY(1.0_8)
        IF (H2_initial <= 0.0_8) H2_initial = TINY(1.0_8)
    END SUBROUTINE get_initial_guesses


    LOGICAL FUNCTION newtonRaphsonSystemSolver(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out, tol, maxit)
        ! @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
        !
        ! This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
        ! of the system of non-linear equations F1 and F2. It iteratively refines the guesses
        ! for H1_Hrms and H2_Hrms until the functions F1 and F2 are sufficiently close to zero.
        !
        ! @param Htr_Hrms_val The normalized transitional wave height (constant for this solve).
        ! @param H1_Hrms_out Output: The converged normalized scale parameter of the first Weibull distribution.
        ! @param H2_Hrms_out Output: The converged normalized scale parameter of the second Weibull distribution.
        ! @param tol The desired tolerance for convergence (maximum absolute value of F1 and F2).
        ! @param maxit The maximum number of iterations allowed.
        ! @return `.TRUE.` if the solver successfully converges, `.FALSE.` otherwise.
        REAL(KIND=8), INTENT(IN) :: Htr_Hrms_val, tol
        REAL(KIND=8), INTENT(OUT) :: H1_Hrms_out, H2_Hrms_out
        INTEGER, INTENT(IN) :: maxit
        REAL(KIND=8) :: f1_val, f2_val
        REAL(KIND=8) :: J11, J12, J21, J22
        REAL(KIND=8) :: dH1, dH2
        INTEGER :: iter

        CALL get_initial_guesses(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out)

        DO iter = 0, maxit - 1
            f1_val = F1(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)
            f2_val = F2(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)

            IF (ABS(f1_val) < tol .AND. ABS(f2_val) < tol) THEN
                newtonRaphsonSystemSolver = .TRUE.
                RETURN
            END IF

            ! Calculate the Jacobian matrix elements using central finite differences.
            J11 = (F1(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F1(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_8 * JACOBIAN_DX)
            J12 = (F1(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F1(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_8 * JACOBIAN_DX)
            J21 = (F2(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F2(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_8 * JACOBIAN_DX)
            J22 = (F2(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F2(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_8 * JACOBIAN_DX)

            CALL solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2)

            H1_Hrms_out = H1_Hrms_out + dH1
            H2_Hrms_out = H2_Hrms_out + dH2

            IF (H1_Hrms_out <= 0.0_8) H1_Hrms_out = TINY(1.0_8)
            IF (H2_Hrms_out <= 0.0_8) H2_Hrms_out = TINY(1.0_8)
        END DO

        WRITE(*,*) "Newton-Raphson system solver failed to converge for Htr_Hrms =", Htr_Hrms_val, " after ", maxit, " iterations."
        newtonRaphsonSystemSolver = .FALSE.
    END FUNCTION newtonRaphsonSystemSolver

END PROGRAM carvalho2025_table