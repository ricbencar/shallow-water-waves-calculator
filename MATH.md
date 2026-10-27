# Mathematics of Rayleigh and Composite Weibull Wave-Height Distributions

*Mathematical foundations, probability laws, moment normalization and individual-wave height statistics*

## 1. Introduction and mathematical definitions

The Rayleigh distribution and the Composite Weibull Distribution (CWD) describe statistical properties of individual wave heights within a specified sea state.

The Rayleigh distribution provides a theoretical reference for approximately linear, narrow-banded Gaussian waves. The Composite Weibull Distribution, developed by Battjes and Groenendijk (2000), represents the modification of individual-wave height statistics caused predominantly by depth-induced breaking on shallow foreshores.

The main quantities are the cumulative probability, exceedance probability, probability density, root-mean-square height, exceedance thresholds and conditional means of the highest waves.

Let $H$ denote an individual zero-crossing wave height. Its cumulative distribution function is

$$ F(h)=\Pr(H\leq h),\qquad h\geq 0. $$

The exceedance or survival function is

$$ Q(h)=\Pr(H>h)=1-F(h). $$

Where the distribution is differentiable, its probability density is

$$ f(h)=\frac{dF(h)}{dh}=-\frac{dQ(h)}{dh}. $$

Probability normalization requires

$$ \int_0^\infty f(h)\,dh=1. $$

The root-mean-square individual wave height is defined by

$$ H_{rms}=\sqrt{\operatorname E[H^2]}. $$

For a prescribed exceedance probability $1/N$, the threshold $H_N$ satisfies

$$ Q(H_N)=\frac1N. $$

The average of the highest fraction $1/N$ is

$$ H_{1/N}=\operatorname E[H\mid H>H_N]. $$

These two quantities must be distinguished. In particular,

$$ H_{1/N}=N\int_{H_N}^{\infty}h f(h)\,dh. $$

Integration by parts gives an equivalent expression:

$$ H_{1/N}=H_N+N\int_{H_N}^{\infty}Q(h)\,dh. $$

Consequently,

$$ H_{1/N}>H_N. $$

This distinction is fundamental when determining characteristic wave heights for coastal engineering applications.

## 2. Derivation of the Rayleigh distribution

### 2.1 Gaussian representation of the sea surface

Under linear random-wave theory, the surface elevation can be represented by a superposition of sinusoidal components:

$$ \eta(t)=\sum_{j=1}^{J}A_j\cos(\omega_jt+\phi_j). $$

The phases $\phi_j$ are assumed independently and uniformly distributed between zero and $2\pi$.

The expected surface elevation is zero, and the variance is

$$ m_0=\operatorname{Var}(\eta)=\frac12\sum_{j=1}^{J}A_j^2. $$

In spectral form,

$$ m_0=\int_0^\infty S_\eta(f)\,df. $$

The spectral significant wave height is defined as

$$ H_{m0}=4\sqrt{m_0}. $$

For a narrow-banded Gaussian sea, the elevation can be expressed approximately in terms of two slowly varying Gaussian quadratures:

$$ \eta(t)=X(t)\cos(\omega_pt)-Y(t)\sin(\omega_pt). $$

Under the circular Gaussian assumption,

$$ X,Y\sim\mathcal N(0,m_0). $$

The corresponding joint density is

$$ f_{X,Y}(x,y)=\frac1{2\pi m_0}\exp\left[-\frac{x^2+y^2}{2m_0}\right]. $$

Introduce the wave-envelope amplitude $A$ and phase $\phi$:

$$ X=A\cos\phi,\qquad Y=A\sin\phi. $$

The Jacobian of this transformation is

$$ \left|\frac{\partial(X,Y)}{\partial(A,\phi)}\right|=A. $$

Hence,

$$ f_{A,\phi}(a,\phi)=\frac{a}{2\pi m_0}\exp\left[-\frac{a^2}{2m_0}\right]. $$

Integration over the phase yields

$$ f_A(a)=\frac{a}{m_0}\exp\left[-\frac{a^2}{2m_0}\right]. $$

The wave-envelope amplitude therefore follows a Rayleigh distribution.

### 2.2 Transformation from amplitude to individual wave height

Under the narrow-band assumption, the individual crest-to-trough height is approximately twice the envelope amplitude:

$$ H\simeq 2A. $$

Using the transformation of probability densities,

$$ f_H(h)=\frac12 f_A(h/2). $$

Therefore,

$$ \boxed{f_R(h)=\frac{h}{4m_0}\exp\left(-\frac{h^2}{8m_0}\right).} $$

Integration gives the Rayleigh cumulative distribution:

$$ \boxed{F_R(h)=1-\exp\left(-\frac{h^2}{8m_0}\right).} $$

The exceedance probability is

$$ \boxed{Q_R(h)=\exp\left(-\frac{h^2}{8m_0}\right).} $$

This derivation requires both Gaussian quadratures and the approximation relating the wave envelope to the individual height. Gaussianity of the instantaneous surface elevation alone is insufficient to establish an exact Rayleigh distribution for broadband waves.

## 3. Mathematical properties of the Rayleigh distribution

### 3.1 Normalization by root-mean-square height

Define the Rayleigh root-mean-square height:

$$ R_R=\sqrt{8m_0}. $$

The probability distribution becomes

$$ \boxed{F_R(h)=1-\exp\left[-\left(\frac{h}{R_R}\right)^2\right].} $$

The corresponding density is

$$ f_R(h)=\frac{2h}{R_R^2}\exp\left[-\left(\frac{h}{R_R}\right)^2\right]. $$

The Rayleigh distribution therefore depends on a single dimensional scale.

### 3.2 General moments

The moment of order $r$ is

$$ \operatorname E[H^r]=\int_0^\infty h^r f_R(h)\,dh. $$

Introduce the variable

$$ u=\left(\frac{h}{R_R}\right)^2. $$

The integral becomes

$$ \operatorname E[H^r]=R_R^r\int_0^\infty u^{r/2}e^{-u}\,du. $$

Using the complete gamma function,

$$ \Gamma(s)=\int_0^\infty u^{s-1}e^{-u}\,du, $$

the general Rayleigh moment is

$$ \boxed{\operatorname E[H^r]=R_R^r\Gamma\left(1+\frac r2\right),\qquad r>-2.} $$

In particular,

$$ \operatorname E[H]=\frac{\sqrt\pi}{2}R_R. $$

$$ \operatorname E[H^2]=R_R^2. $$

$$ \operatorname{Var}(H)=R_R^2\left(1-\frac\pi4\right). $$

$$ \mathrm{CV}_R=\sqrt{\frac4\pi-1}. $$

### 3.3 Exceedance thresholds

The height exceeded with probability $1/N$ satisfies

$$ \exp\left[-\left(\frac{H_N}{R_R}\right)^2\right]=\frac1N. $$

Taking logarithms gives

$$ \boxed{H_{N,R}=R_R\sqrt{\ln N}.} $$

Expressed in terms of spectral significant wave height,

$$ \boxed{\frac{H_{N,R}}{H_{m0}}=\sqrt{\frac{\ln N}{2}}.} $$

### 3.4 Mean of the highest fraction

The mean height of the highest $1/N$ fraction is

$$ H_{1/N,R}=N\int_{H_{N,R}}^\infty h f_R(h)\,dh. $$

Substitution of the Rayleigh density gives

$$ \boxed{H_{1/N,R}=NR_R\Gamma\left(\frac32,\ln N\right).} $$

Since

$$ \Gamma\left(\frac12,x\right)=\sqrt\pi\,\operatorname{erfc}(\sqrt x), $$

one obtains

$$ \boxed{H_{1/N,R}=R_R\left[\sqrt{\ln N}+\frac{N\sqrt\pi}{2}\operatorname{erfc}(\sqrt{\ln N})\right].} $$

This formula gives the exact Rayleigh coefficient for the mean of the highest fraction, rather than the corresponding exceedance threshold.

### 3.5 Characteristic Rayleigh coefficients

For $R_R=H_{m0}/\sqrt2$:

| Fraction | $H_N/H_{m0}$ | $H_{1/N}/H_{m0}$ |
|---|---:|---:|
| $1/3$ | 0.741152 | 1.001076 |
| $1/10$ | 1.072983 | 1.272734 |
| $1/50$ | 1.398575 | 1.560113 |
| $1/100$ | 1.517427 | 1.668233 |
| $1/250$ | 1.661545 | 1.801017 |
| $1/1000$ | 1.858461 | 1.984836 |

In particular,

$$ H_{1/3,R}=4.0043029478\sqrt{m_0}. $$

This is very close to, but not exactly equal to, $H_{m0}=4\sqrt{m_0}$.

## 4. Weibull distributions and incomplete gamma functions

### 4.1 General Weibull law

The Weibull cumulative distribution has the form

$$ F_W(h)=1-\exp\left[-\left(\frac hb\right)^k\right], $$

where $b>0$ is the scale parameter and $k>0$ the shape exponent.

The corresponding density is

$$ f_W(h)=\frac{k}{b}\left(\frac hb\right)^{k-1}\exp\left[-\left(\frac hb\right)^k\right]. $$

The Rayleigh distribution is the particular case $k=2$.

For an arbitrary Weibull distribution,

$$ \operatorname E[H^r]=b^r\Gamma\left(1+\frac rk\right). $$

The root-mean-square height is therefore

$$ H_{rms}=b\sqrt{\Gamma\left(1+\frac2k\right)}. $$

### 4.2 Lower and upper incomplete gamma functions

The lower incomplete gamma function is

$$ \gamma(s,x)=\int_0^x u^{s-1}e^{-u}\,du. $$

The upper incomplete gamma function is

$$ \Gamma(s,x)=\int_x^\infty u^{s-1}e^{-u}\,du. $$

For $s>0$, they satisfy

$$ \boxed{\gamma(s,x)+\Gamma(s,x)=\Gamma(s).} $$

The derivatives with respect to the second argument are

$$ \frac{\partial\gamma(s,x)}{\partial x}=x^{s-1}e^{-x}. $$

$$ \frac{\partial\Gamma(s,x)}{\partial x}=-x^{s-1}e^{-x}. $$

These functions make it possible to integrate the lower and upper branches of the Composite Weibull Distribution analytically.

## 5. Composite Weibull Distribution

### 5.1 Mathematical formulation

The CWD represents individual wave heights using two Weibull branches separated by a transitional height $H_{tr}$.

Its survival function is

$$ \boxed{Q_C(h)=\begin{cases}\exp[-(h/H_1)^{k_1}],&h\leq H_{tr},\\ \exp[-(h/H_2)^{k_2}],&h>H_{tr}.\end{cases}} $$

The exponents established by Battjes and Groenendijk are

$$ \boxed{k_1=2.0,\qquad k_2=3.6.} $$

The lower branch has Rayleigh form. The larger exponent in the upper branch represents a more rapid decrease in exceedance probability with increasing wave height.

The density in each branch is

$$ f_i(h)=\frac{k_i}{H_i}\left(\frac h{H_i}\right)^{k_i-1}\exp\left[-\left(\frac h{H_i}\right)^{k_i}\right]. $$

### 5.2 Continuity condition

Probability continuity requires

$$ F_1(H_{tr})=F_2(H_{tr}). $$

Equivalently,

$$ \boxed{\left(\frac{H_{tr}}{H_1}\right)^{k_1}=\left(\frac{H_{tr}}{H_2}\right)^{k_2}.} $$

Define the common transition exponent $q$ by

$$ q=\left(\frac{H_{tr}}{H_1}\right)^{k_1}=\left(\frac{H_{tr}}{H_2}\right)^{k_2}. $$

Consequently,

$$ \boxed{H_i=H_{tr}q^{-1/k_i}.} $$

The probability of exceeding the transitional wave height is

$$ \boxed{Q_C(H_{tr})=e^{-q}.} $$

The reciprocal exceedance count at the transition is

$$ N_{tr}=e^q. $$

The transition is not a deterministic maximum wave height. Both branches together define an unbounded statistical distribution.

### 5.3 Density discontinuity

Immediately below the transition,

$$ f_C(H_{tr}^-)=\frac{k_1q}{H_{tr}}e^{-q}. $$

Immediately above,

$$ f_C(H_{tr}^+)=\frac{k_2q}{H_{tr}}e^{-q}. $$

Therefore,

$$ \boxed{\frac{f_C(H_{tr}^+)}{f_C(H_{tr}^-)}=\frac{k_2}{k_1}=1.8.} $$

This jump is an intrinsic feature of the empirical two-branch representation. It does not introduce a discrete probability mass at the transition.

## 6. Composite Weibull moments and RMS normalization

### 6.1 General moment

The general CWD moment is

$$ \operatorname E[H^r]=\int_0^{H_{tr}}h^rf_1(h)\,dh+\int_{H_{tr}}^\infty h^rf_2(h)\,dh. $$

Using $u=(h/H_i)^{k_i}$ in each branch gives

$$ \boxed{\operatorname E[H^r]=H_1^r\gamma\left(1+\frac r{k_1},q\right)+H_2^r\Gamma\left(1+\frac r{k_2},q\right).} $$

For the normalized height $Z=H/H_{rms}$, define

$$ \tau=\frac{H_{tr}}{H_{rms}},\qquad a_i=\frac{H_i}{H_{rms}}. $$

Then

$$ \operatorname E[Z^r]=a_1^r\gamma\left(1+\frac r{k_1},q\right)+a_2^r\Gamma\left(1+\frac r{k_2},q\right). $$

### 6.2 Second-moment constraint

Because $Z$ is normalized by its RMS,

$$ \boxed{\operatorname E[Z^2]=1.} $$

Therefore,

$$ a_1^2\gamma\left(1+\frac2{k_1},q\right)+a_2^2\Gamma\left(1+\frac2{k_2},q\right)=1. $$

Continuity gives

$$ a_1=\tau q^{-1/k_1},\qquad a_2=\tau q^{-1/k_2}. $$

Substitution reduces the normalization to

$$ \boxed{\tau^2U(q)=1,} $$

where

$$ U(q)=q^{-2/k_1}\gamma\left(1+\frac2{k_1},q\right)+q^{-2/k_2}\Gamma\left(1+\frac2{k_2},q\right). $$

For $k_1=2$ and $k_2=3.6$,

$$ \boxed{U(q)=q^{-1}\gamma(2,q)+q^{-5/9}\Gamma(14/9,q).} $$

Since

$$ \gamma(2,q)=1-(1+q)e^{-q}, $$

the complete RMS closure is

$$ \boxed{\tau^2\left[\frac{1-(1+q)e^{-q}}q+q^{-5/9}\Gamma(14/9,q)\right]=1.} $$

This is an exact single-variable formulation of the CWD scale-normalization problem.

## 7. Existence and uniqueness of the normalized CWD

Let $Y$ be a unit-rate exponentially distributed random variable:

$$ \Pr(Y>y)=e^{-y},\qquad y\geq0. $$

The normalized height has the representation

$$ Z=\tau\begin{cases}(Y/q)^{1/k_1},&Y\leq q,\\(Y/q)^{1/k_2},&Y>q.\end{cases} $$

Its second moment is

$$ \operatorname E[Z^2]=\tau^2U(q). $$

For fixed $Y>0$, the transformed height decreases strictly as $q$ increases. Hence $U(q)$ is continuous and strictly decreasing.

Moreover,

$$ \lim_{q\to0^+}U(q)=+\infty. $$

$$ \lim_{q\to\infty}U(q)=0. $$

Therefore, for every $\tau>0$, the equation

$$ \tau^2U(q)=1 $$

has exactly one positive solution.

The RMS-normalized Composite Weibull Distribution is therefore mathematically unique for a prescribed positive transitional height and fixed positive Weibull exponents.

### 7.1 Derivative of the normalization function

Define

$$ U_1=q^{-2/k_1}\gamma\left(1+\frac2{k_1},q\right). $$

$$ U_2=q^{-2/k_2}\Gamma\left(1+\frac2{k_2},q\right). $$

Differentiation gives

$$ \boxed{U'(q)=-\frac2q\left(\frac{U_1}{k_1}+\frac{U_2}{k_2}\right)<0.} $$

Implicit differentiation of $\tau^2U(q)=1$ yields

$$ \frac{dq}{d\tau}=-\frac{2U(q)}{\tau U'(q)}>0. $$

Increasing the normalized transitional height therefore moves the transition to progressively smaller exceedance probabilities.


### 7.2 Newton–Raphson solution in logarithmic coordinates

The normalization condition is solved for each prescribed positive normalized transition

$$
\tau=\frac{H_{tr}}{H_{rms}}>0.
$$

Introduce the logarithmic unknown

$$
z=\ln q,\qquad q=e^z>0.
$$

Writing $U(q)=U_1(q)+U_2(q)$ as in Section 7.1, the logarithmic residual is

$$
\boxed{R(z;\tau)=2\ln\tau+\ln\left[U_1(e^z)+U_2(e^z)\right].}
$$

The normalization equation is equivalent to $R(z;\tau)=0$.

#### Explicit starting approximation

The initial value for the numerical iteration is

$$
\boxed{
q_0(\tau)=1.002054709\,\tau^2
\tanh\left(\tau^{1.652259933}\right)
\exp\left[-\tanh(\tanh\tau)\exp(-2.146235748\,\tau)\right].
}
$$

The starting logarithmic iterate is

$$
z_0=\ln q_0(\tau).
$$

All constants in this expression are fixed numerical values. The approximation provides an initial iterate; it does not replace the exact RMS closure.

#### Analytical derivative

The derivatives of the two moment terms with respect to $z$ are

$$
\frac{dU_1}{dz}=-\frac{2}{k_1}U_1+q e^{-q},
$$

$$
\frac{dU_2}{dz}=-\frac{2}{k_2}U_2-q e^{-q}.
$$

The $q e^{-q}$ contributions cancel. Consequently,

$$
\frac{dU}{dz}=-\frac{2}{k_1}U_1-\frac{2}{k_2}U_2,
$$

and

$$
\boxed{R'(z;\tau)=-\frac{(2/k_1)U_1+(2/k_2)U_2}{U_1+U_2}.}
$$

For $k_1=2$ and $k_2=3.6$, this simplifies to

$$
\boxed{R'(z;\tau)=-\frac{U_1+\frac{5}{9}U_2}{U_1+U_2}.}
$$

Both moment contributions are positive when $q>0$. Thus

$$
-1<R'(z;\tau)<-\frac{5}{9}.
$$

The derivative is everywhere finite, negative and bounded away from zero over the positive domain.

#### Newton iteration and numerical acceptance

Newton–Raphson advances the logarithmic unknown according to

$$
\boxed{z_{n+1}=z_n-\frac{R(z_n;\tau)}{R'(z_n;\tau)},\qquad q_n=e^{z_n}.}
$$

Equivalently, in terms of the two moment contributions,

$$
z_{n+1}=z_n+R(z_n;\tau)
\frac{U_1+U_2}{U_1+\frac59U_2}.
$$

The numerical implementation uses a maximum of 40 Newton steps and accepts a root when

$$
\boxed{|R(z_n;\tau)|\leq5\times10^{-13}.}
$$

It rejects non-finite iterates and non-positive or non-finite starting values. After convergence, the two scale parameters are reconstructed from $q=e^z$, and the original second-moment and continuity conditions must separately satisfy

$$
|R_m|\leq10^{-10},\qquad |R_c|\leq10^{-10}.
$$

The gamma-series and continued-fraction implementations use a relative tolerance of $2.5\times10^{-14}$ and a maximum of 10,000 iterations per special-function evaluation. No numerical differentiation or solution of a two-variable system is required.

### 7.3 Analytical boundary at zero transition

The original equation contains powers of $q$ and cannot be evaluated by direct substitution of $\tau=q=0$. Its positive-side limit is well defined. Using the leading term $U(q)\sim\Gamma(14/9)q^{-5/9}$ gives

$$
\boxed{q(\tau)\sim\left[\Gamma(14/9)\tau^2\right]^{9/5}\quad(\tau\to0^+).}
$$

Hence $q(0)=0$ in the sense of continuous extension. The normalized scale parameters satisfy

$$
\lim_{\tau\to0^+}a_2=\frac{1}{\sqrt{\Gamma(14/9)}},\qquad
\lim_{\tau\to0^+}a_1=+\infty.
$$

At $\tau=0$, $z=\ln q$ is undefined. This boundary value is therefore obtained from the analytic limit rather than from Newton iteration. For the point calculator, positive $H_{m0}$, $d$, and $M$ give a strictly positive $\tau$.

## 8. CWD exceedance thresholds

The CWD exceedance threshold satisfies

$$ Q_C(H_N)=\frac1N. $$

Taking logarithms gives

$$ -\ln Q_C(H_N)=\ln N. $$

### 8.1 Threshold below the transition

When $\ln N\leq q$,

$$ \boxed{H_N=H_1(\ln N)^{1/k_1}.} $$

For $k_1=2$,

$$ H_N=H_1\sqrt{\ln N}. $$

### 8.2 Threshold above the transition

When $\ln N>q$,

$$ \boxed{H_N=H_2(\ln N)^{1/k_2}.} $$

With $k_2=3.6$,

$$ H_N=H_2(\ln N)^{5/18}. $$

The two expressions coincide at $N=e^q$.

## 9. Mean heights of the highest wave fractions

### 9.1 Mean involving both Weibull branches

If $\ln N<q$, the threshold lies below $H_{tr}$ and the selected population crosses the transition:

$$ H_{1/N}=N\int_{H_N}^{H_{tr}}hf_1(h)\,dh+N\int_{H_{tr}}^\infty hf_2(h)\,dh. $$

Using $u=(h/H_i)^{k_i}$ gives

$$ \boxed{H_{1/N}=NH_1[\Gamma(1+1/k_1,\ln N)-\Gamma(1+1/k_1,q)]+NH_2\Gamma(1+1/k_2,q).} $$

### 9.2 Mean entirely in the upper branch

When $\ln N\geq q$,

$$ \boxed{H_{1/N}=NH_2\Gamma\left(1+\frac1{k_2},\ln N\right).} $$

For $k_2=18/5$,

$$ \boxed{H_{1/N}=NH_2\Gamma(23/18,\ln N).} $$

The two conditional-mean expressions are continuous at $\ln N=q$.

### 9.3 Upper-tail asymptotic behaviour

For large $x$ and fixed $s$,

$$ \Gamma(s,x)\sim x^{s-1}e^{-x}\left[1+\frac{s-1}{x}+\frac{(s-1)(s-2)}{x^2}+\cdots\right]. $$

For a Weibull branch of exponent $k$,

$$ \frac{H_{1/N}}{H_N}\sim1+\frac1{k\ln N}+\frac{(1/k)(1/k-1)}{(\ln N)^2}+\cdots. $$

Consequently,

$$ \boxed{\frac{H_{1/N}}{H_N}\longrightarrow1\qquad(N\to\infty).} $$

This limiting relationship concerns the mathematical Weibull tail and does not imply unlimited validity of extreme extrapolation in physical seas.

## 10. Physical parameterization of the CWD

### 10.1 Spectral variance and RMS wave height

The spectral variance is related to the local spectral significant height by

$$ m_0=\left(\frac{H_{m0}}4\right)^2. $$

Battjes and Groenendijk adopted

$$ \boxed{H_{rms}=\left(2.69+3.24\frac{\sqrt{m_0}}d\right)\sqrt{m_0}.} $$

Define

$$ \mu=\frac{\sqrt{m_0}}d=\frac{H_{m0}}{4d}. $$

Then

$$ \frac{H_{rms}}d=(2.69+3.24\mu)\mu. $$

### 10.2 Transitional wave height

The empirical transitional wave height is

$$ \boxed{H_{tr}=(0.35+5.8\tan\alpha)d.} $$

For a slope $1:M$,

$$ \tan\alpha=\frac1M. $$

Therefore,

$$ H_{tr}=\left(0.35+\frac{5.8}{M}\right)d. $$

### 10.3 Governing dimensionless parameter

Combining the preceding relations gives

$$ \boxed{\tau=\frac{H_{tr}}{H_{rms}}=\frac{0.35+5.8/M}{(2.69+3.24\mu)\mu}.} $$

This expression establishes the relationship between the physical sea-state parameters and the normalized CWD shape.

For fixed slope, increasing $H_{m0}/d$ decreases $\tau$ and moves a greater proportion of the distribution into the upper Weibull branch.

For fixed $H_{m0}$ and $d$, increasing the foreshore slope increases $H_{tr}$ and $\tau$.

These relationships are empirical properties of the shallow-foreshore parameterization, not universal laws of wave transformation.

## 11. Limiting behaviour and the Rayleigh relationship

### 11.1 Small transitional height

As $\tau\to0^+$, the upper Weibull branch dominates. The analytic boundary value $q(0)=0$ is defined by the positive-side limit in Section 7.3, not by direct evaluation of the singular scalar equation.

$$ U(q)\sim q^{-2/k_2}\Gamma\left(1+\frac2{k_2}\right). $$

Consequently,

$$ q\sim\left[\tau^2\Gamma\left(1+\frac2{k_2}\right)\right]^{k_2/2}. $$

The normalized upper scale tends to

$$ \boxed{a_2\longrightarrow\left[\Gamma\left(1+\frac2{k_2}\right)\right]^{-1/2}.} $$

The limiting distribution is a pure Weibull law with exponent $k_2$ and unit second moment.

### 11.2 Large transitional height

As $\tau\to\infty$, the lower Weibull branch dominates.

For $k_1=2$,

$$ q\sim\tau^2. $$

The normalized lower scale tends to unity:

$$ a_1\longrightarrow1. $$

Therefore, for every fixed normalized height $z$,

$$ \boxed{Q_C(H_{rms}z)\longrightarrow e^{-z^2}.} $$

The RMS-normalized CWD thus approaches the Rayleigh distribution as the transition moves into the far upper tail.

### 11.3 Spectral and individual-wave normalization

The empirical CWD relation tends to

$$ H_{rms}\longrightarrow2.69\sqrt{m_0} $$

as the degree of saturation approaches zero.

Ideal narrow-band Rayleigh theory gives

$$ R_R=\sqrt8\sqrt{m_0}. $$

Their limiting ratio is

$$ \boxed{\frac{H_{rms}}{R_R}\longrightarrow\frac{2.69}{\sqrt8}\simeq0.9511.} $$

Hence, convergence of the normalized CWD shape to Rayleigh does not imply identical dimensional wave-height statistics at fixed spectral significant height.

## 12. Worked calculations

### 12.1 Shallow-foreshore case

Consider

$$ H_{m0}=2.5\ \mathrm{m},\qquad d=5.0\ \mathrm{m},\qquad M=100. $$

The spectral standard deviation is

$$ \sqrt{m_0}=\frac{2.5}{4}=0.625\ \mathrm{m}. $$

The degree of saturation is

$$ \mu=\frac{0.625}{5}=0.125. $$

The empirical RMS wave height is

$$ H_{rms}=(2.69+3.24\times0.125)\times0.625=1.934375\ \mathrm{m}. $$

The transitional wave height is

$$ H_{tr}=(0.35+5.8/100)\times5=2.0400\ \mathrm{m}. $$

Therefore,

$$ \tau=\frac{2.0400}{1.934375}=1.0546042. $$

The scalar RMS equation gives approximately

$$ q=0.8312678123. $$

For this case, the explicit initial approximation gives $q_0\simeq0.8303316764$. Newton–Raphson reaches the stated scalar root in two steps, evaluated with the analytical derivative and the logarithmic residual.

The normalized Weibull scales are

$$ a_1=1.1566954161,\qquad a_2=1.1101551380. $$

The dimensional scales are

$$ H_1=2.237483\ \mathrm{m},\qquad H_2=2.147456\ \mathrm{m}. $$

The transition probability is

$$ Q_C(H_{tr})=e^{-q}\simeq0.435497. $$

Because $\ln3>q$, the selected characteristic thresholds all lie in the upper branch.

| $N$ | Threshold $H_N$ (m) | Mean $H_{1/N}$ (m) |
|---:|---:|---:|
| 3 | 2.204297 | 2.591137 |
| 10 | 2.707321 | 2.970068 |
| 50 | 3.136753 | 3.329627 |
| 100 | 3.282156 | 3.456711 |
| 250 | 3.451838 | 3.607696 |
| 1000 | 3.673445 | 3.808553 |

These are statistics of the continuous CWD itself.

### 12.2 Case with mixed-branch integration

Retain $H_{m0}=2.5$ m and slope $1:100$, but adopt $d=8.0$ m.

The characteristic values become

$$ H_{rms}=1.839453125\ \mathrm{m},\qquad H_{tr}=3.2640\ \mathrm{m}. $$

The corresponding normalized transition and transition exponent are

$$ \tau\simeq1.7744404332,\qquad q\simeq3.0756452504. $$

The dimensional Weibull scales are

$$ H_1\simeq1.861153\ \mathrm{m},\qquad H_2\simeq2.388979\ \mathrm{m}. $$

Here,

$$ \ln3<q,\qquad\ln10<q,\qquad\ln50>q. $$

The means of the highest one-third and one-tenth must therefore integrate both Weibull branches. The smaller fractions lie entirely in the upper branch.

| $N$ | Integration domain | $H_{1/N}$ (m) |
|---:|---|---:|
| 3 | Both branches | 2.604212 |
| 10 | Both branches | 3.247634 |
| 50 | Upper branch | 3.704107 |
| 100 | Upper branch | 3.845484 |
| 250 | Upper branch | 4.013451 |
| 1000 | Upper branch | 4.236898 |

## 13. Mathematical and physical limitations

The Rayleigh law is derived under linear, Gaussian, narrow-banded assumptions. It does not exactly describe every broadband or nonlinear sea state.

The CWD is an empirical representation of individual wave-height distributions on shallow foreshores. Its mathematical definition is internally consistent, but physical applicability depends on the assumptions underlying the empirical parameterization.

The original shallow-foreshore investigations were based primarily on relatively simple, gradually varying laboratory profiles, with slopes approximately between $1:20$ and $1:250$.

Particular caution is necessary for extensive constant-depth beds, abrupt bathymetric transitions, barred profiles and situations involving substantial wave-current interaction or strong multidirectional effects.

The mathematical CWD has no finite upper endpoint. The transition height is not itself a breaker-height criterion.

Finally, the distribution of wave heights alone does not determine wave periods, crest elevations, wave-induced velocities or full hydrodynamic loads.

## 14. Conclusions

The Rayleigh distribution follows from a Gaussian wave-envelope construction under narrow-band conditions. Its moments, exceedance thresholds and highest-fraction means admit analytical expressions involving the complete and incomplete gamma functions.

The Composite Weibull Distribution extends this mathematical structure through two Weibull branches whose exponents represent different portions of the individual-wave population.

Continuity at the transition gives one common exponent $q$, while RMS normalization determines this parameter uniquely. The entire normalized distribution is consequently governed by the normalized transitional height $\tau$.

The essential CWD statistics—probabilities, quantiles, moments and averages of the highest wave fractions—follow from the same two-branch law without introducing further independent scale parameters.

The distinction between spectral significant height, individual-wave RMS, exceedance threshold and upper-fraction mean remains essential throughout:

$$ \boxed{H_{m0},\quad H_{rms},\quad H_N,\quad H_{1/N}\quad\mathrm{are\ distinct\ statistical\ quantities}.} $$

## References

1. Battjes, J. A., and Groenendijk, H. W. (2000). *Wave height distributions on shallow foreshores*. Coastal Engineering, 40, 161–182. DOI: 10.1016/S0378-3839(00)00007-7.
2. Groenendijk, H. W. (1998). *Shallow foreshore wave height statistics*. Delft University of Technology.
3. Groenendijk, H. W., and van Gent, M. R. A. (1998). *Shallow foreshore wave height statistics: A predictive model for the probability of exceedance of wave heights*. WL | Delft Hydraulics, Report H3351.
4. Longuet-Higgins, M. S. (1952). *On the statistical distribution of the heights of sea waves*. Journal of Marine Research, 11(3), 245–266.
5. Longuet-Higgins, M. S. (1980). *On the distribution of the heights of sea waves: Some effects of nonlinearity and finite bandwidth*. Journal of Geophysical Research, 85(C3), 1519–1523.
6. Caires, S., and van Gent, M. R. A. (2012). *Wave height distribution in constant and finite depths*. Proceedings of the 33rd International Conference on Coastal Engineering.
7. NIST Digital Library of Mathematical Functions. Chapters 5 and 8: Gamma functions, incomplete gamma functions, recurrence relations and asymptotic expansions.
