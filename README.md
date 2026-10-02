# Wave Height Distributions on Shallow Foreshores

This document describes a multi-language engineering calculator for the local distribution of individual wave heights on shallow foreshores, based primarily on the **Composite Weibull Distribution** of Battjes and Groenendijk (2000).

The project computes characteristic individual-wave statistics from three local inputs:

- spectral significant wave height, $H_{m0}$;
- local still-water depth, $d$;
- foreshore slope written as $1:M$, with $\tan\alpha=1/M$.

The same computational model is provided as a C++ command-line program, a native Windows C++ graphical interface, a Fortran command-line program, a MATLAB function, an interactive Jupyter notebook, and an Android app.

The implementation includes the deep-water convergence treatment discussed by Caires and Van Gent (2012): the Battjes-Groenendijk result is not permitted to exceed the corresponding Rayleigh value, and the program switches directly to Rayleigh statistics when the normalized transitional wave height is sufficiently large.

---

## Scope of the calculator

The calculator predicts the **local short-term distribution of individual zero-crossing wave heights** at a point on a shallow, gently sloping foreshore. It returns:

- free-surface variance $m_0$;
- root-mean-square individual wave height $H_{rms}$;
- transitional wave height $H_{tr}$;
- normalized transition $\widetilde H_{tr}=H_{tr}/H_{rms}$;
- Composite Weibull scale parameters $H_1$ and $H_2$ when the Battjes-Groenendijk branch is active;
- means of the highest $1/3$, $1/10$, $1/50$, $1/100$, $1/250$, and $1/1000$ fractions of the individual waves;
- diagnostic ratios relative to $H_{1/3}$;
- the selected distribution branch: `B&G` or `Rayleigh`.

The program is a **point model**. It does not propagate a spectrum, solve wave transformation along a profile, calculate wave setup, calculate breaking dissipation, or resolve individual waves in time. The local $H_{m0}$, depth, and slope must already be known from measurements or from an appropriate wave-transformation model.

---

## Historical development

### Chronological overview

| Period | Development | Relevance to this software |
|---|---|---|
| 1880s | Rayleigh developed the probability distribution now bearing his name in the context of random vibrations and wave-like phenomena. | Provides the one-parameter reference distribution for linear random-wave heights. |
| 1952 | Longuet-Higgins connected Gaussian narrow-band sea-surface theory to the Rayleigh distribution of crest-to-trough wave heights. | Establishes the deep-water statistical baseline and the standard exceedance formulas. |
| 1978 | Battjes and Janssen formulated a random-wave breaking dissipation model based on an energy balance and bore-type dissipation. | Clarified that depth-induced breaking acts statistically on a population of random waves and selectively removes energy from the largest waves. |
| 1979–1980 | Goda and Longuet-Higgins documented finite-bandwidth and nonlinearity effects in measured wind-wave statistics. | Explains why the empirical deep-water value of $H_{rms}/\sqrt{m_0}$ used by Battjes and Groenendijk differs from the ideal narrow-band value. |
| 1980s–1990s | Several shallow-water distributions were proposed, including transformed-Rayleigh, Glukhovskiy-type, modified Glukhovskiy, and Weibull formulations. | Demonstrated the need for depth-dependent upper-tail models, but single-shape distributions did not fully reproduce the observed change between low and high waves. |
| 1998 | Groenendijk's MSc thesis and the WL \| Delft Hydraulics H3351 report developed and validated the Composed Weibull point model. | Established the two-branch distribution, its calibration database, the slope-dependent transition, and the calculation recipe. |
| 2000 | Battjes and Groenendijk published the Composite Weibull Distribution in *Coastal Engineering*. | Provides the principal scientific model implemented by this repository. |
| 2012 | Caires and Van Gent examined finite-depth and constant-depth behaviour and clarified deep-water convergence and flat-bottom limitations. | Motivates the Rayleigh switch, Rayleigh caps, and explicit warnings for constant-depth applications. |
| Present implementation | The same equations are implemented in C++, Fortran, MATLAB, Python/Jupyter, and the Android app. | Provides reproducible cross-language calculations and engineering reporting. |

### From Gaussian surface elevation to individual wave heights

The historical transition from linear deep-water statistics to the Composite Weibull model is best understood as a sequence of increasingly realistic descriptions.

A linear, stationary, narrow-band random sea may be represented by a Gaussian surface-elevation process. Under those assumptions, the slowly varying envelope amplitude is Rayleigh-distributed, and the corresponding zero-crossing wave heights are also described by a Rayleigh law. This result is not merely an empirical curve fit: it follows from the statistics of two independent Gaussian quadrature components.

Finite spectral bandwidth, nonlinear crest-trough asymmetry, triad interactions, wave breaking, and dissipation progressively weaken that theoretical foundation. On a shallow foreshore, the principal distortion is not uniform across all wave heights. Smaller waves remain comparatively close to the Rayleigh population, whereas the largest waves are preferentially limited by breaking. The observed distribution therefore changes curvature in the upper tail. This is the central empirical observation behind the two-branch Composite Weibull representation.

### Relation to random-wave transformation models

The calculator is intentionally narrower in scope than a spectral wave-transformation model. Historical models such as Battjes-Janssen calculate the spatial evolution of wave energy and breaking dissipation. The present software starts after that transformation problem has been solved: it accepts the **local** $H_{m0}$, local depth, and local slope and reconstructs the corresponding local distribution of individual wave heights.

The distinction is fundamental:

- a transformation model determines how $m_0$ changes along a profile;
- the Composite Weibull point model determines how the local energy is distributed among individual wave heights at one selected point;
- a structural, run-up, or overtopping method subsequently uses one or more of those local wave-height statistics.

The three stages should not be conflated.


### Rayleigh distribution and early random-wave statistics

The Rayleigh distribution originated in nineteenth-century probability and vibration theory. Longuet-Higgins (1952) established its central role in ocean-wave statistics by showing that crest-to-trough wave heights in a narrow-banded linear Gaussian sea follow a Rayleigh distribution.

For deep-water or weakly depth-limited conditions, the Rayleigh cumulative distribution is

$$
F_R(H)=1-\exp\left[-\left(H/H_{rms}\right)^2\right]
$$


The corresponding exceedance probability is

$$
P(H>h)=\exp\left[-\left(h/H_{rms}\right)^2\right]
$$


This distribution is completely defined by one scale parameter. Consequently, all characteristic wave heights are related by fixed ratios.

### Random-wave breaking models

Battjes and Janssen (1978) introduced a physically based energy-dissipation model for random waves breaking on a beach. Their work combined an energy balance, bore-type dissipation, and a depth-limited representation of the breaking-wave population. Although that model is not the Composite Weibull model implemented here, it established essential concepts used in later shallow-water random-wave modelling:

- depth-induced breaking selectively affects the largest waves;
- the local depth limits the upper part of the wave-height distribution;
- wave-energy dissipation and wave setup can be treated through averaged conservation equations;
- random breaking requires a statistical description rather than a single deterministic breaker height.

### Modified shallow-water distributions

During the following decades, several alternatives to the Rayleigh distribution were proposed for finite and shallow depths, including Glukhovskiy-type, modified Glukhovskiy, Weibull, and transformed-Rayleigh models. Klopman (1996) developed a modified Glukhovskiy formulation intended to account for depth limitation using local wave energy and depth.

These single-form distributions improved some shallow-water predictions but did not reproduce the observed abrupt change in slope between the lower and upper portions of measured wave-height distributions on Rayleigh probability paper.

### Groenendijk thesis and Delft Hydraulics report, 1998

Groenendijk's 1998 Delft University of Technology thesis, followed by the WL | Delft Hydraulics report H3351 by Groenendijk and Van Gent, developed the **Composed Weibull Distribution**. The terminology in the thesis and report is commonly written *Composed Weibull*; the later journal paper uses *Composite Weibull*. Both names refer to the same two-branch concept.

The model was calibrated and validated using laboratory data for plane foreshores with slopes from approximately $1:20$ to $1:250$. The principal observation was that measured distributions contain:

- a lower-wave region that remains approximately Rayleigh-shaped;
- an upper-wave region in which breaking suppresses the tail much more strongly;
- a distinct transitional height separating the two regions.

### Battjes and Groenendijk, 2000

Battjes and Groenendijk (2000) published the model in *Coastal Engineering*. The paper fixed the Weibull exponents at

$$
k_1=2.0
$$


$$
k_2=3.6
$$


and parameterized the transitional height using local depth and foreshore slope. With a continuity condition and a second-moment constraint, the entire normalized distribution is controlled by the single parameter $\widetilde H_{tr}$.

### Caires and Van Gent, 2012

Caires and Van Gent (2012) examined implementation details and the use of the Battjes-Groenendijk distribution in constant and finite depths. Two conclusions are directly relevant to this software:

1. When dimensional wave heights are recovered using the Battjes-Groenendijk $H_{rms}$ parameterization, the calculated values can overshoot the Rayleigh limit before approaching an incorrect lower asymptote in deep water.
2. On shallow **flat bottoms**, which are outside the original calibration domain, the Battjes-Groenendijk distribution with a nominal $1:250$ slope underestimated measured high-wave quantiles by approximately 7% to 15% on average.

The present implementation therefore applies a Rayleigh switch and Rayleigh caps, and it explicitly treats constant-depth flat-bottom applications as outside the recommended model domain.

---

## Physical basis

### Why Rayleigh statistics change on a shallow foreshore

The Rayleigh model relies on an approximately linear, Gaussian, narrow-banded sea surface. These assumptions progressively deteriorate as waves propagate into shallow water.

**Shoaling and nonlinear interactions.** Interaction with the seabed increases nonlinearity. Triad interactions transfer energy between harmonics, producing sharper crests, flatter troughs, and non-Gaussian surface elevations.

**Selective depth-induced breaking.** The highest waves break first. Breaking removes energy preferentially from the upper tail of the individual-wave distribution, while smaller waves may continue to propagate with less modification.

**Spatial adaptation.** Breaking and statistical reshaping do not occur instantaneously at the depth where a breaker criterion is first reached. The wave field requires a finite propagation distance to adapt. This spatial lag explains why the foreshore slope appears in the transitional-height parameterization.

### Two statistical populations

The Composite Weibull model represents the observed distribution using two Weibull branches:

- the first branch represents lower, substantially non-breaking waves;
- the second branch represents higher, breaking or recently broken waves.

The branches meet at $H_{tr}$, the transitional wave height.

---

## Notation

| Symbol | Meaning | Units |
|---|---|---:|
| $H$ | individual zero-crossing wave height | m |
| $H_{m0}$ | spectral significant wave height, $4\sqrt{m_0}$ | m |
| $m_0$ | zero-order spectral moment or free-surface variance | m² |
| $H_{rms}$ | root-mean-square individual wave height | m |
| $d$ | local still-water depth | m |
| $M$ | denominator of a $1:M$ foreshore slope | - |
| $\tan\alpha$ | local foreshore slope, $1/M$ | - |
| $H_{tr}$ | transitional wave height | m |
| $H_1$, $H_2$ | Weibull scale parameters | m |
| $k_1$, $k_2$ | Weibull shape exponents | - |
| $\widetilde H$ | wave height normalized by $H_{rms}$ | - |
| $H_N$ | height exceeded with probability $1/N$ | m |
| $H_{1/N}$ | mean height of the highest $1/N$ fraction of waves | m |
| $\gamma(a,x)$ | unnormalized lower incomplete gamma function | - |
| $\Gamma(a,x)$ | unnormalized upper incomplete gamma function | - |

A tilde denotes normalization with $H_{rms}$:

$$
\widetilde H=H/H_{rms}
$$


The notation $H_N$ and $H_{1/N}$ must not be confused:

- $H_N$ is an exceedance threshold;
- $H_{1/N}$ is the mean of all waves above that threshold.

The software reports $H_{1/3}$, $H_{1/10}$, $H_{1/50}$, $H_{1/100}$, $H_{1/250}$, and $H_{1/1000}$. It does not report a deterministic maximum wave height.

---

## Probability definitions and statistical conventions

### Cumulative distribution, exceedance probability, and density

For a non-negative individual wave height $H$, the cumulative distribution function is

$$
F_H(h)=P(H\leq h), \qquad h\geq0.
$$


The survival or exceedance function is

$$
Q_H(h)=P(H>h)=1-F_H(h).
$$


Where the distribution is differentiable, the probability density is

$$
f_H(h)=\frac{dF_H(h)}{dh}.
$$


The density integrates to unity:

$$
\int_0^\infty f_H(h)\,dh=1.
$$


### Exceedance level $H_N$

The level $H_N$ is exceeded, on average, once in every $N$ individual waves. It is defined by

$$
Q_H(H_N)=\frac{1}{N},
$$


or equivalently,

$$
F_H(H_N)=1-\frac{1}{N}.
$$


$H_N$ is a quantile or threshold. It is **not** the average of the highest waves.

### Mean of the highest $1/N$ fraction, $H_{1/N}$

The mean of the highest $1/N$ fraction is the conditional mean above $H_N$:

$$
H_{1/N}=E\left[H\mid H>H_N\right].
$$


Because $P(H>H_N)=1/N$,

$$
H_{1/N} = N\int_{H_N}^{\infty}h\,f_H(h)\,dh.
$$


This definition is used throughout the code. For example:

- $H_{1/3}$ is the arithmetic mean of the highest one-third of the individual waves;
- $H_{1/100}$ is the arithmetic mean of the highest one percent;
- $H_{1/1000}$ is the arithmetic mean of the highest 0.1 percent.

Neither $H_N$ nor $H_{1/N}$ is a deterministic storm maximum. A storm maximum additionally depends on the number of waves, dependence between successive waves, sea-state duration, and nonstationarity.

### Spectral and wave-by-wave quantities

The zero-order spectral moment is

$$
m_0=\int_0^\infty S_\eta(f)\,df,
$$


where $S_\eta(f)$ is the variance-density spectrum of free-surface elevation. The spectral significant wave height is

$$
H_{m0}=4\sqrt{m_0}.
$$


$H_{m0}$ is obtained from the spectrum, while $H_{1/3}$ is obtained from an ordered sample or a probability distribution of individual zero-crossing wave heights. They are close under many deep-water conditions but are not identical by definition and need not remain close in shallow water.

---

## Rayleigh distribution

### Cumulative distribution and density

The Rayleigh cumulative distribution is

$$
F_R(H)=1-\exp\left[-\left(H/H_{rms}\right)^2\right]
$$


The probability density is

$$
f_R(H)=\frac{2H}{H_{rms}^2}\exp\left[-\left(H/H_{rms}\right)^2\right]
$$


### Exceedance wave height

The wave height exceeded once, on average, in every $N$ waves is obtained from $P(H>H_N)=1/N$:

$$
H_N=H_{rms}\sqrt{\ln N}
$$


### Mean of the highest fraction

For a Rayleigh distribution, the mean of the highest $1/N$ fraction is

$$
H_{1/N}=N H_{rms}\Gamma\left(\frac{3}{2},\ln N\right)
$$


For the narrow-band relation $H_{rms}=\sqrt{8m_0}$ and $H_{m0}=4\sqrt{m_0}$, the exact ratios used by the software are:

| Statistic | Rayleigh ratio to $H_{m0}$ |
|---|---:|
| $H_{1/3}$ | 1.001075736951740 |
| $H_{1/10}$ | 1.272734273369137 |
| $H_{1/50}$ | 1.560113379974762 |
| $H_{1/100}$ | 1.668233372358517 |
| $H_{1/250}$ | 1.801017222497626 |
| $H_{1/1000}$ | 1.984835590575388 |

These constants are used both in the direct Rayleigh branch and as upper caps on dimensional Battjes-Groenendijk results.

### Rayleigh moments and scale relations

For the Rayleigh formulation used here,

$$
F_R(H)=1-\exp\left[-\left(\frac{H}{H_{rms}}\right)^2\right],
$$


the general non-central moment is

$$
E\left[H^r\right] = H_{rms}^{\,r} \Gamma\left(1+\frac{r}{2}\right).
$$


For $r=2$,

$$
E\left[H^2\right]=H_{rms}^2,
$$


which confirms that the scale appearing in the selected Rayleigh form is the root-mean-square individual wave height.

Under ideal narrow-band Gaussian theory,

$$
H_{rms}=\sqrt{8m_0},
$$


and therefore

$$
H_{rms}=\frac{H_{m0}}{\sqrt{2}}.
$$


The empirical Battjes-Groenendijk relationship used elsewhere in the calculator is intentionally different because it incorporates the observed finite-bandwidth dependence of individual wave heights.

### Derivation of the Rayleigh mean-high-wave expression

The threshold exceeded once in $N$ waves satisfies

$$
\exp\left[-\left(\frac{H_N}{H_{rms}}\right)^2\right]=\frac{1}{N}.
$$


Hence,

$$
\left(\frac{H_N}{H_{rms}}\right)^2=\ln N,
$$


and

$$
H_N=H_{rms}\sqrt{\ln N}.
$$


Substitution into the conditional-mean integral gives

$$
H_{1/N} = N\int_{H_N}^{\infty} h\, \frac{2h}{H_{rms}^2} \exp\left[-\left(\frac{h}{H_{rms}}\right)^2\right]dh.
$$


Using the change of variable

$$
u=\left(\frac{h}{H_{rms}}\right)^2,
$$


the result is

$$
H_{1/N} = N H_{rms} \Gamma\left(\frac{3}{2},\ln N\right).
$$


This exact expression is the source of the Rayleigh constants tabulated above.

---

## Composite Weibull Distribution

### Normalized variables

The governing equations are solved in terms of wave heights normalized by $H_{rms}$:

$$
\widetilde H=\frac{H}{H_{rms}}, \qquad \widetilde H_1=\frac{H_1}{H_{rms}}, \qquad \widetilde H_2=\frac{H_2}{H_{rms}}, \qquad \widetilde H_{tr}=\frac{H_{tr}}{H_{rms}}.
$$


Normalization removes the absolute dimensional scale from the distribution-shape problem. Once the normalized solution is known, dimensional values are recovered by multiplication by $H_{rms}$.

### Survival function

The exceedance form is often more convenient than the CDF:

$$
Q(H)=1-F(H).
$$


For the Composite Weibull Distribution,

$$
Q(H)= \begin{cases} \exp\left[-\left(\dfrac{H}{H_1}\right)^{k_1}\right], & H\leq H_{tr}, \\[8pt] \exp\left[-\left(\dfrac{H}{H_2}\right)^{k_2}\right], & H>H_{tr}. \end{cases}
$$


Continuity of the CDF automatically implies continuity of the exceedance probability at the transition.

### Probability associated with the transition

Define

$$
x_1= \left(\frac{\widetilde H_{tr}}{\widetilde H_1}\right)^{k_1}, \qquad x_2= \left(\frac{\widetilde H_{tr}}{\widetilde H_2}\right)^{k_2}.
$$


The continuity condition requires $x_1=x_2$. Denoting their common value by $x_{tr}$,

$$
x_{tr}=x_1=x_2.
$$


The exceedance probability at the transition is

$$
P(H>H_{tr})=\exp(-x_{tr}),
$$


and the nominal recurrence count associated with the transition is

$$
N_{tr}=\exp(x_{tr}).
$$


For $N<N_{tr}$, the exceedance threshold lies below the transition and the mean-high-wave integral includes both branches. For $N\geq N_{tr}$, the threshold lies in the upper branch.

### Scale-parameter relation implied by continuity

The two scale parameters are not independent. From

$$
\left(\frac{H_{tr}}{H_1}\right)^{k_1} = \left(\frac{H_{tr}}{H_2}\right)^{k_2},
$$


one may write

$$
H_2 = H_{tr} \left(\frac{H_1}{H_{tr}}\right)^{k_1/k_2}.
$$


The implementation nevertheless solves for both normalized scales simultaneously. This preserves the direct two-equation formulation and provides an explicit continuity residual for convergence checking.

### General moment of the Composite Weibull Distribution

For any real moment order $r>-k_1$, the normalized piecewise moment can be expressed as

$$
E\left[\widetilde H^r\right] = \widetilde H_1^{\,r} \gamma\left(1+\frac{r}{k_1},x_1\right) + \widetilde H_2^{\,r} \Gamma\left(1+\frac{r}{k_2},x_2\right).
$$


The first term integrates the lower branch from zero to the transition. The second term integrates the upper branch from the transition to infinity.

For $r=0$, the expression represents total probability and equals one when continuity is satisfied. For $r=1$, it gives the normalized mean individual wave height. For $r=2$, it gives the second moment used to impose the $H_{rms}$ normalization.

### Piecewise cumulative distribution

The dimensional cumulative distribution is

$$
F(H)=\begin{cases}1-\exp\left[-\left(H/H_1\right)^{k_1}\right],&H\leq H_{tr},\\1-\exp\left[-\left(H/H_2\right)^{k_2}\right],&H>H_{tr}.\end{cases}
$$


The fixed shape exponents are

$$
k_1=2.0
$$


$$
k_2=3.6
$$


Because $k_1=2$, the lower branch has Rayleigh form. The larger exponent $k_2=3.6$ produces a more rapidly decaying upper tail.

### Piecewise probability density

The corresponding density is

$$
f(H)=\begin{cases}\frac{k_1}{H_1}\left(H/H_1\right)^{k_1-1}\exp\left[-\left(H/H_1\right)^{k_1}\right],&H\leq H_{tr},\\\frac{k_2}{H_2}\left(H/H_2\right)^{k_2-1}\exp\left[-\left(H/H_2\right)^{k_2}\right],&H>H_{tr}.\end{cases}
$$


The CDF is continuous at $H_{tr}$, but the density is generally discontinuous because $k_1\neq k_2$. This is an accepted empirical simplification: the integral statistics required for engineering calculations remain well behaved.

### Continuity condition

Continuity at the transition requires

$$
F_1(H_{tr})=F_2(H_{tr})
$$


which is equivalent to

$$
\left(H_{tr}/H_1\right)^{k_1}=\left(H_{tr}/H_2\right)^{k_2}
$$


In normalized form,

$$
\left(\widetilde H_{tr}/\widetilde H_1\right)^{k_1}=\left(\widetilde H_{tr}/\widetilde H_2\right)^{k_2}
$$


### Second-moment normalization

The normalized distribution must satisfy

$$
E\left[\widetilde H^2\right]=1
$$


Defining

$$
x_1=\left(\widetilde H_{tr}/\widetilde H_1\right)^{k_1}
$$


$$
x_2=\left(\widetilde H_{tr}/\widetilde H_2\right)^{k_2}
$$


the second-moment condition is

$$
\widetilde H_1^2\gamma\left(1+\frac{2}{k_1},x_1\right)+\widetilde H_2^2\Gamma\left(1+\frac{2}{k_2},x_2\right)=1
$$


The code solves the equivalent residual equation

$$
F_1=\sqrt{\widetilde H_1^2\gamma\left(1+\frac{2}{k_1},x_1\right)+\widetilde H_2^2\Gamma\left(1+\frac{2}{k_2},x_2\right)}-1=0
$$


The continuity residual is

$$
F_2=x_1-x_2=0
$$


Once $\widetilde H_{tr}$ is known, these two equations determine $\widetilde H_1$ and $\widetilde H_2$.

---

## Local physical parameterization

### Spectral variance

The input spectral significant wave height is converted to free-surface variance using

$$
m_0=\left(H_{m0}/4\right)^2
$$


### Root-mean-square individual wave height

The implementation uses the Battjes-Groenendijk empirical relation

$$
H_{rms}=\left(2.69+3.24\frac{\sqrt{m_0}}{d}\right)\sqrt{m_0}
$$


Equivalently,

$$
\frac{H_{rms}}{\sqrt{m_0}}=2.69+3.24\frac{\sqrt{m_0}}{d}
$$


The coefficient 2.69 represents the broad-banded deep-water limit adopted by Battjes and Groenendijk following Goda's analysis of wind-wave data. It differs from the narrow-band Rayleigh relation $H_{rms}/\sqrt{m_0}=\sqrt{8}\approx2.828427$.

This distinction is the reason dimensional Battjes-Groenendijk values require the deep-water convergence treatment described later.

### Foreshore slope

For an input slope $1:M$,

$$
\tan\alpha=1/M
$$


The program expects the positive denominator $M$, not the tangent itself.

### Transitional wave height

The dimensional transitional wave height is

$$
H_{tr}=\left(0.35+5.8\tan\alpha\right)d
$$


For a slope entered as $1:M$,

$$
H_{tr}=\left(0.35+5.8/M\right)d
$$


The normalized transition controlling the CWD shape is

$$
\widetilde H_{tr}=H_{tr}/H_{rms}
$$


A relatively small $\widetilde H_{tr}$ places more of the distribution in the breaking-controlled upper branch. A large $\widetilde H_{tr}$ moves the transition into the far tail and causes the normalized CWD to approach Rayleigh statistics.

### Degree of saturation

A useful dimensionless measure of local wave energy relative to depth is

$$
\mu=\frac{\sqrt{m_0}}{d}.
$$


Because $H_{m0}=4\sqrt{m_0}$,

$$
\mu=\frac{H_{m0}}{4d}.
$$


The empirical $H_{rms}$ relation can then be written as

$$
\frac{H_{rms}}{d} = \left(2.69+3.24\mu\right)\mu.
$$


Combining this expression with the transition-height formula gives a direct dimensionless expression for the governing transition:

$$
\widetilde H_{tr} = \frac{0.35+5.8/M} {\left(2.69+3.24\mu\right)\mu}.
$$


This equation shows explicitly that the normalized distribution shape is governed by two local dimensionless descriptors:

1. the energy-to-depth ratio $\mu$;
2. the foreshore slope denominator $M$.

For fixed slope, increasing $H_{m0}/d$ decreases $\widetilde H_{tr}$ and places more probability in the breaking-controlled upper branch. For fixed local energy and depth, a steeper foreshore increases $H_{tr}$ and therefore increases $\widetilde H_{tr}$, reflecting the finite spatial distance required for the wave population to adapt to breaking.

---

## Statistical wave-height calculations

### Wave height exceeded once in every $N$ waves

The candidate threshold from the lower branch is

$$
\widetilde H_{N,1}=\widetilde H_1\left(\ln N\right)^{1/k_1}
$$


The final exceedance height is

$$
\widetilde H_N=\begin{cases}\widetilde H_1\left(\ln N\right)^{1/k_1},&\widetilde H_{N,1}<\widetilde H_{tr},\\\widetilde H_2\left(\ln N\right)^{1/k_2},&\widetilde H_{N,1}\geq\widetilde H_{tr}.\end{cases}
$$


### Mean of the highest $1/N$ fraction

When $\widetilde H_N<\widetilde H_{tr}$, the averaging integral crosses both Weibull branches:

$$
\widetilde H_{1/N}=N\widetilde H_1\left[\Gamma\left(1+\frac{1}{k_1},\ln N\right)-\Gamma\left(1+\frac{1}{k_1},x_1\right)\right]+N\widetilde H_2\Gamma\left(1+\frac{1}{k_2},x_2\right)
$$


When $\widetilde H_N\geq\widetilde H_{tr}$, the averaging integral lies entirely in the upper branch:

$$
\widetilde H_{1/N}=N\widetilde H_2\Gamma\left(1+\frac{1}{k_2},\ln N\right)
$$


The dimensional result is

$$
H_{1/N}=\widetilde H_{1/N}H_{rms}
$$


The implementation evaluates these expressions for $N=3$, $10$, $50$, $100$, $250$, and $1000$.

---

### Unified branch criterion

The lower-branch candidate threshold is

$$
\widetilde H_{N,1} = \widetilde H_1(\ln N)^{1/k_1}.
$$


The same decision may be written using $x_{tr}$:

$$
\ln N < x_{tr} \quad\Longleftrightarrow\quad \widetilde H_N<\widetilde H_{tr}.
$$


Thus the branch is selected from a comparison between the requested exceedance level $\ln N$ and the transition exponent $x_{tr}$.

### Derivation of the mean-high-wave integral

For a Weibull branch with scale $H_s$ and exponent $k$,

$$
f(H) = \frac{k}{H_s} \left(\frac{H}{H_s}\right)^{k-1} \exp\left[-\left(\frac{H}{H_s}\right)^k\right].
$$


The substitution

$$
u=\left(\frac{H}{H_s}\right)^k
$$


gives

$$
H=H_su^{1/k}, \qquad dH=\frac{H_s}{k}u^{1/k-1}du.
$$


Consequently,

$$
\int Hf(H)\,dH = H_s\int u^{1/k}e^{-u}\,du.
$$


This integral is expressed by incomplete gamma functions. The two formulas used by the program follow by splitting the conditional-mean integral at $H_{tr}$ whenever the exceedance threshold is below the transition.

---

## Incomplete gamma functions

The complete gamma function is

$$
\Gamma(a)=\int_0^\infty t^{a-1}e^{-t}\,dt
$$


The unnormalized lower incomplete gamma function is

$$
\gamma(a,x)=\int_0^x t^{a-1}e^{-t}\,dt
$$


The unnormalized upper incomplete gamma function is

$$
\Gamma(a,x)=\int_x^\infty t^{a-1}e^{-t}\,dt
$$


They satisfy

$$
\Gamma(a,x)=\Gamma(a)-\gamma(a,x)
$$


The C++ and Fortran implementations evaluate the lower incomplete gamma function using:

- a convergent series for $x<a+1$;
- a continued fraction for larger $x$;
- `tgamma` or its language equivalent for the complete gamma function;
- logarithmic gamma evaluation where useful for numerical stability.

The MATLAB implementation follows the same algorithmic structure. The notebook uses SciPy special functions.

### Normalized and unnormalized definitions

Some libraries return the regularized functions

$$
P(a,x)=\frac{\gamma(a,x)}{\Gamma(a)}, \qquad Q(a,x)=\frac{\Gamma(a,x)}{\Gamma(a)}.
$$


The governing equations in this repository use the **unnormalized** incomplete gamma functions. An implementation that substitutes $P$ or $Q$ without multiplying by $\Gamma(a)$ will produce incorrect moments and scale parameters.

### Numerical evaluation strategy

The compiled implementations use a hybrid algorithm:

1. For smaller arguments, a power-series representation is used for the lower incomplete gamma function.
2. For larger arguments, a continued-fraction representation is used for the upper tail.
3. The complementary identity

$$
\Gamma(a,x)=\Gamma(a)-\gamma(a,x)
$$


connects the two forms.
4. Iteration stops when the relative increment is below the local tolerance.
5. Invalid arguments or failure to converge are propagated as calculation errors rather than silently accepted.

This split avoids the severe cancellation and slow convergence that would occur if one numerical representation were used over the entire argument range.

---

## Nonlinear solution for $\widetilde H_1$ and $\widetilde H_2$

### Newton-Raphson system

Let

$$
\mathbf{F}(\mathbf{x})=\begin{bmatrix}F_1(\widetilde H_1,\widetilde H_2)\\F_2(\widetilde H_1,\widetilde H_2)\end{bmatrix}
$$


with

$$
\mathbf{x}=\begin{bmatrix}\widetilde H_1\\\widetilde H_2\end{bmatrix}
$$


At iteration $i$, the Newton correction is obtained from

$$
\mathbf{J}(\mathbf{x}^{(i)})\Delta\mathbf{x}=-\mathbf{F}(\mathbf{x}^{(i)})
$$


followed by

$$
\mathbf{x}^{(i+1)}=\mathbf{x}^{(i)}+\Delta\mathbf{x}
$$


The Jacobian is approximated by centered finite differences. For example,

$$
\frac{\partial F_i}{\partial x_j}\approx\frac{F_i(x_j+\delta)-F_i(x_j-\delta)}{2\delta}
$$


The production constants are:

| Parameter | Value |
|---|---:|
| residual tolerance | $10^{-12}$ |
| Jacobian increment $\delta$ | $10^{-8}$ |
| local gamma tolerance | $10^{-16}$ |
| maximum Newton iterations | 100 |

The C++ and Fortran versions solve the $2\times2$ Newton system directly. The MATLAB function implements the same Newton approach. The Jupyter notebook uses `scipy.optimize.fsolve` with the same residual equations and the same empirical starting values.

### Explicit $2\times2$ Newton correction

With

$$
\mathbf J = \begin{bmatrix} J_{11} & J_{12}\\ J_{21} & J_{22} \end{bmatrix}, \qquad \mathbf F = \begin{bmatrix} F_1\\ F_2 \end{bmatrix},
$$


the determinant is

$$
D=J_{11}J_{22}-J_{12}J_{21}.
$$


The correction solving $\mathbf J\Delta\mathbf x=-\mathbf F$ is

$$
\Delta \widetilde H_1 = \frac{-F_1J_{22}+F_2J_{12}}{D},
$$


$$
\Delta \widetilde H_2 = \frac{-J_{11}F_2+J_{21}F_1}{D}.
$$


A nearly singular Jacobian is rejected when $|D|$ is extremely small. After each correction, non-positive scale parameters are replaced by a small positive tolerance before the next residual evaluation.

### Convergence criterion

The nonlinear solution is accepted only when both governing residuals satisfy

$$
|F_1|<10^{-12}
$$


and

$$
|F_2|<10^{-12}.
$$


Requiring both conditions is important. A small moment residual alone does not guarantee continuity, and a small continuity residual alone does not guarantee the prescribed $H_{rms}$ normalization.

### Approximate initialization formulas

Good initial estimates materially improve convergence of the coupled nonlinear system. Let

$$
x=\widetilde H_{tr}=H_{tr}/H_{rms}
$$


The first initial estimate is

$$
\widetilde H_1^{(0)}=\frac{\tanh\left(6.739139344110821x-0.01265095590917251\right)^{-0.6551633251836707}}{\tanh\left[\sinh\left(0.6947756601426412x+0.7908718490781483\right)\right]^{5.484052848550241}}
$$


The second initial estimate is

$$
\widetilde H_2^{(0)}=1.059259665431797+\frac{0.2059286860468916x}{1+3.865701948059343x^{-3.479682433107255}}
$$


These equations are **numerical initialization approximations**, not additional physical equations in the Battjes-Groenendijk model. They provide starting values only; the final values are obtained by satisfying the continuity and second-moment residuals.

After initialization, non-positive trial values are replaced by a small positive number before the next residual evaluation.

The approximations are smooth regressions of previously solved scale parameters as functions of $\widetilde H_{tr}$. Their purpose is numerical rather than physical:

- they place the initial iterate close to the solution branch;
- they reduce the number of Newton iterations;
- they reduce the risk of crossing into non-positive scale parameters;
- they improve consistency between the compiled and notebook implementations.

The approximations must not be substituted for the final solution. The converged values are always those satisfying the exact continuity and second-moment equations.

---

## Deep-water convergence and overshoot prevention

### Origin of the dimensional inconsistency

The normalized Composite Weibull solution approaches the normalized Rayleigh distribution as $\widetilde H_{tr}$ becomes large. However, dimensional results are obtained by multiplying by the empirical Battjes-Groenendijk $H_{rms}$.

In the deep-water limit, that parameterization tends to

$$
H_{rms}=2.69\sqrt{m_0}
$$


whereas the narrow-band Rayleigh relation used to derive the exact $H_{1/N}/H_{m0}$ limits is

$$
H_{rms}=\sqrt{8m_0}\approx2.828427\sqrt{m_0}
$$


Caires and Van Gent showed that direct dimensionalization can therefore overshoot the Rayleigh prediction and later converge to a value below it. The program prevents this non-physical behavior with two safeguards.

### Safeguard 1: Rayleigh branch

When

$$
\widetilde H_{tr}>2.75
$$


the program bypasses the nonlinear Composite Weibull solution and directly assigns the exact Rayleigh $H_{1/N}/H_{m0}$ values listed above.

In this branch, $H_1$ and $H_2$ are not required and should be treated as not applicable, even if a particular language implementation displays zero or `NaN` placeholders.

### Safeguard 2: dimensional caps

When the Battjes-Groenendijk branch is active, every dimensional mean-high-wave statistic is capped as

$$
H_{1/N}=\min\left(H_{1/N}^{BG},C_NH_{m0}\right)
$$


where $C_N$ is the corresponding exact Rayleigh ratio.

After capping, the software back-calculates the reported normalized values:

$$
\widetilde H_{1/N}=H_{1/N}/H_{rms}
$$


This ensures that dimensional values, normalized values, and diagnostic ratios remain mutually consistent in the report.

### Why a normalized Rayleigh limit is not sufficient

The CWD shape does approach Rayleigh as $\widetilde H_{tr}\rightarrow\infty$. That statement concerns normalized ratios based on the $H_{rms}$ used inside the CWD. It does not by itself guarantee that dimensional values obtained from the empirical finite-bandwidth $H_{rms}$ parameterization approach the exact Rayleigh values based on $H_{m0}$.

The correction is therefore applied in dimensional space. This distinction prevents an apparently converged normalized distribution from producing an inconsistent dimensional upper tail.

### Behaviour at the switch

The threshold

$$
\widetilde H_{tr}=2.75
$$


is an implementation criterion for deep-water convergence, not a physical breaker index. Immediately above the threshold, the direct Rayleigh constants are used. At or below the threshold, the CWD is solved and each final dimensional statistic is capped individually.

Because the caps are applied statistic by statistic, the selected CWD tail may coincide with the Rayleigh limit for one $N$ while remaining below it for another. The final reported normalized values are recalculated from the capped dimensional values so that all output columns remain internally consistent.

---

## Calculation sequence

For each valid input set, the software performs the following operations:

1. Validate that $H_{m0}$, $d$, and $M$ are finite and positive.
2. Calculate $m_0=(H_{m0}/4)^2$.
3. Calculate $H_{rms}$ from the Battjes-Groenendijk empirical relation.
4. Calculate $\tan\alpha=1/M$.
5. Calculate $H_{tr}=(0.35+5.8\tan\alpha)d$.
6. Calculate $\widetilde H_{tr}=H_{tr}/H_{rms}$.
7. Select the Rayleigh branch when $\widetilde H_{tr}>2.75$.
8. Otherwise, initialize and solve for $\widetilde H_1$ and $\widetilde H_2$.
9. Calculate $\widetilde H_N$ and $\widetilde H_{1/N}$ for each requested $N$.
10. Convert the normalized results to metres.
11. Apply the Rayleigh dimensional caps.
12. Recalculate normalized values from the final capped dimensional results.
13. Calculate diagnostic ratios relative to $H_{1/3}$.
14. Display the results and write `report.txt` where supported.

### Reference pseudocode

```text
read Hm0, d, M
validate Hm0 > 0, d > 0, M > 0 and all values finite

m0        = (Hm0 / 4)^2
sqrt_m0   = sqrt(m0)
Hrms      = (2.69 + 3.24 * sqrt_m0 / d) * sqrt_m0
tan_alpha = 1 / M
Htr       = (0.35 + 5.8 * tan_alpha) * d
Htr_tilde = Htr / Hrms

if Htr_tilde > 2.75:
    distribution = Rayleigh
    assign exact Rayleigh H1/N = C_N * Hm0
else:
    distribution = B&G
    evaluate empirical initial guesses
    solve continuity and second-moment residuals
    for N in {3, 10, 50, 100, 250, 1000}:
        determine the branch containing H_N
        integrate the corresponding upper tail
        dimensionalize with Hrms
        cap the result at C_N * Hm0
        recalculate the normalized final value

calculate ratios H1/N / H1/3
format the report
write report.txt
```

### Computational complexity

The calculation is small. Each case requires only a two-variable nonlinear solve and a limited number of incomplete-gamma evaluations. Runtime is normally negligible compared with spectral wave propagation, CFD, Boussinesq, or phase-resolving models. The primary implementation concern is numerical correctness and cross-language consistency rather than computational cost.

---

## Inputs

### $H_{m0}$: spectral significant wave height

Enter the local spectral significant wave height in metres. The value should correspond to the same location and water level as the input depth.

### $d$: local water depth

Enter the positive local still-water depth in metres. The model does not add tides, storm surge, wave setup, or sea-level rise automatically.

### $M$: slope denominator

Enter the denominator of the local foreshore slope $1:M$.

Examples:

| Foreshore slope | Input $M$ | $\tan\alpha$ |
|---|---:|---:|
| $1:20$ | 20 | 0.050000 |
| $1:50$ | 50 | 0.020000 |
| $1:100$ | 100 | 0.010000 |
| $1:250$ | 250 | 0.004000 |

Use a representative local or averaged foreshore slope consistent with the point-model assumptions. Do not enter a slope angle in degrees.

---

## Outputs

### Distribution type

- `B&G`: the Composite Weibull equations were solved.
- `Rayleigh`: the normalized transition exceeded 2.75 and the direct Rayleigh branch was used.

### Scale parameters

$H_1$ and $H_2$ are mathematical scale parameters of the two Weibull branches. They are not breaker heights and should not be interpreted as independent physical wave statistics.

### Mean-high-wave statistics

The program reports:

- $H_{1/3}$: mean of the highest one-third;
- $H_{1/10}$: mean of the highest one-tenth;
- $H_{1/50}$: mean of the highest one-fiftieth;
- $H_{1/100}$: mean of the highest one-hundredth;
- $H_{1/250}$: mean of the highest one-two-hundred-and-fiftieth;
- $H_{1/1000}$: mean of the highest one-thousandth.

Both normalized values $H/H_{rms}$ and dimensional values in metres are reported.

### Diagnostic ratios

The report includes

$$
H_{1/10}/H_{1/3}
$$


$$
H_{1/50}/H_{1/3}
$$


$$
H_{1/100}/H_{1/3}
$$


$$
H_{1/250}/H_{1/3}
$$


$$
H_{1/1000}/H_{1/3}
$$


These ratios describe the upper-tail shape independently of the absolute wave-height scale.

---

## Worked example

For

- $H_{m0}=2.5$ m;
- $d=5.0$ m;
- slope $1:100$;

then

$$
m_0=(2.5/4)^2=0.390625\ \text{m}^2
$$


$$
H_{rms}=1.9344\ \text{m}
$$


$$
H_{tr}=(0.35+5.8/100)\times5=2.0400\ \text{m}
$$


$$
\widetilde H_{tr}=2.0400/1.9344=1.0546
$$


Because $\widetilde H_{tr}<2.75$, the `B&G` branch is used. The C++ implementation gives approximately:

| Quantity | Result |
|---|---:|
| $\widetilde H_1$ | 1.1567 |
| $\widetilde H_2$ | 1.1102 |
| $H_{1/3}$ | 2.5027 m |
| $H_{1/10}$ | 2.9701 m |
| $H_{1/50}$ | 3.3296 m |
| $H_{1/100}$ | 3.4567 m |
| $H_{1/250}$ | 3.6077 m |
| $H_{1/1000}$ | 3.8086 m |
| $H_{1/1000}/H_{1/3}$ | 1.5218 |

Small differences in the last displayed digit can occur between language implementations because of formatting and special-function evaluation.

### Rayleigh-branch example

For

- $H_{m0}=2.5$ m;
- $d=20.0$ m;
- slope $1:100$;

the intermediate values are

$$
m_0=\left(\frac{2.5}{4}\right)^2=0.390625\ \mathrm{m}^2,
$$


$$
H_{rms} = \left(2.69+3.24\frac{0.625}{20}\right)0.625 = 1.74453125\ \mathrm{m},
$$


$$
H_{tr} = \left(0.35+\frac{5.8}{100}\right)20 = 8.16\ \mathrm{m},
$$


and

$$
\widetilde H_{tr} = \frac{8.16}{1.74453125} \approx4.6772.
$$


Because $\widetilde H_{tr}>2.75$, the direct Rayleigh branch is selected. The corresponding dimensional statistics are:

| Quantity | Result |
|---|---:|
| $H_{1/3}$ | 2.502689 m |
| $H_{1/10}$ | 3.181836 m |
| $H_{1/50}$ | 3.900283 m |
| $H_{1/100}$ | 4.170583 m |
| $H_{1/250}$ | 4.502543 m |
| $H_{1/1000}$ | 4.962089 m |

This example also demonstrates that the direct Rayleigh branch uses exact ratios to $H_{m0}$ rather than dimensionalizing the empirical Battjes-Groenendijk $H_{rms}$.

---

## Sensitivity and dimensional interpretation

### Dependence on $H_{m0}/d$

The degree of saturation is proportional to $H_{m0}/d$. Increasing local wave height at fixed depth, or reducing depth at fixed wave height, moves the distribution farther from Rayleigh. The largest changes generally occur in the high-N statistics because breaking primarily suppresses the upper tail.

### Dependence on slope

At fixed $H_{m0}$ and $d$, decreasing $M$ means a steeper foreshore. The transition formula then gives a larger $H_{tr}$ and a larger $\widetilde H_{tr}$. In the point-model parameterization, this leaves more of the distribution in the lower Rayleigh-shaped branch. The slope dependence should be interpreted as an empirical representation of spatial adaptation, not as a universal statement that steeper seabeds always produce larger physical waves.

### Dimensional scaling

For fixed dimensionless conditions, the normalized quantities remain unchanged and dimensional wave heights scale with $H_{rms}$. However, changing $H_{m0}$ at fixed depth also changes the degree of saturation, so practical scaling is not exactly linear unless the ratios $H_{m0}/d$ and the slope are preserved.

### Upper-tail sensitivity

As $N$ increases, $H_{1/N}$ samples progressively rarer waves. These values are more sensitive to:

- the selected distribution;
- numerical precision in incomplete gamma functions;
- the deep-water cap;
- model-domain violations;
- measurement sampling uncertainty;
- the number of waves available for empirical validation.

Consequently, $H_{1/1000}$ should not be treated as equally well constrained as $H_{1/3}$ merely because both are produced by the same deterministic formula.

---

## Model applicability

### Original calibration domain

The underlying laboratory database included plane shallow foreshores with slopes approximately between

$$
1:20\ \text{and}\ 1:250
$$


The 1998 validation report identified the main validated conditions as:

- degree of saturation approximately $0<\sqrt{m_0}/d<0.15$;
- equivalently, approximately $0<H_{m0}/d<0.6$ because $H_{m0}=4\sqrt{m_0}$;
- straight, parallel depth contours;
- simple, gradually varying foreshore geometry;
- local wave conditions that have had sufficient propagation distance to adapt to the depth.

The source database included slopes $1:20$, $1:30$, $1:50$, $1:100$, and $1:250$ and predominantly JONSWAP-type laboratory spectra.

### Point-model assumption

The distribution is assumed to depend primarily on the local $m_0$, $d$, and $\tan\alpha$, rather than on the full transformation history. This assumption is most defensible on a long, simple, monotonic foreshore.

### Steep foreshores

For slopes steeper than approximately $1:20$, the transitional-height relation becomes less well supported and the spatial lag of breaking becomes increasingly important. A Rayleigh estimate may provide a conservative alternative, but project-specific data or a more suitable model should be preferred.

### Very mild and flat bottoms

A flat bottom is not equivalent to a $1:250$ foreshore. On a constant-depth bed:

- waves are not continuously forced to break by a shoreward depth reduction;
- local wind input may be important;
- the wave-height distribution may remain closer to Rayleigh than the B&G model predicts;
- applying $M=250$ can underestimate high-wave statistics.

Caires and Van Gent found mean underestimation of approximately 7% to 15% for high-wave quantities on shallow flat bottoms. For such cases, a Rayleigh distribution is generally the safer default unless a validated constant-depth model or site-specific measurements are available.

### Bars, troughs, and rapidly varying bathymetry

The point model does not explicitly retain wave-field memory. On barred profiles, the distribution shoreward of a bar may gradually reform toward Rayleigh as waves enter deeper water, then depart from Rayleigh again if depth-induced breaking resumes. The simple local formula cannot fully represent this non-local evolution.

Use caution for:

- offshore bars and troughs;
- abrupt dredged channels;
- reef edges;
- submerged structures;
- rapidly changing slopes;
- strongly two-dimensional bathymetry;
- locations dominated by refraction, diffraction, or wave-current interaction.

### Other limitations

The calculator does not explicitly include:

- wave period or spectral shape as an input;
- directional spreading;
- wave-current interaction;
- reflected waves;
- infragravity-wave separation;
- local wind growth;
- wave grouping;
- nonstationarity;
- individual-wave maximum statistics based on storm duration;
- uncertainty propagation.

The model should not be used outside its calibration domain without engineering judgement, sensitivity analysis, and comparison against alternative models or measurements.

### Calibration evidence and interpretation of the stated range

The calibration range is not a rigid mathematical boundary at which the equations suddenly become invalid. It identifies the experimental region in which the parameters were fitted and checked. Predictions near the edge of the range should receive more scrutiny than predictions near the central portion of the database.

The original data represent controlled laboratory conditions. Field applications may include directional spreading, currents, irregular bathymetry, wind input, reflection, and mixed sea-swell systems that were not independently parameterized. Agreement of the input ratios with the laboratory range is therefore necessary but not sufficient evidence of applicability.

### Recommended engineering checks

Before adopting a result, verify:

1. that $H_{m0}$ and $d$ refer to the same location and water level;
2. that the slope represents the adapting foreshore rather than an isolated local cell;
3. that the bathymetry is sufficiently monotonic for a point-model interpretation;
4. that the selected statistic is the one required by the downstream design method;
5. that Rayleigh, CWD, and any project-specific measured or modelled distributions have been compared;
6. that the result is not being interpreted as a deterministic maximum;
7. that sensitivity to water level, wave height, slope, and the 2.75 switch has been examined.

---

## Engineering interpretation

Accurate upper-tail wave statistics are relevant to:

- armour stability and wave loading;
- revetment and dike design;
- run-up and overtopping assessments;
- crest-level studies;
- reliability analysis;
- interpretation of shallow-water physical-model data.

The calculator supplies wave-height statistics only. It does not replace the hydraulic formula, numerical model, physical model, design standard, or partial-safety framework required for the final engineering application.

For overtopping or structural calculations, verify which wave-height definition the selected design method requires. $H_{m0}$, $H_{1/3}$, an exceedance quantile, and the mean of the highest fraction are not interchangeable in shallow water.

---

## Project files

| File | Purpose |
|---|---|
| `shallow-water-waves_cli.cpp` | portable C++ command-line implementation |
| `shallow-water-waves_gui.cpp` | native Windows C++ graphical interface |
| `shallow-water-waves_cli.f90` | Fortran command-line implementation |
| `shallow_water_waves.m` | MATLAB function returning a results structure |
| `shallow-water-waves.ipynb` | Python/Jupyter implementation, tables, and plots |
| `android.apk` | Android Studio project for the Weibull Waves app |
| `README.md` | theory, equations, limitations, build instructions, and usage |

All implementations use the same physical parameterization, target statistics, Rayleigh limits, and nonlinear governing equations.

---

## C++ command-line program

### Compile with GCC or MinGW-w64

Portable build:

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic \
    -Wconversion -Wsign-conversion \
    shallow-water-waves_cli.cpp \
    -o shallow-water-waves_cli
```

Optimized static Windows build, where the required static libraries are available:

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
    -Wconversion -Wsign-conversion \
    -static -static-libgcc -static-libstdc++ \
    shallow-water-waves_cli.cpp \
    -o shallow-water-waves_cli.exe
```

### Run with command-line arguments

Linux or macOS:

```bash
./shallow-water-waves_cli 2.5 5.0 100
```

Windows:

```bat
shallow-water-waves_cli.exe 2.5 5.0 100
```

The argument order is:

```text
Hm0  depth  slope_denominator
```

### Run interactively

```bash
./shallow-water-waves_cli
```

When the three arguments are omitted, the program prompts for each input.

### Output

The formatted report is printed to the terminal and written to

```text
report.txt
```

in the current working directory.

---

## Native Windows C++ GUI

### Compile with MinGW-w64

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
    -Wconversion -Wsign-conversion -municode \
    shallow-water-waves_gui.cpp \
    -o shallow-water-waves_gui.exe \
    -mwindows -static -static-libgcc -static-libstdc++
```

### Run

```bat
shallow-water-waves_gui.exe
```

Enter:

- `Hm0 (m)`;
- `d (m)`;
- `Beach slope m` for a slope $1:m$;

then press **Compute**. The report is displayed in the window and written to `report.txt`.

---

## Weibull Waves Android app

The Android Studio project is in the `android` directory. The app is written in Kotlin with Jetpack Compose and supports Android 6.0 (API 23) and newer.

It accepts the same three local inputs used by the other implementations: $H_{m0}$, water depth $d$, and slope denominator $M$ for a slope $1:M$. It reports the selected distribution branch, calculated parameters, dimensional and normalized mean-high-wave statistics, diagnostic ratios, and a copyable or shareable engineering report. The `Weibull Theory` link opens an in-app explanation of the equations, numerical method, interpretation, and model limits.

### Install `weibull.apk`

1. Copy `weibull.apk` to the phone or download it directly to the phone.
2. If requested, allow the browser or file manager opening the APK to install unknown apps.
3. Open `weibull.apk`, approve the installation, and launch `Weibull Waves`.
4. The install-from-source permission may be disabled again after installation.

If Android reports that the app is not installed while an older copy is present, the installed copy may use a different signing certificate or a newer version code. Uninstalling the old copy normally resolves a certificate conflict, but it also removes that app's locally stored data.

With Android Platform Tools and USB debugging enabled, installation can also be performed with:

```bash
adb install -r weibull.apk
```

---

## Fortran command-line program

### Compile with GNU Fortran

```bash
gfortran -O3 -march=native -std=f2018 -Wall -Wextra -pedantic \
    -Wconversion -static -fno-underscoring \
    shallow-water-waves_cli.f90 \
    -o shallow-water-waves_cli_f
```

When static linking is unavailable, remove `-static`.

### Run

```bash
./shallow-water-waves_cli_f 2.5 5.0 100
```

or interactively:

```bash
./shallow-water-waves_cli_f
```

The program prints the report and writes `report.txt`.

---

## MATLAB function

Place `shallow_water_waves.m` on the MATLAB path or in the current working directory.

### Return the result structure

```matlab
results = shallow_water_waves(2.5, 5.0, 100);
```

### Print the report without assigning an output

```matlab
shallow_water_waves(2.5, 5.0, 100);
```

The function validates all three inputs, returns the full results structure, and writes `report.txt`.

Key structure fields include:

```text
Hm0
d
slopeM
distribution_type
m0
Hrms
tanAlpha
Htr_dim
Htr_tilde
H1_Hrms
H2_Hrms
H1_3_Hrms ... H1_1000_Hrms
H1_dim
H2_dim
H1_3_dim ... H1_1000_dim
ratio_1_10_div_1_3 ... ratio_1_1000_div_1_3
```

---

## Jupyter notebook

### Create a virtual environment

Windows:

```bat
py -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install numpy scipy pandas matplotlib notebook
jupyter notebook shallow-water-waves.ipynb
```

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install numpy scipy pandas matplotlib notebook
jupyter notebook shallow-water-waves.ipynb
```

The notebook contains:

- theoretical background;
- incomplete-gamma functions;
- nonlinear solution with `fsolve`;
- the empirical initialization formulas;
- the main analysis workflow;
- formatted result tables;
- graphical analysis of the distribution and characteristic wave heights.

Run the cells sequentially and modify the input cell before executing the final analysis cells.

---

## Numerical verification

Recommended cross-language verification steps:

1. Use the same $H_{m0}$, $d$, and $M$ in every implementation.
2. Confirm the same `distribution_type`.
3. Compare $m_0$, $H_{rms}$, $H_{tr}$, and $\widetilde H_{tr}$.
4. In the `B&G` branch, compare $\widetilde H_1$ and $\widetilde H_2$.
5. Compare all final dimensional $H_{1/N}$ values after capping.
6. Compare diagnostic ratios using the final capped values.
7. Confirm both nonlinear residuals are below the specified tolerance.

Minor floating-point differences are expected, but engineering outputs should agree to the displayed precision.

### Residual and identity checks

For a converged B&G solution, the following identities should be checked.

Continuity:

$$
\left(\frac{\widetilde H_{tr}}{\widetilde H_1}\right)^{k_1} - \left(\frac{\widetilde H_{tr}}{\widetilde H_2}\right)^{k_2} \approx0.
$$


Second moment:

$$
\widetilde H_1^2 \gamma\left(1+\frac{2}{k_1},x_1\right) + \widetilde H_2^2 \Gamma\left(1+\frac{2}{k_2},x_2\right) \approx1.
$$


Monotonicity:

$$
H_{1/3}<H_{1/10}<H_{1/50}<H_{1/100}<H_{1/250}<H_{1/1000}.
$$


Rayleigh cap:

$$
H_{1/N}\leq C_NH_{m0}.
$$


Dimensional-normalized consistency:

$$
H_{1/N} = \widetilde H_{1/N}H_{rms}.
$$


The same checks should be applied automatically in regression tests wherever practical.

### Reproducibility tolerance

For identical inputs and the same distribution branch, cross-language values should generally agree much more closely than the engineering precision of the final result. Differences in the final few binary or printed digits may arise from:

- the standard-library gamma implementation;
- compiler floating-point optimization;
- `fsolve` versus the explicit Newton iteration;
- decimal output formatting;
- intermediate rounding in notebook tables.

A branch mismatch or a difference that changes engineering conclusions is not a normal floating-point effect and should be investigated.

Useful regression cases should include:

- a strongly depth-limited case with small $\widetilde H_{tr}$;
- a typical shallow-foreshore case;
- a case close to $\widetilde H_{tr}=2.75$;
- a Rayleigh-branch case above the threshold;
- invalid zero, negative, non-numeric, and non-finite inputs.

---

## Common interpretation errors

**Confusing $H_{m0}$ and $H_{1/3}$.** They are approximately equal in deep water but not generally identical on a shallow foreshore.

**Confusing $H_N$ and $H_{1/N}$.** One is a threshold; the other is a conditional mean above that threshold.

**Entering the tangent instead of the slope denominator.** For $1:100$, enter `100`, not `0.01`.

**Using the model on a flat bottom by entering an arbitrarily large $M$.** The flat-bottom process is outside the original formulation and can produce non-conservative high-wave estimates.

**Using offshore $H_{m0}$ with local depth.** The input $H_{m0}$ must be local to the point being analysed.

**Interpreting $H_1$ or $H_2$ as physical breaker heights.** They are Weibull scale parameters.

**Treating $H_{1/1000}$ as the maximum wave in a storm.** Storm maximum statistics also depend on the number of waves and the duration of the sea state.

---

## Core references

1. Battjes, J. A. (1972). *Set-up due to irregular waves*. Proceedings of the 13th International Conference on Coastal Engineering.
2. Battjes, J. A., and Janssen, J. P. F. M. (1978). *Energy loss and set-up due to breaking of random waves*. Proceedings of the 16th International Conference on Coastal Engineering, 569–587.
3. Battjes, J. A., and Groenendijk, H. W. (2000). *Wave height distributions on shallow foreshores*. Coastal Engineering, 40, 161–182.
4. Caires, S., and Van Gent, M. R. A. (2012). *Wave height distribution in constant and finite depths*. Proceedings of the 33rd International Conference on Coastal Engineering.
5. Goda, Y. (1975). *Irregular wave deformation in the surf zone*. Coastal Engineering in Japan, 18, 13–26.
6. Goda, Y. (1979). *A review on statistical interpretation of wave data*. Report of the Port and Harbour Research Institute, 18(1), 5–32.
7. Groenendijk, H. W. (1998). *Shallow foreshore wave height statistics*. MSc thesis, Delft University of Technology.
8. Groenendijk, H. W., and Van Gent, M. R. A. (1998). *Shallow foreshore wave height statistics: A predictive model for the probability of exceedance of wave heights*. WL | Delft Hydraulics Report H3351.
9. Hasselmann, K., Barnett, T. P., Bouws, E., Carlson, H., Cartwright, D. E., Enke, K., Ewing, J. A., Gienapp, H., Hasselmann, D. E., Kruseman, P., Meerburg, A., Müller, P., Olbers, D. J., Richter, K., Sell, W., and Walden, H. (1973). *Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP)*.
10. Klopman, G. (1996). *Extreme wave heights in shallow water*. WL | Delft Hydraulics Report H2486.
11. Longuet-Higgins, M. S. (1952). *On the statistical distribution of the heights of sea waves*. Journal of Marine Research, 11(3), 245–266.
12. Longuet-Higgins, M. S. (1980). *On the distribution of the heights of sea waves: Some effects of nonlinearity and finite bandwidth*. Journal of Geophysical Research, 85(C3), 1519–1523.
13. Naess, A. (1985). *On the distribution of crest-to-trough wave heights*. Ocean Engineering, 12(3), 221–234.
14. Ochi, M. K. (1998). *Ocean Waves: The Stochastic Approach*. Cambridge University Press.
15. Tayfun, M. A. (1990). *Distribution of large wave heights*. Journal of Waterway, Port, Coastal, and Ocean Engineering, 116(6), 686–707.
16. Thornton, E. B., and Guza, R. T. (1982). *Energy saturation and phase speeds measured on a natural beach*. Journal of Geophysical Research, 87(C12), 9499–9508.
17. Thornton, E. B., and Guza, R. T. (1983). *Transformation of wave height distribution*. Journal of Geophysical Research, 88(C10), 5925–5938.
