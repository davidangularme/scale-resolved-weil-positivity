# Configuration Space Temporality v33: Corrections and a Resolution Law

Oct 8, 2026 · @Frédéric David Blum

Frédéric David Blum, Catalyst AI Research, Haifa. Version 33 of the Zenodo record 10.5281/zenodo.18859602. Computations and drafting assisted by Claude (Anthropic).

## Summary

Version 33 withdraws every claim in this archive that the Riemann Hypothesis has been proved, reduced to one assumption, or approached through the Bost–Connes Hamiltonian, the Fisher–Rao action or the rigidity matrix. A re-audit with every number recomputed found four errors; the nine withdrawn results and the eight that stand are listed in Section 2. The physical framework of Configuration Space Temporality (CST) — time as the cost of transitions in configuration space, with a bounded transition rate — is kept, and Section 3 states it in one page with its one validated measurement and its untested predictions.

The new content of this version is a single law, stated in Section 4 and proved in Section 5 as a theorem of time–band limiting. CST's third axiom (a bounded transition rate) is, in quantum mechanics, the family of quantum speed limits (Mandelstam–Tamm, Margolus–Levitin): a system with energy scale ΔE traverses at most a number of order ΔEΔτ/h of mutually orthogonal states in a time Δτ. The archive's cascade number κ = ΔEΔτ/(ħ ln 2) is a count of that kind, expressed in bits. The law says what the trajectory contains beyond the count: for a state whose energy distribution is flat on a band of width Ω\_E, the number of modes of the trajectory distinguishable at level ε is not a count but

```latex
d_\varepsilon(\Delta\tau)=\frac{\Omega_E\,\Delta\tau}{h}+\frac{1}{\pi^{2}}\,\ln c\,\ln\frac{1-\varepsilon}{\varepsilon}+o(\ln c),\qquad c=\frac{\Omega_E\,\Delta\tau}{4\hbar},
```

so that the horizon of a physical process is soft: distinguishability falls past the Landau count Ω\_EΔτ/h at a universal rate fixed by the Landau–Widom theorem, not at a wall, and is cut off super-exponentially beyond twice the Slepian parameter. The law has no free parameter; it is a theorem once the energy distribution of the state is given; for non-flat distributions the count and the plunge are bounded on both sides (Section 5) and the leading term is set by the level set of the distribution; and it is falsifiable by a time-resolved measurement that the archive's existing qubit data do not provide (Section 8). It was found by transplanting to CST the structure of the companion paper *Scale-Resolved Weil Positivity* (Zenodo 10.5281/zenodo.23229271), where the same time–band limiting fixes which zeros of ζ the primes up to λ² determine; Section 7 states what transfers from that arithmetic instance and what does not.

What is new is stated in Section 8 without inflation: in physics, the identification of the cascade number with the Landau count of the trajectory and the soft horizon beyond the quantum speed limits; in mathematics, the identification of the Gram operator of a trajectory with Slepian's time–band limiting operator, which brings the Landau–Widom plunge and an explicit super-exponential cutoff to the effective dimension of a trajectory, a quantity otherwise studied in equilibration theory and as Krylov spread. The law does not bear on gravity, dark matter or the Riemann Hypothesis.

## 2. Corrections

Of the 17 RH-related results in versions 25–32, 3 hold as stated, 3 hold but are trivial or known, 1 is falsified numerically, 1 stays open (Li positivity, equivalent to RH) and 9 are withdrawn. The four errors behind the withdrawals:

**Error 1 — the Fisher pole identity is generic.** The archive's "quantization condition" |I(s)|·|s − ρ|² = 1, with I = (log ζ)″, holds at every simple zero of every analytic function, on or off the critical line (checked on (s − a)(s − b)eˢ with a = 0.8 + 20i: the product is 1.0000000). It states that zeros are simple; it says nothing about where they are. For complex s the weights n⁻ˢ/ζ(s) are not a probability distribution, so I(s) is a Fisher information only for real s > 1.

**Error 2 — the rigidity matrix cannot confine.** Every row of R sums to zero, so the constant vector is an exact null eigenvector and the confinement constant α = μ₁/μ₀ is undefined (μ₀ = −1.5 × 10⁻⁶¹ at 60 digits). The sign is reversed as well: for the log-gas energy the collinear configuration is a transverse maximum (Earnshaw), not a minimum. R is built from ordinates already on the line, so no property of R can be equivalent to RH. The Fisher–Rao action diverges on every path through a zero (1/|s − ρ| is the 3D Coulomb law, not the 2D one).

**Error 3 — exponential dominance fails, and the effect is not arithmetic.** The v31 table (log₁₀ μ₀ = −97.58, α = 0.1427 at N = 45, λ = 3) is reproduced to two decimals, but μ₁/μ₀ grows only from 10⁶·⁴² at N = 45 to 10⁷·⁸¹ at N = 200, so α(N) ≈ c/N → 0. GUE, Poisson and evenly spaced control points of the same density give α = 0.0609, 0.0611, 0.0608 against 0.0608 for the zeros: the tiny eigenvalues come from the cosine basis on an interval shorter than 2π, not from the zeros. The matrix M = 2ΦΦᵀ is a sum over the zeros themselves and is positive by construction; it presupposes RH.

**Error 4 — the operator claims.** The Bost–Connes Hamiltonian has spectrum {log n}, which is why Tr e⁻ᵗᴴ = ζ(t); the regularised determinant of Assumption A vanishes at s = ½ ± i√(log n − ¼), the first at height 0.666, not 14.13. The resolvent's poles are at log n, all real; zeros of a Laplace transform are not spectral values. H\_mod is self-adjoint, so Assumption B's "real part ½" does not apply to it, and the reflection s ↔ 1 − s comes from the archimedean factor, absent from the Bost–Connes partition function. The CMP v25–v27 "complete proof chain" is withdrawn in full; the Krein–Rutman step was already retracted in v29. If submission CIMP-D-26-00460 is still with the journal, the editor is to be told.

| Result | v32 | v33 |
| --- | --- | --- |
| Fisher pole quantization | Proved | Holds, trivial |
| Coulomb reduction of the Fisher–Rao action | Exact | Withdrawn |
| Rigidity matrix R, constant α | Open | Withdrawn |
| α∞ > 0 ⇔ RH (Thm 4.3, Conj 7.2) | Conjecture | Withdrawn |
| Generalized Stieltjes principle | Proved | Holds, vacuous |
| Polyak–Łojasiewicz inequality | Proved | Holds (standard) |
| Exponential isolation table, N ≤ 45 | Verified | Holds |
| Spectral gap Δ(N) > 0 | Verified | Holds (also for random points) |
| Exponential dominance (Conj ED) | Open | Falsified |
| Weil matrix M as "the Weil form" | Assumed | Withdrawn, replaced by the prime-side form |
| Assumption A (det\_reg = C·ξ) | Open | Withdrawn |
| Assumption B (spectral rigidity via KMS) | Proved | Withdrawn |
| Zeros as resolvent poles of H\_BC | Proved | Withdrawn |
| CMP v25–v27 proof chain | Partly retracted | Withdrawn |
| Fisher–Bakry–Émery bridge | Open | Withdrawn |
| Li criterion bridge | Known | Holds, known |
| Li positivity for n > 34 | Open | Open (⇔ RH) |

Files marked superseded: CMP\_Paper\_FDBC\_v25–v27; FDBC\_v2\_CMP\_Submission; FDBC\_v2\_Complete\_Proof\_Chain; FDBC\_v2\_Proof\_Chain\_CORRECTED; FDBC\_v2\_Assumption\_B\_PROVED (both); FDBC\_v2\_Spectral\_Gap\_PROVED; FDBC\_v2\_No\_Gap\_Final; Information\_Geometric\_Confinement\_Riemann\_Zeros\_Blum\_2026; blum-claude-weil-minimiser-2026 (its table stands, its conclusions do not). The full audit with every computation is the document *CST/FDBC Archive v33 — Corrections and a Non-Circular Weil-Form Test* in this record; the positive mathematics that came out of it is the companion paper (Section 7).

## 3. CST in one page

**Axioms.** (1) Configuration uniqueness: each distinct configuration of a physical system is one point of a smooth manifold, the configuration space. (2) Energy-driven transition: a change of configuration costs energy, at a rate governed by a transition potential Φ on pairs of configurations. (3) Bounded rate: no transition rate exceeds a maximum fixed by the energy available. Time is the accumulated cost of transitions; it is not a background parameter.

**Objects.** The tensor of transition G, a metric on configuration space (the Bures metric on quantum configurations, with the Umegaki relative entropy as potential); the transport distance d\_CST, which satisfies the Carlen–Maas axioms of a quantum Wasserstein distance (v18, T2); the Witten Laplacian H\_W = −Δ + |∇W|² − ΔW, whose thermodynamic limit on the configuration torus is a Lindblad generator (v18, T1); and the cascade number κ = ΔE·Δτ/(ħ ln 2), the archive's measure of the irreversibility of a process.

**What is validated.** One measurement: on 768 qubit pairs of six IBM Eagle r3 processors, the spectral gap λ₁ − λ₀ of H\_W(transmon) correlates with two-qubit gate fidelity at r = 0.365 after the bare calibration parameters (f₀₁, α, E\_J/E\_C) are projected out; the signal survived the methodological audit that removed the earlier universality claim (CV = 0.014). This is a diagnostic correlation across calibration snapshots, not a time-resolved measurement; it does not test the law of Section 4 and is not claimed to.

**What is not validated.** The nine predictions of the February 2026 roadmap, among them P2 (cascade criticality: a sharp reversible/irreversible transition at κ = 1), P6 (no dark-matter particles) and P9 (gravitational waves without a graviton pole), remain untested. Section 4 replaces P2: the law gives a smooth plunge of calculable width, not a threshold at κ = 1.

**What this version removes from CST.** Everything that tied CST to the zeros of ζ through the Bost–Connes system (Section 2). CST is a physical framework about transitions and their cost; it has no arithmetic content, and the archive's attempt to give it one was wrong. The relation that survives runs the other way: the harmonic analysis developed for the Weil form (Section 7) applies to any process with a bounded transition rate, and that is where the law comes from.

## 3′. Foundations: where CST stands among timeless formulations

CST's founding thesis — there is no time, only states, and what we call time is a reading of transitions between states — is not an isolated intuition. It is the thesis of three established programmes, each with a precise mathematical procedure, and CST should be placed on them explicitly; doing so costs it nothing and gives it the ground the archive lacked.

**Jacobi geometry and Barbour's relational dynamics.** For a system with configuration q and potential V, Jacobi's principle replaces Newton's equations with geodesics of the metric ds² = 2(E − V)·|dq|²\_m on configuration space, with no time parameter; the arc length s = ∫√(2(E − V))|dq|\_m = ∫ p·dq is the Maupertuis action, and Newtonian time is recovered as the parameter in which the motion is geodesic (the ephemeris time). Barbour and Bertotti \[BB82\] built Newtonian mechanics on this with no absolute time or space (best matching), Barbour, Foster and Ó Murchadha \[BFO02\] showed that general relativity itself is such a theory (the Baierlein–Sharp–Wheeler action), and shape dynamics makes the scale relational as well. CST's tensor of transition G is a Jacobi metric: the cost of a transition is a length, and CST's time is the arc length. In the classical limit the cascade number is the Maupertuis action in bits, κ\_cl = S₀/(ħ ln 2), which is the identification Section 4 makes on the quantum side (Landau count = action budget in units of h).

**Wheeler–DeWitt and Page–Wootters.** In canonical quantum gravity the state of the universe satisfies HΨ = 0: nothing evolves. Page and Wootters \[PW83\] showed that time nonetheless appears as a correlation between a subsystem used as a clock and the rest, conditional on the clock reading; the mechanism has been realised experimentally with entangled photons \[MGV14\]. This is the quantum form of "states, not time", and it is where CST's statement that the Schrödinger equation is the optimal-transport equation on configuration space has to be made precise: the transport is in the clock correlation, not in an external parameter.

**Thermal time.** Connes and Rovelli \[CR94\] showed that a state on an algebra of observables singles out a flow, the modular (Tomita–Takesaki) flow, and proposed that physical time is this flow: time is a property of the state, with temperature as the ratio between thermal time and geometric time. The KMS condition and the Bost–Connes system that the archive used to reach the zeros of ζ are this machinery; this is its correct use in CST, as the origin of the time parameter from the state, and it has nothing to do with the zeros.

**The point that decides what a law can be.** All three formulations are empirically equivalent to their timed counterparts by construction: reparametrisation invariance is the statement that the parameter carries no physics. A "unified law of states" that reproduces known physics is therefore a reformulation, however deep, not a new law. A new law requires a place where the state-based description predicts something the timed one does not. One such place exists in print: Barbour, Koslowski and Mercati \[BKM14\] showed that in the E = 0 Newtonian N-body problem the shape complexity C\_S = ℓ\_rms/ℓ\_mhl has a unique minimum on every solution (the Janus point) and grows on both sides of it, so the arrow of time is a property of the geometry of states, not of a special initial condition.

**A tested question: is the cascade number tied to shape complexity?** If CST's clock is the Jacobi length s, the natural conjecture is a universal relation between C\_S and s. We tested it (janus.py): random E = 0 three-body solutions, integrated from the Janus region in both directions to s ≈ 120–140, with C\_S and s computed along the way. The Janus structure is reproduced (C\_S minima of 0.20–0.26, growth on both branches). The relation C\_S ∝ s^α is not universal: the fitted exponents are 0.92, 0.96, 0.92 on the forward branches and 1.90, 0.81, 1.23 on the backward ones, and the reason is plain — the generic end state is a bound pair plus an escaper, the Jacobi length is then dominated by the pair's kinetic energy (s ≈ 2Kt with K constant), and the exponent depends on which asymptotic state the solution reaches. There is no law of the form "complexity grows as a fixed power of the CST clock"; what the computation shows is that the CST clock is the ephemeris time of the bound pair, which is Barbour's conclusion. The negative result is recorded so that it is not tried again as a positive one.

## 4. The Resolution Law

**Setting.** A quantum system with Hamiltonian H is prepared in a state ψ whose energy distribution μ (the spectral measure of H in ψ) is supported in an interval of width Ω\_E. The trajectory is ψ(t) = e^{−iHt/ħ}ψ for 0 ≤ t ≤ Δτ. The question the law answers is: how many distinguishable states does this trajectory contain?

**Why the answer is a time–band limiting problem.** The overlap of two points of the trajectory is the survival amplitude, ⟨ψ(t)|ψ(s)⟩ = ∫ e^{−iE(t−s)/ħ} dμ(E) = A(t − s), the Fourier transform of the energy distribution. A is band-limited to |ω| ≤ Ω\_E/2ħ (after centring the energy), and it is observed on a window of length Δτ. The Gram operator of the trajectory on L²\[0, Δτ\], (Gf)(t) = ∫₀^{Δτ} A(t − s) f(s) ds, is therefore a time–band limiting operator in the sense of Slepian, Landau and Pollak, with Slepian parameter

```latex
c=\frac{\Omega_E}{2\hbar}\cdot\frac{\Delta\tau}{2}=\frac{\Omega_E\,\Delta\tau}{4\hbar}.
```

**The law (flat energy distribution).** Let the energy distribution be uniform on its band. Normalise G by W/π, W = Ω\_E/2ħ; its eigenvalues λ₀ > λ₁ > … are then Slepian's concentration eigenvalues, in (0, 1), and d\_ε(Δτ) = #{n : λ\_n > ε} is the number of modes of the trajectory distinguishable at level ε. Then

```latex
d_\varepsilon(\Delta\tau)=\underbrace{\frac{\Omega_E\,\Delta\tau}{h}}_{\text{Landau count}}\;+\;\underbrace{\frac{1}{\pi^{2}}\ln\!\Big(\frac{\Omega_E\Delta\tau}{4\hbar}\Big)\ln\frac{1-\varepsilon}{\varepsilon}}_{\text{soft horizon}}\;+\;o(\ln c),
```

and beyond the plunge the eigenvalues fall super-exponentially: for every K ≥ 2c,

```latex
\sum_{n\ge K}\lambda_n\;\le\;\sum_{k\ge K}\frac{2}{\pi}\,\frac{c^{2k+1}}{\big((2k+1)!!\big)^{2}},
```

(Section 5, where the bound is proved for every K and every bounded energy density, with the density's maximum as a prefactor). For a non-flat distribution with density bounded between m and M on the band, the count is sandwiched between the flat counts at levels ε/m and ε/M (Theorem (b)); the leading term is in general (Δτ/2πħ) × the measure of the energy set where the density exceeds ε (Widom's Szegő-type theorem), which reduces to Ω\_EΔτ/h when ε is below the density minimum, and the logarithmic term is set by the jumps of the density at the band edges and at interior discontinuities.

**Reading.** The leading term is the number of degrees of freedom of a signal of bandwidth Ω\_E/h observed for Δτ — Landau's 2WT theorem \[LP62\] — and it equals the action budget of the process in units of Planck's constant h. This is Axiom 3 made quantitative, and it is of the same form as the quantum speed limits, which bound the number of successive orthogonal states on the trajectory: for the flat band starting at the ground state the Margolus–Levitin bound is 2Ω\_EΔτ/h and the Mandelstam–Tamm bound is 4σ\_EΔτ/h ≈ 1.15 Ω\_EΔτ/h (σ\_E = Ω\_E/√12). The two kinds of count differ in meaning: the speed limits count orthogonalisations along the curve (a two-level system saturates Margolus–Levitin while its trajectory never leaves a plane), the Landau count is the ε-rank of the span of the trajectory, that is, of the time-averaged state (1/Δτ)∫|ψ(t)⟩⟨ψ(t)|dt. The second term is what CST's axioms did not contain: distinguishability does not stop at the count. There is a transition region of width (2/π²) ln c · ln((1 − ε)/ε) modes in which states are partially distinguishable, and the region grows only logarithmically with the action budget. This is the soft horizon. It has no adjustable constant, and it is the same plunge that governs, in the companion paper, which zeros of ζ the primes up to λ² determine.

**The cascade number.** With ΔE identified with Ω\_E, κ = ΔEΔτ/(ħ ln 2) = (2π/ln 2) × (Landau count): κ counts the same quantity in a different unit, 9.06 κ-units per distinguishable state. For the flat band the first orthogonal state on the trajectory (the first zero of the survival amplitude) appears at a Landau count of 1, i.e. κ ≈ 9.06, not at κ = 1; and the second mode of the trajectory becomes distinguishable through a plunge, not a threshold. Prediction P2 of the roadmap (a sharp transition at κ = 1) is therefore withdrawn and replaced by the law. The archive's cascade time Δτ = ħ ln 2/ΔE (22 ps for a 5 GHz transmon) is the time at which the Landau count reaches ln 2/2π ≈ 0.11: the trajectory is still inside its first mode, and nothing discrete happens there.

## 5. Theorem: the ε-dimension of a quantum trajectory

**Definitions.** Let H be self-adjoint on a Hilbert space, ψ a unit vector with spectral measure μ supported in \[E₀, E₀ + Ω\_E\], and ψ(t) = e^{−iHt/ħ}ψ. For Δτ > 0 let W = Ω\_E/2ħ, T = Δτ/2, c = WT, and let μ̃ be the centred, rescaled measure on \[−W, W\]. The Gram operator G\_μ on L²\[−T, T\] has kernel A(t − s) = ∫ e^{−iω(t−s)} dμ̃(ω); it is the Toeplitz operator with symbol μ̃ compressed to the window. When μ̃ is the uniform distribution on \[−W, W\], G\_μ = (π/W)·P\_T B\_W P\_T, where P\_T B\_W P\_T is Slepian's time–band limiting operator (time-limit to \[−T, T\], band-limit to \[−W, W\]), whose eigenvalues λ₀ > λ₁ > … lie in (0, 1) and sum to 2c/π. We normalise G\_μ by the same factor W/π and write λ\_n(μ) for the eigenvalues of (W/π)G\_μ; for a density ρ = dμ̃/dω this is the compression to the window of the multiplier 2Wρ(−ω) (the reflection is immaterial below), so the eigenvalues lie in (0, sup 2Wρ) and sum to 2c/π for every probability measure μ̃. The ε-dimension of the trajectory is d\_ε(μ, Δτ) = #{n : λ\_n(μ) > ε}; it is the ε-rank of the time-averaged state (1/Δτ)∫₀^{Δτ}|ψ(t)⟩⟨ψ(t)|dt, whose nonzero eigenvalues are those of G\_μ/Δτ. For a discrete spectrum with N\_lev populated levels the rank is at most N\_lev, so d\_ε ≤ N\_lev for every Δτ and the statements below describe the regime in which the window does not resolve the levels, Δτ × (level spacing)/h ≪ 1, equivalently N\_lev ≫ Ω\_EΔτ/h.

**Theorem.** (a) *Flat energy distribution.* If μ̃ is uniform on \[−W, W\], then for every fixed ε ∈ (0, 1), as c → ∞,

```latex
d_\varepsilon=\frac{2c}{\pi}+\frac{1}{\pi^{2}}\ln c\,\ln\frac{1-\varepsilon}{\varepsilon}+o(\ln c).
```

(b) *Bounded density.* If m ≤ 2Wρ(ω) ≤ M on \[−W, W\], then for every ε and every Δτ,

```latex
d_{\varepsilon/m}^{\mathrm{flat}}\;\le\;d_\varepsilon(\mu,\Delta\tau)\;\le\;d_{\varepsilon/M}^{\mathrm{flat}},
```

so for ε < m the Landau count and the logarithmic plunge are bounded on both sides by flat counts; to leading order the sandwich has width (1/π²) ln c · \[ln(M/m) + ln((1 − ε/M)/(1 − ε/m))\] modes. The bound is not the actual behaviour: by Widom's Szegő-type theorem for compressed multipliers the leading term of d\_ε(μ) is (T/π)·|{ω : 2Wρ(ω) > ε}|, which equals 2c/π exactly when ε < m, and the logarithmic term is (ln c/2π²)\[ln((b₊ − ε)/ε) + ln((b₋ − ε)/ε)\] with b± the values of 2Wρ at the band edges, plus one such term for each interior jump of ρ; a density vanishing continuously at the edges has no logarithmic term. For m = 0 the lower bound is vacuous: a distribution supported on two half-bands has d\_ε ≈ c/π for all ε < 2, half the Landau count, although its trace is still 2c/π.

(c) *Hard horizon.* For every μ as in (b) and every K ≥ 1,

```latex
\sum_{n\ge K}\lambda_n(\mu)\;\le\;M\sum_{k\ge K}\frac{2}{\pi}\frac{c^{2k+1}}{\big((2k+1)!!\big)^{2}}\;\le\;M\,e^{-1.5\,c}\quad\text{for }K\ge 2c,
```

so no mode beyond index 2c is distinguishable at any level above e^{−1.5c}.

**Proof.** (a) is the theorem of Landau and Widom \[LW80\] on the eigenvalue distribution of time and frequency limiting, since (W/π)G\_μ = P\_T B\_W P\_T exactly. (b) The operator (W/π)G\_μ is P\_T F⁻¹ M\_{2Wρ} F P\_T, and m·1\_{\[−W,W\]} ≤ 2Wρ ≤ M·1\_{\[−W,W\]} as multipliers gives m·P\_T B\_W P\_T ≤ (W/π)G\_μ ≤ M·P\_T B\_W P\_T in the operator order; the eigenvalue counting function is monotone in the operator order (Weyl's inequality), which is the two-sided bound. (c) By the Ky Fan minimum principle, the sum of all eigenvalues beyond the K-th of a positive compact operator is the minimum, over K-dimensional subspaces V, of the trace of the operator compressed to V^⊥. Take V = polynomials of degree < K on \[−T, T\]; the orthonormal Legendre polynomials p\_k, k ≥ K, span V^⊥ and the trace is Σ\_{k≥K} ⟨p\_k, (W/π)G\_μ p\_k⟩ ≤ M Σ\_{k≥K} ℓ\_k, where ℓ\_k = ⟨p\_k, P\_T B\_W P\_T p\_k⟩ is the in-band energy of p\_k. Since the Fourier transform of p\_k(t) = √((2k+1)/2T) P\_k(t/T) is √(2T(2k+1)) (−i)^k j\_k(ωT) with j\_k the spherical Bessel function, ℓ\_k = (2k+1)(2/π)∫₀^c j\_k(x)² dx (so that Σ\_k ℓ\_k = 2c/π), and the bound |J\_ν(x)| ≤ (x/2)^ν/Γ(ν + 1) (ν ≥ −½, x > 0) gives j\_k(x) ≤ x^k/(2k + 1)!!, hence ℓ\_k ≤ (2/π)c^{2k+1}/((2k+1)!!)². With (2k+1)!! ≥ √2·((2k+1)/e)^{k+½} this is at most (1/π)(ec/(2k+1))^{2k+1}, which for 2k+1 ≥ 4c is at most (1/π)(e/4)^{4c} < e^{−1.54c}; the terms decrease geometrically with ratio below (e/4)² = 0.46, so the sum is below 0.6 e^{−1.54c} < e^{−1.5c}, and K ≥ 2c gives 2K + 1 > 4c. ∎

**Remarks.** (1) The theorem is an assembly of known results — Landau–Widom, Widom's Szegő-type theorem, Weyl monotonicity, Ky Fan, a Bessel inequality — around one identification: the Gram operator of a quantum trajectory is a time–band limiting operator with the energy distribution as symbol. The quantity itself is not new: the time-averaged state and its effective dimension 1/Tr ω² are standard in equilibration theory \[Re08, LPSW09\], and the span of the trajectory is the Krylov space of the survival amplitude, whose occupation is the spread complexity of \[BCMW22\]. What we have not found written is the identification with Slepian's operator, which is what brings the Landau count, the Landau–Widom plunge and the explicit cutoff of part (c) to that dimension; the claim of novelty is restricted to that step. The quantum speed limits \[MT45, ML98, DC17\] bound orthogonalisations along the curve and do not describe the plunge. (2) Part (c) is the hard horizon of the law and is fully rigorous for every c and every bounded density; part (a) is asymptotic, with an o(ln c) error that is not uniform in ε \[MRS23\], and Section 6 measures how far the finite-c counts sit from it. (3) Nothing in the proof uses CST; the theorem is a statement of quantum mechanics. CST's contribution is to have singled out the cascade number as the quantity to look at, and the theorem tells what that quantity actually measures.

## 6. Numerical verification

The prolate eigenvalues λ\_n(c) were computed by Bouwkamp's tridiagonal Legendre method (prolate.py, mpmath, 30–50 digits; the sum Σλ\_n = 2c/π is reproduced to 15 digits at every c). The table gives the measured count d\_ε = #{n : λ\_n > ε} against the Landau–Widom value 2c/π + (1/π²) ln c · ln((1 − ε)/ε), and the hard-horizon bound of Theorem (c) at K = 2c against the true tail sum.

| c | Landau count 2c/π | d₀.₅ (LW) | d₁₀⁻² (LW) | d₁₀⁻⁴ (LW) | d₁₀⁻⁶ (LW) | d₁₀⁻¹⁰ (LW) | tail bound at K = 2c | true tail Σ\_{n≥2c}λ\_n |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 6.37 | 6 (6.4) | 9 (7.4) | 10 (8.5) | 12 (9.6) | 15 (11.7) | 3.1 × 10⁻⁷ | 3.9 × 10⁻¹⁰ |
| 30 | 19.10 | 19 (19.1) | 22 (20.7) | 24 (22.3) | 26 (23.9) | 30 (27.0) | 2.9 × 10⁻²⁰ | 5.1 × 10⁻²⁴ |
| 100 | 63.66 | 64 (63.7) | 67 (65.8) | 70 (68.0) | 72 (70.1) | 77 (74.4) | 7.2 × 10⁻⁶⁶ | 1.7 × 10⁻⁷¹ |

Three things are visible. The Landau count is exact at ε = ½ to within one mode at every c. The plunge is real and logarithmic: between ε = ½ and ε = 10⁻¹⁰ the count grows by 9 modes at c = 10 and 13 at c = 100, while the count itself grows tenfold. The Landau–Widom asymptotic written with ln c underestimates the plunge at finite c by 2–3 modes already at ε = 10⁻², and the discrepancy does not grow with c; written with ln 4c (the logarithm of the product of the two interval lengths, which is how Landau and Widom parametrise it, the difference being o(ln c)) it reproduces every entry of the table to within one mode down to ε = 10⁻⁶ (c = 10: 8.1, 9.8, 11.5, 15.0; c = 100: 66.5, 69.3, 72.1, 77.6). The asymptotic has an unbounded o(ln c) error and no uniformity in ε \[MRS23\]; at c of order 10 the computed counts, not the formula, are the prediction. The hard horizon holds with room: the bound e^{−1.5c} at K = 2c exceeds the true tail by three to five orders of magnitude, and the true tail is already 10⁻⁷¹ at c = 100, so beyond twice the Slepian parameter nothing is distinguishable at any physical level.

For the arithmetic instance of Section 7 (c = ΩL = 55.5 at λ = 2, Ω = 80): even modes fall below 10⁻²⁰ at index 56 and below 10⁻³⁰ at index 64 (Landau count 35.3), and the Legendre tail bound below degree 110 is 1.08 × 10⁻⁴²; that bound is what makes the certified theorem of the companion paper a finite computation.

## 7. The arithmetic instance

The companion paper *Scale-Resolved Weil Positivity* (Zenodo 10.5281/zenodo.23229271; code at github.com/davidangularme/scale-resolved-weil-positivity) studies the Weil quadratic form W\_λ restricted to test functions supported in \[−log λ, log λ\], which depends only on the prime powers up to λ². Under RH it equals Σ\_γ |f̂(γ)|² over the zeros of ζ: the "trajectory" is the set of zeros read as frequencies, the "window" is the support, and the same time–band limiting decides what is resolved. The dictionary:

| Section 4–5 (a quantum trajectory) | Companion paper (the Weil form) |
| --- | --- |
| window length Δτ | 2L = 2 log λ |
| band Ω\_E/ħ, flat or bounded density | zero density (1/2π) log(γ/2π), unbounded, increasing |
| Landau count Ω\_EΔτ/h | resolution height H(λ) = 2πλ², where the density reaches L/π |
| λ\_n(μ) > ε: mode distinguishable | zero below H reconstructed from the primes to exp(−cH) |
| soft horizon, width (2/π²) ln c · ln(1/ε) | range law: a quarter digit lost per unit height above H, horizon near 2H |
| hard horizon e^{−1.5c} beyond K = 2c | the Legendre tail bound 10⁻⁴² that closes Theorem 2 |
| Gram operator = Toeplitz with symbol μ | W\_λ = Toeplitz with the Weil symbol σ\_λ(r) = Re ψ(¼ + ir/2) − log π − 2ΣΛ(n)n^{−1/2}cos(r log n) |

**What transfers.** The structure: a count fixed by the window and the band, a logarithmic plunge, a super-exponential tail, and the fact that positivity or distinguishability is decided by a finite block plus a tail that time–band limiting controls. The certified theorem of the companion paper (W₂ ≥ 1.38 × 10⁻¹³ on the whole window) is the hard-horizon mechanism of Theorem (c) applied to an arithmetic symbol.

**What does not transfer.** The numbers. The zero density grows like log γ, so the arithmetic plunge is exponential in height (a quarter digit per unit) rather than logarithmic in mode index, and the horizon sits at about 2H rather than at a logarithmic distance from the count; the resolution height 2πλ² is a property of ζ, not of any physical system. In the other direction, nothing in the law says anything about where the zeros are: the Weil form's positivity for all λ is RH, and the law only describes what a finite window can see of it. The archive's earlier claim that CST's spectrum *is* the zeros is withdrawn (Section 2) and the law does not revive it.

**What the arithmetic instance taught the physics.** Two things that we would not have looked for otherwise. First, the minimiser of the Weil form is Riemann's Φ times a smooth window — the finite-scale object is the global object tapered, which is what the trajectory's prolate modes are to the full orbit. Second, positivity on a window pins the data generating the symbol to many digits (eleven digits of Λ(2) at λ = 2, certified): the physical counterpart is that the Gram spectrum of a trajectory over a window determines the energy distribution to a precision set by the plunge, a rigidity statement that Section 8 turns into a measurement protocol.

## 8. What is new, and how to falsify it

**In physics.** (1) The cascade number is identified: κ = (2π/ln 2) × (action budget in units of h) = (2π/ln 2) × (Landau count of the trajectory), a quantity of the same kind as the quantum speed limits (which, for the flat band, are 1.15 to 2 times the Landau count) but measuring the ε-rank of the span of the trajectory rather than orthogonalisations along it; the archive's reading of κ = 1 as a reversibility threshold (P2) is withdrawn. (2) The speed limits have a penumbra: beyond the Landau count a trajectory still contains, for a flat band, about (1/π²) ln(4c) · ln((1 − ε)/ε) modes distinguishable at level ε, and none beyond twice the Slepian parameter at any level above e^{−1.5c}. For an action budget of 10⁶ h and ε = 10⁻³ this is about 10 extra modes on a count of 10⁶; for a budget of 10 h (c = 15.7) the computed count at ε = 10⁻³ is 14 against a Landau count of 10, four extra modes, a 40 % effect. The penumbra is where the law is visible. (3) Time is not granular at the cascade time ħ ln 2/ΔE (roadmap prediction on discrete time): at that time the trajectory has a Landau count of 0.11 and sits inside its first mode.

**In mathematics.** The identification of the Gram operator of a quantum trajectory with a Toeplitz operator whose symbol is the energy distribution, and the counting law that follows (Theorem, Section 5), with an explicit non-asymptotic hard-horizon bound. The ingredients are classical and the quantity (the effective dimension of the time-averaged state, the Krylov span of the trajectory) is studied elsewhere; what is new, to our knowledge, is the identification with Slepian's operator and what it brings (Section 5, Remark 1). The open mathematical question it leaves is the same one the companion paper leaves on the arithmetic side: a quantitative inequality for the plunge at finite c (a Remez or Logvinenko–Sereda type bound on Paley–Wiener spaces), which would turn the o(ln c) into an explicit error.

**The measurement.** The law is a statement about a trajectory whose energy distribution has many levels inside the band: the ε-rank is at most the number of populated levels, so a qubit (rank 2) or a few-level transmon cannot show a Landau count above its level number, and the regime of the law is N\_lev ≫ Ω\_EΔτ/h, the window shorter than the inverse level spacing. The natural systems are a harmonic mode in a broad coherent or thermal state (a superconducting cavity, a trapped-ion motional mode), a spin ensemble, or a many-body quantum simulator, where tens to hundreds of levels are populated and the survival amplitude A(t) = ⟨ψ(0)|ψ(t)⟩ is accessible by Ramsey or Loschmidt-echo interferometry. Protocol: (i) prepare a state with known populations, so that μ is known, Ω\_E is the spread of the populated levels and N\_lev is their number (for a cavity with populated levels spanning 5 GHz, the Landau count reaches 1 at Δτ = 0.2 ns, 10 at 2 ns, 100 at 20 ns, so N\_lev of a few hundred is needed to see the plunge at c of order 100); (ii) measure A(t) on a grid t₁ < … < t\_m in \[0, Δτ\]; (iii) form the Hermitian Toeplitz matrix A(t\_i − t\_j), normalise by W/π and the grid spacing, and compute its eigenvalues; (iv) count the eigenvalues above ε for ε = ½, 10⁻¹, 10⁻², and repeat for several Δτ so that c varies by a factor of ten or more. The law predicts the count at ε = ½ equal to Ω\_EΔτ/h within one mode, and an excess at small ε that grows like the logarithm of c, with the computed values of Section 6 as the prediction for a flat band: at ε = 10⁻² the excess is 3 modes at c = 10, 30 and 100 (its logarithmic growth is below one mode over that range), at ε = 10⁻⁶ it is 6, 7 and 8 modes, so the scaling test needs small ε and c spanning at least a decade. Decoherence does not lower these counts, it raises them: a decay factor on A(t) leaves A(0) = 1 and the trace unchanged, and its spectral tails add modes at small ε — an excess growing linearly in c for Gaussian dephasing, and far larger for exponential dephasing (at c = 100, ε = 10⁻⁴, a Gaussian factor exp(−t²/T²) gives 72 modes against 70, an exponential factor exp(−t/2T) gives hundreds). The discriminator is therefore the c-scaling of the excess: logarithmic is the law, linear or faster is dephasing, and a sharp cutoff is neither. The law is falsified if, after dephasing is accounted for by its measured scaling, the ε = ½ count departs from Ω\_EΔτ/h by more than one mode or the small-ε excess does not grow logarithmically. The same data determine μ from the Gram spectrum (the Toeplitz inverse problem), which is the physical rigidity statement of Section 7.

**The existing qubit data.** The 768-pair diagnostic of v18 measures a correlation between a Witten-Laplacian gap and gate fidelity across calibration snapshots; it contains no time series and cannot test the law. It is kept as what it is, a diagnostic, and the law is a new prediction, not a reinterpretation of that result.

**What is not claimed.** No statement about gravity, dark matter, gravitational waves, or the Standard Model (roadmap items P6, P9 and the derivations of Part II are untouched and untested). No statement about the Riemann Hypothesis. The law is a theorem of quantum mechanics about what a finite window of a trajectory can distinguish, which CST's third axiom had stated in its sharp half only.

## 9. Reproducibility, record text, references

**Code and data (v33\_cst\_code.zip).** prolate.py computes the prolate eigenvalues λ\_n(c) by Bouwkamp's tridiagonal Legendre method; lw\_check.py produces the table of Section 6 (run: python3 lw\_check.py, about two minutes; c = 300 needs a truncation of 360 modes per parity); janus.py runs the three-body complexity test of Section 3′ (run: python3 janus.py SEED 3). The audit computations of Section 2 and the companion paper's scripts are in the repository github.com/davidangularme/scale-resolved-weil-positivity (code/, with the ball-arithmetic certificates). Requirements: Python 3, mpmath, python-flint, numpy, scipy.

**Zenodo record text for this version.** *Title:* Configuration Space Temporality — Complete Research Archive, v33: Corrections and a Resolution Law. *Description:* Version 33 withdraws every claim of this archive that the Riemann Hypothesis has been proved, reduced to one assumption, or approached through the Bost–Connes Hamiltonian, the Fisher–Rao action or a rigidity matrix; a re-audit with every number recomputed found four errors, documented in this version with the status of all 17 RH-related claims. The physical framework (time as the cost of transitions in configuration space, bounded transition rate, cascade number) is kept and restated. The new content is one law: the number of states that a quantum trajectory of duration Δτ and energy spread Ω\_E makes distinguishable at level ε equals, for a flat energy band, the action budget Ω\_EΔτ/h (the Landau count, a quantity of the same kind as the quantum speed limits, which identifies the cascade number) plus a logarithmic plunge (1/π²) ln(Ω\_EΔτ/4ħ) ln((1 − ε)/ε) set by the band edges, with a rigorous super-exponential cutoff beyond twice the Slepian parameter for every bounded energy density. The law is proved as a theorem of time–band limiting, verified numerically, and comes with a falsification protocol on transmon qubits. It replaces the earlier prediction of a sharp reversibility threshold at κ = 1. The law was found through the companion preprint *Scale-Resolved Weil Positivity* (10.5281/zenodo.23229271), where the same harmonic analysis fixes which zeros of ζ the primes up to λ² determine; nothing here bears on the Riemann Hypothesis. Earlier files remain for the record; the superseded ones are listed in Section 2. *Keywords:* Configuration Space Temporality; emergent time; quantum speed limit; Margolus–Levitin; time–band limiting; prolate spheroidal functions; Landau–Widom; Bures metric; transmon qubits; Weil positivity. *Related identifiers:* 10.5281/zenodo.23229271 (is supplemented by); github.com/davidangularme/scale-resolved-weil-positivity (is supplemented by).

**References.**

- &#91;LW80\] H. J. Landau, H. Widom, Eigenvalue distribution of time and frequency limiting, J. Math. Anal. Appl. 77 (1980), 469–481.
- &#91;SP61\] D. Slepian, H. O. Pollak, Prolate spheroidal wave functions, Fourier analysis and uncertainty I, Bell Syst. Tech. J. 40 (1961), 43–63.
- &#91;La67\] H. J. Landau, Necessary density conditions for sampling and interpolation of certain entire functions, Acta Math. 117 (1967), 37–52.
- &#91;ML98\] N. Margolus, L. B. Levitin, The maximum speed of dynamical evolution, Physica D 120 (1998), 188–195.
- &#91;MT45\] L. Mandelstam, I. Tamm, The uncertainty relation between energy and time in non-relativistic quantum mechanics, J. Phys. (USSR) 9 (1945), 249–254.
- &#91;DC17\] S. Deffner, S. Campbell, Quantum speed limits: from Heisenberg's uncertainty principle to optimal quantum control, J. Phys. A 50 (2017), 453001.
- &#91;CM14\] E. A. Carlen, J. Maas, An analog of the 2-Wasserstein metric in non-commutative probability under which the fermionic Fokker–Planck equation is gradient flow for the entropy, Comm. Math. Phys. 331 (2014), 887–926.
- &#91;Bl26a\] F. D. Blum, Scale-Resolved Weil Positivity: Reconstruction of Zeta Zeros from Finitely Many Primes, a Detection Law, and a Certified Block, Zenodo, 10.5281/zenodo.23229271 (2026).
- &#91;Bl26b\] F. D. Blum, CST/FDBC Archive v33 — Corrections and a Non-Circular Weil-Form Test (the full audit; in this record).
- &#91;Zh26\] X. Zhu, Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau–Widom decay law, arXiv:2608.24827 (2026).
- &#91;CST18\] F. D. Blum, Configuration Space Temporality v18 (Physical Review D draft) and What Does CST Predict? (February 2026 roadmap), both in this record.
- &#91;PT21\] D. J. Platt, T. S. Trudgian, The Riemann hypothesis is true up to 3·10¹², Bull. Lond. Math. Soc. 53 (2021), 792–797.
- &#91;LP62\] H. J. Landau, H. O. Pollak, Prolate spheroidal wave functions, Fourier analysis and uncertainty III: the dimension of the space of essentially time- and band-limited signals, Bell Syst. Tech. J. 41 (1962), 1295–1336.
- &#91;Wi60\] H. Widom, A theorem on translation kernels in n dimensions, Trans. Amer. Math. Soc. 94 (1960), 170–180 (Szegő-type eigenvalue distribution of compressed multipliers).
- &#91;MRS23\] F. Marceca, J. L. Romero, M. Speckbacher, Eigenvalue estimates for Fourier concentration operators on two domains, [arXiv:2301.09616](https://arxiv.org/abs/2301.09616) (2023).
- &#91;Re08\] P. Reimann, Foundation of statistical mechanics under experimentally realistic conditions, Phys. Rev. Lett. 101 (2008), 190403.
- &#91;LPSW09\] N. Linden, S. Popescu, A. J. Short, A. Winter, Quantum mechanical evolution towards thermal equilibrium, Phys. Rev. E 79 (2009), 061103.
- &#91;BCMW22\] V. Balasubramanian, P. Caputa, J. M. Magan, Q. Wu, Quantum chaos and the complexity of spread of states, Phys. Rev. D 106 (2022), 046007, [arXiv:2202.06957](https://arxiv.org/abs/2202.06957).
- &#91;BB82\] J. B. Barbour, B. Bertotti, Mach's principle and the structure of dynamical theories, Proc. R. Soc. Lond. A 382 (1982), 295–306.
- &#91;BFO02\] J. Barbour, B. Z. Foster, N. Ó Murchadha, Relativity without relativity, Class. Quantum Grav. 19 (2002), 3217–3248.
- &#91;BKM14\] J. Barbour, T. Koslowski, F. Mercati, Identification of a gravitational arrow of time, Phys. Rev. Lett. 113 (2014), 181101.
- &#91;PW83\] D. N. Page, W. K. Wootters, Evolution without evolution: dynamics described by stationary observables, Phys. Rev. D 27 (1983), 2885–2892.
- &#91;MGV14\] E. Moreva, G. Brida, M. Gramegna, V. Giovannetti, L. Maccone, M. Genovese, Time from quantum entanglement: an experimental illustration, Phys. Rev. A 89 (2014), 052122.
- &#91;CR94\] A. Connes, C. Rovelli, Von Neumann algebra automorphisms and time–thermodynamics relation in generally covariant quantum theories, Class. Quantum Grav. 11 (1994), 2899–2917.
