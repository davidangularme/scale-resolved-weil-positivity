# CST/FDBC Archive v33 — Corrections and a Non-Circular Weil-Form Test

Oct 7, 2026 · @Frédéric David Blum

## Summary

Version 33 withdraws every claim in this archive that the Riemann Hypothesis has been proved or reduced to a single open assumption. A re-audit of the core documents, with every number recomputed, found four errors that break all three proposed routes to RH. The one positive result is new and non-circular: the Weil quadratic form built from the primes 2, 3, 5 and 7 alone, in the archive's own cosine basis, has a minimiser whose Fourier transform vanishes at the first 17 zeta zeros to between 8 and 36 digits. A new section, "Scale-resolved Weil positivity", states exactly what positivity at a finite scale λ can say about zeros (everything below height 2πλ², almost nothing above), verifies it numerically, and sets the first provable targets beyond the known λ = √2. A further section measures the block itself: the positivity margin at scale λ is below 10⁻³⁶ at λ = 3 and below 10⁻⁷⁴ at λ = 5, so positivity pins the prime weights to tens of digits and no perturbative route exists. Three new objects follow: pseudo-zero measures (Krein extension), the prime-driven Krein system (induction on primes done exactly), and the compressed prime mass (which, it turned out, cannot move the certification frontier: see "The way"). The last section measures the range law: the primes up to λ² fix the zeros up to 2πλ² exactly (the first zero to 105 digits from the primes up to 25), and their information fades to nothing near twice that height. The minimiser itself turns out to be Riemann's Φ cut to the window \[−log λ, log λ\] by a smooth taper, which gives the range law the form of an equation. The work was then written up as a separate paper, "Scale-Resolved Weil Positivity", with four theorems: ball-arithmetic positivity of W\_2 on the whole window \[1/2, 2\] with the margin enclosed between 1.38 × 10⁻¹³ and 7.6 × 10⁻¹³; a certified rigidity statement (the window alone pins Λ(2) to eleven digits and log 2 to twelve); and, under RH, coercivity at every scale and localisation of the minimiser's zeros. A literature search made at the end found that the positivity theorem had been obtained independently and more generally (support up to L = 0.8) by Xuefeng Zhu in August 2026, with the same reduction and the same resolution height; the paper cites him and states what is ours.

The four corrections:

1. The Fisher pole identity |I(s)| · |s − ρ|² = 1 holds at a simple zero of any analytic function, on or off the critical line. It carries no information about RH.
2. The rigidity matrix R has an exact zero eigenvalue, so the confinement constant α is undefined. The transverse Hessian of the log-gas energy is −R, so the critical line is a transverse maximum, not a minimum.
3. The "exponential dominance" of the Weil minimiser is false. Extended from N = 45 to N = 200, α(N) falls as roughly 1/N, from 0.143 to 0.039. GUE, Poisson and evenly spaced points give the same values to within 1 %.
4. Assumption A asks an operator with spectrum {log n} to have eigenvalues 1/4 + γₙ². It cannot hold as stated. The CMP v25–v27 "complete proof chain" is withdrawn.

The physical CST programme (time as optimal transport on the Bures manifold, the IBM transmon analysis) is outside the scope of this audit. Its own February 2026 status report already records what held and what fell there.

## Status of every RH-related claim

Of the 17 RH-related results audited, 3 hold as stated, 3 hold but are trivial or already known, 1 is falsified numerically, 1 stays open and 9 are withdrawn. The "v32 status" column is the archive's own label.

| Result | v32 status | v33 status | Reason |
| --- | --- | --- | --- |
| Fisher pole quantization, abs(I) · abs(s − ρ)² = 1 | Proved | Holds, trivial | True at any simple zero of any analytic function (Correction 1) |
| Coulomb reduction of the Fisher–Rao action | Exact | Withdrawn | 1/r is the 3D Coulomb law, not 2D; the action diverges on any path through a zero |
| Rigidity matrix R and constant α | Open, α ≈ 0.14 | Withdrawn | R has an exact zero eigenvalue, so μ₁/μ₀ = ∞ (Correction 2) |
| α∞ > 0 ⇔ RH (Thm 4.3, Conj 7.2) | Conjecture | Withdrawn | R is built from ordinates already on the line; it cannot see an off-line zero |
| Generalized Stieltjes principle | Proved | Holds, vacuous | Shows the symmetry axis is a critical path, not that zeros lie on it |
| Polyak–Łojasiewicz inequality | Proved | Holds | Standard Rayleigh-quotient bound for any symmetric matrix |
| Exponential isolation table, N ≤ 45 | Verified | Holds | Reproduced to two decimals in v33 |
| Spectral gap Δ(N) > 0 (W1) | Verified, N ≤ 45 | Holds | Checked to N = 200; same for random points (Correction 3) |
| Exponential dominance μ₀/Δ → 0 (Conj ED, W2) | Open, numerical | Falsified | α(N) ≈ c/N, from 0.143 to 0.039 (Correction 3) |
| Weil matrix M as "the Weil form" | Assumed | Withdrawn | M is the zero-side sum over 500 zeros, positive by construction; replaced by the prime-side W (New result) |
| Assumption A: det\_reg(H\_BC − s(1 − s)) = C · ξ(s) | Open | Withdrawn | Spectrum of H\_BC is {log n}, not {1/4 + γₙ²} (Correction 4) |
| Assumption B: spectral rigidity via Galois and KMS | Proved | Withdrawn | H\_mod is self-adjoint with real spectrum; KMS does not give s ↔ 1 − s (Correction 4) |
| Zeros of ζ as poles of the resolvent of H\_BC (CMP v27, §8) | Proved | Withdrawn | The resolvent's poles are at log n (Correction 4) |
| CMP v25–v27 complete proof chain | Partly retracted | Withdrawn | Krein–Rutman step retracted in v29; the files still claimed a proof |
| Fisher–Bakry–Émery bridge (Path B) | Open | Withdrawn | Petz (1996) concerns finite matrix spaces; no curvature bound for type III₁ states is cited |
| Li criterion bridge (Thm 5.1) | Known | Holds, known | Li (1997); equivalent to RH, not progress toward it |
| BDC: Li positivity for n > 34 (O3) | Open | Open | Equivalent to RH; unchanged |

## Correction 1: the Fisher pole identity is generic

The identity that v32 calls a quantization condition counts the multiplicity of a zero, for any analytic function. The archive's Fisher information is the second logarithmic derivative of ζ:

```latex
I(s) = \frac{\zeta''(s)}{\zeta(s)} - \left(\frac{\zeta'(s)}{\zeta(s)}\right)^{2} = \frac{d^{2}}{ds^{2}} \log \zeta(s)
```

If f is analytic with a zero of order m at ρ, then f(s) = (s − ρ)ᵐ g(s) with g(ρ) ≠ 0, and

```latex
\frac{d^{2}}{ds^{2}} \log f(s) = -\frac{m}{(s-\rho)^{2}} + \frac{d^{2}}{ds^{2}} \log g(s)
```

So |I(s)| · |s − ρ|² → m at every zero, on or off the critical line. v33 checked this on f(s) = (s − a)(s − b)eˢ with a deliberately off-line zero at a = 0.8 + 20i: the product equals 1.0000000 at |s − ρ| = 10⁻⁸. At the first zeta zero it also equals 1.0, as v31 reported. The "unit information charge" is the statement that zeros are simple; it would hold just as well if off-line zeros existed.

Two further corrections to §2 of the v32 paper:

- For complex s, the weights pₙ(s) = n⁻ˢ/ζ(s) are complex numbers, not a probability distribution. I(s) is a Fisher information only on the real axis s > 1.
- The Bost–Connes system has a unique KMS state for 0 < β ≤ 1. For β > 1 the extremal KMS states form a family (spontaneous symmetry breaking), the reverse of what §2.1 states.

The remark "off the critical line I(s) remains bounded" (v31, §4.5) is true of every point that is not a zero, wherever the zeros are. It is not evidence for RH.

## Correction 2: the rigidity matrix cannot confine

The confinement constant α of v32 is undefined, and the matrix it is built on points the wrong way. Every row of the rigidity matrix sums to zero by construction:

```latex
\sum_{m} R_{nm} = \sum_{k \neq n} \frac{1}{(\gamma_n-\gamma_k)^{2}} - \sum_{m \neq n} \frac{1}{(\gamma_n-\gamma_m)^{2}} = 0
```

So the constant vector is an eigenvector with eigenvalue exactly 0, and μ₁/μ₀ is infinite for every N. With the first 45 zeros at 60 digits, v33 finds μ₀ = −1.5 × 10⁻⁶¹ (zero to working precision) and μ₁ = 0.0186. The value α ≈ 0.14 quoted in v32 is not a property of R at all: it was taken from the Weil-matrix table of v31, a different matrix (Correction 3).

The sign is also reversed. For the two-dimensional log-gas energy E = −Σ log|zₙ − zₘ| with zₙ = 1/2 + xₙ + i(γₙ + yₙ), the Hessian along the line (in y) is +R, and the Hessian across it (in x) is −R. v33 confirmed this by finite differences: the transverse entry at (0, 0) is −0.0493864, the longitudinal one +0.0493864, and R₀₀ = 0.0493864. This is forced: log|z| is harmonic, so the two Hessians sum to zero (Earnshaw's theorem). The collinear configuration is a transverse maximum of the energy; it does not confine zeros to the line.

Two more points close this route:

- R is assembled from the ordinates γₙ of zeros already placed on the line. An off-line zero has no entry in it, so no property of R can be equivalent to RH. Theorem 4.3 of v32 is withdrawn.
- Near a zero, √|I(s)| behaves like 1/|s − ρ|. That is the 3D Coulomb law, not the 2D one (log|s − ρ|), and its integral along any path through ρ diverges. The Fisher–Rao action S\[Γ\] is infinite for every path the paper considers.

## Correction 3: exponential dominance fails, and the effect is not arithmetic

The v31 table is correct, but its reading is not: μ₁/μ₀ grows slowly or not at all, never exponentially, so α(N) → 0. v33 reproduced the v31 table to two decimals (log₁₀ μ₀ = −97.58 and α = 0.1427 at N = 45, λ = 3), then extended it to N = 200 and ran the archive's own test P2 with three control point sets of the same density.

&#91;embedded content: Computed in v33 with Arb (python-flint 0.9), 2.6N + 80 digits, first 500 zeros or 500 control points, λ = 3\]

The ratio log₁₀(μ₁/μ₀) rises only from 6.42 at N = 45 to 7.81 at N = 200, so α(N) falls from 0.143 to 0.039. At λ = 5 the ratio does not grow at all: 6.61 at N = 80, 6.65 at N = 120 and 6.58 at N = 160. The quoted α ≈ 0.14 is 6.4/45, an artefact of stopping at N = 45. Conjecture ED and route W2 are falsified.

Test P2 has a clear answer: GUE, Poisson and evenly spaced points give α = 0.0609, 0.0611 and 0.0608 at N = 120, against 0.0608 for the zeros. The tiny eigenvalues (μ₀ ≈ 10⁻²·²ᴺ) come from the basis, not from the points: on an interval of length 2 log 3 ≈ 2.2 < 2π, the cosines cos(kt) are nearly linearly dependent. Any point set of this density gives the same spectrum. Prediction P3 (a different convergence rate for zeta and GUE) becomes moot, since both limits are 0.

One further problem: the matrix M = 2ΦΦᵀ is a sum over the zeros themselves, taking every ordinate γ real. It is positive semidefinite by construction for any real point set, so it presupposes RH and cannot test it. The New result below replaces it with the prime-side Weil form.

## Correction 4: the operator claims, and withdrawal of the proof chain

The operator-theoretic core of FDBC v2 cannot hold as stated, because the Bost–Connes Hamiltonian has the wrong spectrum. The archive itself records that H\_BC acts with spectrum {log n : n ≥ 1}, which is why Tr e⁻ᵗᴴ = ζ(t) (CMP v27, Corollary 3.3; FDBC\_v2\_B\_Prime\_Resolvent).

**Assumption A.** The corrected form (FDBC\_v2\_NC2\_Assumption\_A\_CORRECTED) asks H\_BC to have eigenvalues exactly {1/4 + γₙ²}. In the representation where its spectrum is {log n}, the regularised determinant vanishes where s(1 − s) = log n:

```latex
s = \tfrac12 \pm i\sqrt{\log n - \tfrac14}, \qquad n \ge 2
```

The first such point sits at height √(log 2 − 1/4) ≈ 0.666, not at 14.1347. In the GNS representation of the KMS₁ state the factor is of type III₁, so the modular operator's spectrum is the whole half-line \[0, ∞) and has no discrete part to match. Assumption A is therefore withdrawn, not merely open.

**Zeros as resolvent poles (CMP v27, §8).** The resolvent (H\_BC − z)⁻¹ has its poles at z = log n, all real. The zeros of ζ(t) = Tr e⁻ᵗᴴ are zeros of a Laplace transform of the spectral measure, not spectral values. The claim is withdrawn.

**Assumption B (spectral rigidity).** H\_mod is self-adjoint, so its spectrum is real; the statement that "the real part of any spectral value is 1/2" does not apply to it. The reflection s ↔ 1 − s of ζ comes from the archimedean Γ-factor (Poisson summation), which the Bost–Connes partition function ζ(β) does not contain; the archive gives no derivation of it from the KMS condition. The claimed proof is withdrawn.

**The CMP v25–v27 papers.** Their abstracts announce "a complete proof chain for the Riemann Hypothesis". The Krein–Rutman step was already retracted in v29, and the steps above fail independently. These files stay in the record for history but are marked superseded and withdrawn. If submission CIMP-D-26-00460 is still with Communications in Mathematical Physics, the editor should be told.

Files marked superseded in v33: CMP\_Paper\_FDBC\_v25, v26, v27; FDBC\_v2\_CMP\_Submission; FDBC\_v2\_Complete\_Proof\_Chain; FDBC\_v2\_Proof\_Chain\_CORRECTED; FDBC\_v2\_Assumption\_B\_PROVED and \_v2; FDBC\_v2\_Spectral\_Gap\_PROVED; FDBC\_v2\_No\_Gap\_Final; Information\_Geometric\_Confinement\_Riemann\_Zeros\_Blum\_2026; blum-claude-weil-minimiser-2026 (its table stands, its conclusions do not).

## New result: the Weil form built from primes recovers the zeros

Computed from the prime powers 2, 3, 4, 5, 7 and 8 and the archimedean place alone, the Weil quadratic form in the archive's cosine basis has a minimiser whose Fourier transform vanishes at the first 17 zeta zeros. No zero of ζ enters the computation.

**Construction.** Take gₖ(t) = cos(kt) on \[−L, L\] with L = log λ and λ = 3, for k = 0, …, N. For each pair, h = gⱼ ∗ gₖ is supported in \[−2L, 2L\], so only prime powers n ≤ λ² = 9 enter. Wⱼₖ is the right-hand side of Weil's explicit formula for that h:

```latex
W_{jk} = \hat h(\tfrac{i}{2}) + \hat h(-\tfrac{i}{2}) - h(0)\log\pi - \gamma_E\, h(0) + \int_0^{\infty} \frac{e^{-t}h(0) - e^{-t/4}h(t/2)}{1-e^{-t}}\,dt - 2\sum_{n \le 9} \frac{\Lambda(n)}{\sqrt n}\, h(\log n)
```

The archimedean term is written in time-domain form. v33 checked this implementation on a Gaussian test function: the prime side and the sum over the first 500 zeros agree to 25 digits (0.2693081132190328436967269).

**Checks against the archive's matrix.** At N = 10, W − M is positive semidefinite (smallest eigenvalue 1.6 × 10⁻⁵⁷, largest 0.027), as expected if M is the 500-zero truncation of W. W is positive definite at N = 25, 45 and 80. Exponential dominance fails for W as well: log₁₀(μ₁/μ₀) = 5.12, 6.38 and 6.87.

&#91;embedded content: Computed in v33: prime-side Weil matrix W, λ = 3, 2.6N + 80 digits; zeros located by bisection and refined, compared with Arb values of the zeta zeros\]

The agreement is best at the first zero and falls off with height. Going from N = 25 to N = 45 adds about 9 digits; going on to N = 80 adds only 1 to 2 more, so at λ = 3 accuracy is now limited by the primes available, not the basis. At N = 25 the transform also has one spurious zero, at 45.99, between the 8th and 9th zeta zeros.

**Why it works, and what it does not show.** If RH holds, W equals the sum of ĝⱼ(γ) ĝₖ(γ) over all zeros. Its minimiser makes that sum tiny, so its transform must nearly vanish at each low zero. This is the phenomenon studied by Connes and Consani (ζ-cycles) and by Connes, Consani and Moscovici (prolate wave operators); v33 reproduces it in a cosine basis. It is a numerical observation, not a proof: positivity at one λ says nothing about RH.

## New direction: scale-resolved Weil positivity

The block is now precise, and it is not a missing idea but a missing scale. Weil positivity at scale λ (test functions supported in \[1/λ, λ\]) is a finite statement about the primes up to λ². v33 shows, by computation and by a short derivation, exactly what such a statement can and cannot say about a zero off the line: everything below the height 2πλ², almost nothing above it. RH is verified up to height 3 × 10¹², so new information about zeros needs λ ≳ 7 × 10⁵, far beyond any certifiable computation. Positivity at small λ remains a real theorem target, and nothing beyond λ = √2 is proved.

**The object.** For λ > 1 and L = log λ, let W\_λ be the prime-side Weil form of the New result, now viewed as an operator on all of L²\[−L, L\] (even functions). Primes enter one at a time as λ² crosses a prime power. Weil's criterion: RH ⇔ W\_λ ≥ 0 for every λ. Known unconditionally: W\_λ ≥ 0 for λ ≤ √2, where no prime contributes (Connes and Consani 2021, Theorem 1, with two vanishing conditions on the test function). To our knowledge no λ > √2 is proved.

**The detection principle.** Write Φ(r) = ∫ f(t) e^{irt} dt for the transform of a real even f supported in \[−L, L\]. A zero on the line at height γ contributes 2Φ(γ)² ≥ 0 to W\_λ(f). A hypothetical zero at 1/2 + δ + iγ (with its three partners) contributes instead

```latex
4\,\mathrm{Re}\,\Phi(\gamma - i\delta)^{2} = 4\big[(\mathrm{Re}\,\Phi)^{2} - (\mathrm{Im}\,\Phi)^{2}\big], \qquad \mathrm{Im}\,\Phi(\gamma - i\delta) = \int f(t)\,\sinh(\delta t)\,\sin(\gamma t)\,dt
```

which is negative when the imaginary part dominates. For small δ it equals 4Φ(γ)² − 4δ²Φ′(γ)² + O(δ⁴). So W\_λ detects the zero exactly when some f can vanish at the on-line zeros near γ while keeping Φ′(γ) large. A function of exponential type L can do that only where the zeros are sparser than its Nyquist density L/π, that is below the height H(λ) = 2πλ², where the zero density (1/2π) log(γ/2π) first reaches L/π.

v33 tested this by taking the first 2000–3000 true zeros, moving the one nearest a chosen height off the line by δ, and computing the smallest eigenvalue of the form in the cosine basis (negative means positivity fails). The results follow the prediction sharply.

&#91;embedded content: Computed in v33 (detect.py): first 2000–3000 zeros, cosine basis with N = 200–250, 100 digits; one zero moved to 1/2 + 0.3 + iγ\]

At λ = 5 (H = 157) a zero at height 99 moved by only δ = 0.02 already makes the form indefinite, with smallest eigenvalue −3.6 × 10⁻³, scaling as δ² and stable under changing the number of zeros (800 to 3000) and the basis size (120 to 200). At λ = 3 (H = 56) the same zero is invisible to within 10⁻³⁰. Above H(λ) the only remaining mechanism is exponential: the negative term can reach ≈ (∫₀ᴸ φ sinh δt)² against a positive on-line contribution ≈ (log γ) ∫₀ᴸ φ², so detection needs roughly e^{2δL} ≳ 2δ log γ. This gains only like (log log γ)/(2 log λ), which is why finite λ never approaches RH except through H(λ) → ∞.

**The certification ladder.** Proving W\_λ ≥ 0 for a given λ > √2 is a finite problem. On the Fourier side the archimedean part behaves like log(|r|/2π), while the prime part is a sum of shifts by log n with weights Λ(n)/√n, of total norm at most 2S(λ), S(λ) = Σ\_{n≤λ²} Λ(n)/√n. So the form is automatically positive on frequencies |r| > 2π e^{2S(λ)}, and positivity reduces to a certified finite block of size K(λ) ≈ (2L/π) · 2π e^{2S(λ)} in a prolate basis, plus a rigorous lower bound on the archimedean operator outside that block (prolate theory). The cost is the number in the last column.

| λ | Primes entering | H(λ) = 2πλ² | Zeros below H | S(λ) | Block size K(λ) |
| --- | --- | --- | --- | --- | --- |
| √2 | none | 12.6 | 0 | 0 | 1 (proved, Connes–Consani) |
| 2 | 2, 3 (and 4) | 25.1 | 2 | 1.47 | ≈ 50 |
| 3 | up to 9 | 56.5 | 12 | 3.54 | ≈ 5 × 10³ |
| 4 | up to 16 | 100.5 | 29 | 5.15 | ≈ 2 × 10⁵ |
| 5 | up to 25 | 157 | 56 | 7.48 | ≈ 2 × 10⁷ |
| 10 | up to 100 | 628 | 361 | 16.9 | ≈ 4 × 10¹⁵ |

λ = 2 is a few days of work and would be the first positivity theorem with primes included. λ = 3 is feasible. λ = 4 is not, with this bound. The wall is S(λ) ≈ 2λ: the archimedean part grows only logarithmically in frequency while the prime mass grows linearly in λ, and short Dirichlet polynomials Σ Λ(n) n^{−1/2−ir} do reach their maximum at some frequencies below e^{cλ}. Any certification that does not use the cancellation in that polynomial frequency by frequency pays e^{2S(λ)}.

**What this means for the goal.** Three statements, each checkable:

1. Positivity at scale λ sees every off-line zero below H(λ) = 2πλ² and essentially none above. RH to height 3 × 10¹² makes scales below λ ≈ 7 × 10⁵ uninformative about zeros; such scales are 10¹² times beyond the certification wall.
2. The only way finite scales reach RH is H(λ) → ∞, and the only way to certify positivity at large λ is to use the cancellation in Σ\_{n≤λ²} Λ(n) n^{−1/2−ir} at every frequency up to e^{cλ}. That cancellation is a statement about zeros up to that height. The circle closes here, not in any particular formalism.
3. A new input must therefore be a positive structure for W\_λ that does not pass through frequency-by-frequency control: a factorisation W\_λ = T\_λ\* T\_λ + (explicitly controlled remainder) in which the prime shifts are absorbed. This is what the Connes–Consani semi-local programme seeks for the places {∞, 2, 3, …}; nothing in this archive supplies it, and v33 makes no claim to.

**Next computations, in order.** (a) Certify W\_λ ≥ 0 at λ = 2 with ball arithmetic: a 50 × 50 block in a prolate basis plus the archimedean tail bound. (b) Measure the detection threshold δ\*(λ, γ) above H(λ) and compare with e^{2δL} ≈ 2δ log γ. (c) Extend (a) to λ = 3. Each is a theorem or a measurement; none is RH.

## Where it blocks, measured, and three new objects

The block is now a number, not a description: the positivity margin at scale λ is super-exponentially small, so Weil positivity is a rigid constraint rather than an open condition, and no perturbative or inductive route can cross it. Three objects follow from that measurement. Each is new relative to the archive, each is defined from primes alone, and each gives a definite computation.

**1. The margin law and the rigidity of the primes.** Define m(λ) = inf W\_λ(f)/‖f‖² over even f supported in \[−L, L\]. v33 computed it from the prime side (Galerkin, so these are upper bounds), together with the first-order tolerance εₙ on each prime weight: replacing Λ(n) by (1 + ε)Λ(n) makes the form indefinite as soon as |ε| > εₙ. Galerkin indefiniteness implies true indefiniteness, so the tolerances are rigorous upper bounds.

| λ | Basis size N | Margin m(λ), upper bound | Tolerance on Λ(2) | Tolerance on the last prime power present |
| --- | --- | --- | --- | --- |
| 2 | 45 (converged) | 10⁻¹²·¹ | 10⁻¹⁰·⁶ | n = 3: 10⁻⁷·⁵ |
| 3 | 45 (N = 30 gave 10⁻³¹) | ≤ 10⁻³⁶·⁶ | ≤ 10⁻³⁵·³ | n = 8: ≤ 10⁻¹²·³ |
| 4 | 30 | ≤ 10⁻⁴⁰·⁸ | ≤ 10⁻³⁹·⁸ | n = 13: ≤ 10⁻¹⁶·⁰ |
| 5 | 60 (N = 30 gave 10⁻⁴⁸) | ≤ 10⁻⁷⁴·² | ≤ 10⁻⁷³·¹ | n = 23: ≤ 10⁻¹⁷·⁵ |

At λ = 2 the margin is converged in N; at λ ≥ 3 it keeps falling as the basis grows, so the true margins are smaller than shown. The prime most recently entered is the least constrained, since the minimiser's autocorrelation at a lag near 2L is small; the prime 2 is the most constrained. All values are from the prime side (margin.py); none uses a zero of ζ.

The margin is the leakage of the minimiser onto the zeros above H(λ): the minimiser's zeros are Gauss-quadrature nodes (object 2), its transform vanishes at the resolved zeros to tens of digits, and what is left is Σ over unresolved zeros of |Φ(γ)|². So m(λ) is controlled by exactly the zeros that scale λ cannot see. The consequence is quantitative and final for inductive strategies, including the archive's exponential-dominance route: positivity at the next scale is not positivity at this scale plus a small perturbation. The next prime enters with weight Λ(n)/√n ≈ 1/λ against a slack below 10⁻³⁶ at λ = 3 and below 10⁻⁷⁴ at λ = 5, and the weight of the prime 2 is pinned to 35 digits by positivity at λ = 3 and to 73 digits at λ = 5.

**2. Pseudo-zero measures: positivity as a moment problem.** By Krein's extension theorem, W\_λ ≥ 0 holds if and only if there is a positive even measure ν on ℝ whose Fourier transform agrees with the prime kernel K\_λ (archimedean kernel plus the impulses −Λ(n)n^{−1/2} at ±log n) on \[−2L, 2L\]. Let E\_λ be the set of such measures. It is convex, it depends on the primes up to λ² only, and RH is the statement that E\_λ is non-empty for every λ; under RH the zero measure Σ δ\_γ lies in all of them. The Galerkin minimiser is the orthogonal function of this moment problem and its zeros are the Gauss-quadrature nodes, which is why they sit on the true zeros: the nodes of a quadrature for a discrete measure converge to its atoms, and H(λ) is where the degrees of freedom run out. The new computable object is the canonical (Krein, maximum-entropy) extension ν\_λ, a positive measure on ℝ built from primes ≤ λ², which under RH converges to the zero measure. Its atoms below H(λ) are the pseudo-zeros already observed; its behaviour above H(λ) has not been looked at.

**3. The prime-driven Krein system: induction on primes done exactly.** Krein's continuous analogue of the Schur algorithm turns a positive-definite kernel on \[0, 2L\] into a first-order system whose coefficient is the kernel itself. Here the coefficient is the archimedean kernel between prime powers and an impulse of mass Λ(n)/√n at each t = log n. Positivity on \[0, 2L\] is regularity of the system up to time 2L; each prime power is one Schur step with a reflection coefficient sₙ; RH is the statement |sₙ| < 1 for every n, with the archimedean flow in between never blowing up. This is the induction on primes that the archive tried to do with eigenvalue ratios, done exactly: the state carried from one prime to the next is the whole Weyl function, not one number. The first questions are computational: the values of sₙ for n ≤ 25, whether 1 − |sₙ| tracks m(λ), and whether the archimedean flow alone fails at some t\* > log 2, which would mean that the prime impulses are needed to keep the system regular.

**4. Certification with the compressed prime mass.** The wall in the previous section used the crude norm 2S(λ). The relevant quantity is the norm S\_eff(λ) of the prime shift operator compressed to \[−L, L\], which v33 computed:

| λ | S(λ) | S\_eff(λ) | Ratio | Block size with S | Block size with S\_eff |
| --- | --- | --- | --- | --- | --- |
| 2 | 1.12 | 0.45 | 0.40 | 26 | 7 |
| 3 | 3.17 | 1.14 | 0.36 | 2.5 × 10³ | 43 |
| 4 | 4.97 | 1.82 | 0.37 | 1.2 × 10⁵ | 210 |
| 5 | 7.16 | 2.48 | 0.35 | 1.1 × 10⁷ | 910 |
| 10 | 16.9 | 5.54 | 0.33 | 4 × 10¹⁵ | 6 × 10⁵ |

The certified-positivity frontier moves from λ = 3 to λ ≈ 5. What is needed to make it rigorous: a certified bound on the compressed norm (a Boas–Kac inequality for positive-definite functions of compact support gives one), the archimedean tail bound from prolate theory, and arithmetic at the precision of the margin, 10⁻⁵⁰ at λ = 5, which ball arithmetic provides.

**Where a proof would have to come from.** In these terms RH is an a priori lower bound m(λ) > 0 for all λ. From the zero side, m(λ) is a sampling constant of the zero set for functions of exponential type L: it needs the on-line zeros above H(λ) to sample such functions, and it fails only if a window of width about 1/L at some height holds a cluster of off-line zeros that outweighs the on-line ones there. Zero-density estimates control off-line zeros on average over long ranges, not in windows of width 1/L. That is the exact gap, and it is the same gap seen from the prime side: a structural reason why the von Mangoldt impulses never drive a reflection coefficient to modulus 1. Connes and Consani supply that reason for the archimedean part alone, through the compression of the scaling action to Sonin space. The open problem is the same factorisation with the shifts included. v33 does not solve it; it gives three handles on it that did not exist in this archive: ν\_λ, sₙ and m(λ).

## The range law: what the primes up to λ² know about the zeros

The primes up to λ² determine the zeros up to H(λ) = 2πλ² to the working precision, and beyond H their information fades at a fixed rate, about one decimal digit per four units of height, reaching nothing near 2H. The sub-function of scale λ is therefore exact on its range and partial above it, and the only correction that makes the next range exact is the next range of primes.

v33 measured this by locating every real zero of the minimiser transform up to the basis limit and matching it with the true zeros (the maximum-entropy extension of the prime kernel has spectral density proportional to 1/|Ψ|², so these zeros are its atoms).

&#91;embedded content: Computed in v33 (recon.py): prime-side W, cosine basis with N = 60 (λ = 2) and N = 160 (λ = 3), zeros of the minimiser transform located on (0, N\] and compared with Arb values of the zeta zeros\]

- λ = 2 (primes 2 and 3 only, H = 25): the first zero to 10.5 digits, the third (at 25.0) to 6.2, then one digit lost per zero; nothing reliable above 45.
- λ = 3 (primes up to 9, H = 56.5): the first zero to 33.6 digits, the twelfth (at 56.4, the resolution height) to 15.8, the 29th (at 98.8) to 2.4; above 110 the zeros found no longer match. A straight-line fit gives digits ≈ 30.5 − 0.25γ, crossing zero at 123 = 2.2 H.
- λ = 5 (primes up to 25, H = 157): all 38 zeros below 120 recovered, none spurious, none missed; the first zero to 104.7 digits, the 13th (at 59.3) to 77, the 38th (at 118.8) to 43.8, against 120-digit reference zeros. The basis limit (N = 120) sits below H here, so the fall-off near 120 is partly the edge of the basis.

Where the basis reaches well past H (λ = 2 and 3) the slope is the same, 0.23 and 0.25 digits per unit of height, and the information ends near 2H. The accuracy at the first zero grows like H: 10.5, 33.6 and 104.7 digits at λ = 2, 3 and 5, about 0.4 to 0.7 times H(λ), that is an error of order exp(−cλ²). Two consequences. First, a range-by-range correspondence: the zeros in (2πλ², 2πλ′²\] are fixed by the primes in (λ², λ′²\], and positivity at scale λ′ is exactly the consistency check between the two ranges. Second, the count matches sampling theory: the number of zeros below H, about 2λ² log λ, equals the number of degrees of freedom of a function supported in \[−L, L\] seen at the zero density, while the number of primes used, about λ²/(2 log λ), is smaller by a factor 4 (log λ)². The information is in the digits of log p, not in the number of primes.

This is reconstruction, not detection: a displaced zero is seen only below H (previous section), but the positions of on-line zeros are recovered well beyond H. Nothing here bears on whether zeros above the reconstruction horizon are on the line.

## The scale-λ Riemann equation

The minimiser is not an arbitrary function: it is Riemann's Φ cut to the window \[−L, L\] by a smooth taper, and its transform is the Riemann Ξ-function times a positive envelope. This turns the range law into an equation. Recall Riemann's representation

```latex
\Xi(r) = \xi(\tfrac12 + ir) = \int_{-\infty}^{\infty} \Phi(u)\, e^{iru}\, du, \qquad \Phi(u) = \sum_{n \ge 1} \big(2\pi^{2} n^{4} e^{9u/2} - 3\pi n^{2} e^{5u/2}\big)\, e^{-\pi n^{2} e^{2u}}
```

Φ is built from the integers through the theta function and decays doubly exponentially: Φ(1.1) ≈ 1.4 × 10⁻⁹ while Φ(0) = 0.45. v33 compared the prime-side minimiser f\_λ(t) with Φ on \[0, L\] at λ = 3 (N = 45), normalised at t = 0:

| t | 0 | 0.18 | 0.37 | 0.55 | 0.73 | 0.92 | 1.01 | 1.10 = L |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| f\_λ(t) / (c Φ(t)) | 1 | 0.975 | 0.871 | 0.621 | 0.241 | 0.014 | 0.0004 | 10⁻⁹ |

So f\_λ = c Φ w\_λ with w\_λ a smooth even taper, 1 near the origin and vanishing at the edge of the window. On the Fourier side, G\_λ = Ψ\_λ/Ξ is positive and smooth all the way to the horizon: at λ = 3 (N = 60) its logarithm rises slowly and regularly from −47.6 at r = 0 to −44.2 at r = 57 ≈ H, with no sign change, and the first sign changes appear at r = 61, 62, 65, where the reconstruction stops. The envelope is close to log-quadratic, log₁₀ G ≈ −47.6 + 0.001 r².

**The equation.** The minimiser solves the Euler–Lagrange equation of the prime-side form, (K\_λ ∗ f)(t) = m(λ) f(t) on \[−L, L\], with K\_λ the archimedean kernel plus the impulses at ±log n. Substituting f = Φ w gives an equation for the taper alone:

```latex
\big(K_\lambda * (\Phi\, w_\lambda)\big)(t) = m(\lambda)\, \Phi(t)\, w_\lambda(t), \qquad |t| \le \log \lambda
```

with m(λ) the margin of the previous section, below 10⁻³⁶ at λ = 3. Read as a statement about arithmetic: the primes up to λ², through K\_λ, reproduce the theta-function side Φ on the window \[−log λ, log λ\], and the error is governed by the tail of Φ beyond the window, which is of order exp(−πλ²) = exp(−H/2). The measured accuracy at the first zero, exp(−(1.0 to 1.5) H), is of that form. The scale λ plays the role the archive gave to time: as it grows, the window grows, the horizon H(λ) = 2πλ² moves out, and everything inside is frozen to precision exp(−cH).

**What the equation does and does not say.** It makes the reconstruction law a theorem target: under RH it is a statement in sampling theory about how well a function supported in \[−L, L\] can vanish on the zero set, and unconditionally it describes the pseudo-zero measure. It gives the window w\_λ a definition independent of any zero, so its shape (close to a prolate taper) can be studied analytically. It does not bear on RH: the zeros of a windowed Φ being real says nothing about the zeros of Φ itself, since the Pólya and de Bruijn multiplier theorems run in the opposite direction. The first analytic question it poses is whether w\_λ converges, as λ → ∞, to a universal profile in t/L.

## Two tests of the equation

The window has a near-universal core and a scale-dependent edge, and the pseudo-zeros respond to a change in any prime through a single mode. Neither result was predicted; both are stated here as measured.

**Universality of the window.** w\_λ(t) = f\_λ(t)/(cΦ(t)) at the three scales, as a function of s = t/L, normalised to 1 at the origin (basis sizes 45, 45 and 120; the λ = 5 window at N = 60 was not converged and is discarded):

| s = t/L | 0.30 | 0.42 | 0.54 | 0.60 | 0.66 | 0.72 | 0.78 | 0.84 | 0.90 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| λ = 2 | 0.919 | 0.821 | 0.668 | 0.570 | 0.461 | 0.345 | 0.232 | 0.132 | 0.056 |
| λ = 3 | 0.901 | 0.762 | 0.535 | 0.396 | 0.256 | 0.134 | 0.051 | 0.012 | 0.001 |
| λ = 5 | 0.931 | 0.806 | 0.558 | 0.392 | 0.223 | 0.090 | 0.020 | 0.0015 | 1.2 × 10⁻⁵ |

At λ = 3 and λ = 5 the windows agree to within 0.02 up to s = 0.6 and then separate: the edge closes earlier, relative to the support, as λ grows. λ = 2 is gentler throughout. Of three one-parameter shapes, a Gaussian in s, a prolate ground function and the bump exp(−βs²/(1 − s²)), only the bump fits (β = 0.80, 1.76, 2.59; rms of the log 0.12, 0.13, 0.30). So the scaling t/L is not the right one for the edge, and the universal profile, if there is one, needs λ = 4, 7 and 10 to be seen. On the Fourier side G\_λ = Ψ\_λ/Ξ is positive and smooth up to the horizon at every scale, but its curvature is not universal either (convex at λ = 3, nearly flat at λ = 5 in the resolved range).

**Response of the pseudo-zeros to the primes.** v33 computed the Jacobian Sⱼₙ = ∂γⱼ(λ)/∂εₙ, the motion of the j-th reconstructed zero when Λ(n) is replaced by (1 + ε)Λ(n), by finite differences in the linear regime (ε = 10⁻⁴⁵ and 10⁻⁴⁴ agree to 10⁻¹⁰) at λ = 3, N = 45, for the eight zeros below 45. The expected signature of the explicit formula, a factor cos(γⱼ log n), is absent at first order. Instead Sⱼₙ is of rank one to a good approximation: every zero moves in the same direction, with an amplitude A(γⱼ) that grows by about a factor 100 per zero (4 × 10² at γ₁ to 3 × 10¹⁶ at γ₈ for the prime 2, in units of relative change of Λ) and a prime factor B(n) that falls steeply with the lag: 1 for n = 2, 10⁻³ for n = 4, 10⁻⁶·⁵ for n = 5, 10⁻¹⁵ for n = 7, 10⁻²² for n = 8. The reason is first-order perturbation theory: the correction to the minimiser is dominated by the second eigenfunction e₁, whose eigenvalue is the closest, so the response direction is e₁ for every prime, and A(γⱼ) ≈ −Ψ₁(γⱼ)/Ψ₀′(γⱼ), large where the zero is least accurately pinned. The caveat is that the spread of the Galerkin eigenvalues reflects the conditioning of the cosine basis; whether the rank-one response survives in the L²-normalised problem as N → ∞ is the next check.

What these two tests change: the equation f\_λ = cΦw\_λ stands at every scale computed, the shape of w\_λ is not yet understood, and the dynamics of the pseudo-zeros under the addition of primes is, to first order, a two-mode problem (e₀ fixes the zeros, e₁ carries their motion). That reduction is the one new structural fact of this section, and it is the natural starting point for the Krein-system computation of the previous section.

## The way: a staged path, with its first rigorous step taken

The path has three stages. The first is a theorem within reach and v33 has completed its finite part rigorously; the second is a theorem that needs one new inequality; the third is the open heart of the problem, now stated as a per-prime criterion that can be computed.

**Stage A, the first positivity theorem with primes: W\_λ ≥ 0 at λ = 2.** The proof has two parts. The finite part is a ball-arithmetic certificate on a block of test functions; the tail part is an operator bound outside the block.

The finite part is done. v33 rebuilt the prime-side form at λ = 2 (primes 2 and 3, archimedean place, poles) with every entry an Arb ball enclosing the true value: closed forms for the pole, prime and autocorrelation terms, rigorous quadrature for the archimedean integral, and an explicit bound for the removable singularity at the origin. Positivity of the block was then certified by Gershgorin's theorem applied, in ball arithmetic, to QᵀWQ with Q an approximate eigenbasis whose determinant is enclosed away from zero.

**Certified (v33, certify.py).** For λ = 2, the prime-side Weil form is positive definite on span{cos(kt) : 0 ≤ k ≤ 30} ⊂ L²\[−log 2, log 2\]. Every matrix entry is enclosed in a ball of radius below 10⁻¹¹⁸; the smallest eigenvalue of the block is 7.7 × 10⁻⁹⁰ in coefficient coordinates (the cosine basis is badly conditioned on this short interval; in L² terms the margin is about 10⁻¹²), and the Gershgorin lower bound after conjugation by the approximate eigenbasis is positive with room to spare, 7.7 × 10⁻⁹⁰ ± 10⁻⁹⁶. The same certificate holds for the block k ≤ 20 (smallest eigenvalue 9.3 × 10⁻⁶⁰). Two inputs are not themselves machine-checked: the closed forms for the matrix entries, which were validated against the zero side to 10⁻⁵⁷, and one crude derivative bound used at the origin of the archimedean integral.

**The tail part is now done, by a different route than planned, and the result is Theorem 2 of the paper: W\_2(f) ≥ 1.38 × 10⁻¹³ ‖f‖² for every even f supported in \[−log 2, log 2\].** The prolate plan as written above does not work: the prime shift operator, computed in the prolate basis (prolate\_cross.py, Ω = 80, 90 even modes), couples the block to its complement with norm 0.8 at every cut, through the modes of the Landau–Widom transition region, which are neither in-band nor out-of-band; only the fully in-band modes decouple. What works is to put the primes inside the symbol before splitting: in frequency W\_λ is the integral of |f̂|² against σ\_λ(r) = Re ψ(¼ + ir/2) − log π − 2Σ Λ(n)n^{−1/2} cos(r log n), plus the pole term, and σ\_λ ≥ m\_out = Re ψ(¼ + iΩ/2) − log π − 2S for |r| ≥ Ω (positive for Ω > 2πe^{2S} = 59.6 at λ = 2). Replacing σ by m\_out outside \[−Ω, Ω\] gives a smaller form B whose every term except m\_out‖f‖² factors through Slepian's band-limiting operator, so in the Legendre basis the entries decay like the in-band energies, which fall below 10⁻⁴² past degree 110. The 56 × 56 block was certified in Arb (1596 rigorous integrals, 9 minutes): λ\_min ≥ 1.3805 × 10⁻¹³, Schur correction 5 × 10⁻³⁸. The exact margin is thus enclosed: 1.38 × 10⁻¹³ ≤ m(2) ≤ 7.6 × 10⁻¹³. On the same day this was found, a literature search turned up that the identical reduction had been published in August 2026 by Xuefeng Zhu ([arXiv:2608.24827](https://arxiv.org/abs/2608.24827)), who certifies positivity for all support up to L = 0.8 (ours is L = log 2 = 0.693) with the lower bound 8.9 × 10⁻¹⁸, gives certified upper bounds on the margin up to L = 2, the same resolution height 2πe^{2L} = H(λ), and a barrier theorem showing that this kind of certificate cannot pass L ≈ 3.2. Our Theorem 2 is therefore an independent confirmation at L = log 2 with a different implementation, not a first. What is ours and not in Zhu: Theorem 3, a certified rigidity statement — the form becomes indefinite when Λ(2) is raised by a relative 2.6 × 10⁻¹¹ or lowered by 1.8 × 10⁻⁷, or when log 2 is lowered by 2.3 × 10⁻¹² or raised by 2.5 × 10⁻⁸ (eight ball-arithmetic certificates, certify\_rigidity.py); and, under RH, Theorem 4: the form is coercive at every scale by Beurling's sampling theorem, and a small margin localises the minimiser's zeros.

**Stage B, λ = 3 and beyond: blocked by Zhu's barrier.** The compressed norm S\_eff(λ) (0.36 S at λ = 3) cannot replace S in the envelope: the certificate must bound the prime comb pointwise in frequency, where it attains 2S, and the in-band cancellation on the low block, on which the 10⁻¹² margin lives, is destroyed by any operator-norm bound. The next prime, 5, enters at 2L = 1.609 and lifts the cutoff from 119 to 503; λ = 3 (six prime powers, cutoff 3570, margin 10⁻³⁷) needs a Legendre block of several thousand modes. Zhu's Theorem 1.4 shows the cutoff 2πe^{2S} is optimal within such envelopes, so the method stops near L ≈ 3. S\_eff remains the right constant for the operator arguments of Stage C and of the frame identity, not for Stage B.

**Stage C, the heart: a criterion per prime.** In the Krein system of the previous sections, the prime p enters at time log p as an impulse of mass Λ(p)/√p, and the system stays regular, that is positivity survives, exactly when the corresponding reflection coefficient has modulus below 1. Unwinding the Schur step, that condition is an inequality between the weight of the prime and the Christoffel function of the system built from everything before it:

```latex
\frac{\Lambda(p)}{\sqrt p} \;<\; \frac{1}{\kappa_{\lambda}(\log p)}, \qquad \lambda^{2} = p^{-},
```

where κ\_λ(u) is the diagonal of the reproducing kernel of the Weil form at scale λ evaluated at the lag u (the largest value a unit-norm test function's autocorrelation can take at that lag). RH is the statement that this inequality holds at every prime power, with the archimedean flow regular in between. Nothing in it refers to a zero. The margin by which it holds at p is the margin m(λ) of the previous section, which the rigidity measurements already show to be tiny and decreasing. The way, then, is to understand κ\_λ(u): it is a prime-free object (a reproducing-kernel value of a form built from primes below p), and the question becomes why Λ(p)/√p sits inside the interval it allows, at every p, with a margin that closes like exp(−cλ²). That is the same question as the semi-local factorisation of Connes and Consani, but posed one prime at a time, with a quantity on each side that can be computed.

The first computation of Stage C is to evaluate κ\_λ(log p) and Λ(p)/√p for p ≤ 23 and plot the ratio. If the ratio has structure (a smooth function of p, or of log p), there is a conjecture to make about it. If it is erratic, the criterion is true but not explanatory, and the way closes here.

## Context: the OpenAI quasi-Riemann manuscript

The 7/8 result neither uses nor affects anything in this archive. In early October 2026 OpenAI released a manuscript dated 30 September (family 003 of its openai/math collection) claiming that every Dirichlet L-function, including ζ, has no zero with real part above 7/8. A companion gives 11/12 by a second argument, and a third paper excludes Landau–Siegel zeros in the logarithmic range. The 7/8 bound for ζ, Dirichlet and Hecke L-functions over ℚ(√−3) has a Lean formalisation in the repository; outside mathematicians have not yet confirmed the papers.

v33 searched the LaTeX sources of all three papers for every term specific to this programme (Fisher, Bures, KMS, Bost–Connes, Witten, Weil quadratic form, Coulomb, GUE, Li coefficients, Bakry–Émery, Blum, Zenodo). There is no occurrence and no citation. The method is automorphic and analytic: Kubota's metaplectic theory, Patterson's cubic theta series, and sextic large-sieve estimates.

For this programme the result matters in one way. If it is confirmed, any future claim about zeros must be consistent with a zero-free half-plane Re(s) > 7/8. Nothing in v33 depends on it, and the Weil-form computations above are equally valid either way.

## Open problems and next steps

One thread is worth continuing: the prime-side Weil form and its scale ladder (previous section). It is non-circular, it reproduces a phenomenon taken seriously in the literature, and every claim about it can be checked by computation. The first target is a certified proof of W\_λ ≥ 0 at λ = 2. The other routes in the archive should not be revived without a new idea that answers Corrections 1–4.

1. **Accuracy law.** Measure how many zeros the minimiser recovers, and to how many digits, as a function of N and λ (that is, of how many primes enter). v33 has three data points per λ = 3; a grid over λ = 2–20 and N = 20–200 would show whether recovery tracks the number of degrees of freedom, roughly (N + 1) · log λ / π.
2. **Comparison with the prolate picture.** Connes, Consani and Moscovici explain the recovery through prolate spheroidal wave functions. Repeating the computation in a prolate basis instead of cosines would test whether the cosine basis adds or loses anything.
3. **Certified positivity.** With ball arithmetic (Arb), the statement "W restricted to the span of cos(kt), k ≤ N, on \[−log λ, log λ\] is positive definite" can be proved for given N and λ. It does not imply RH, but it is a correct, checkable theorem about primes.
4. **Separate the physics.** Publish the physical CST material (Bures transport, the transmon Witten Laplacian) as its own record, so that its status is judged on its own evidence and not tied to the RH claims withdrawn here.

What v33 does not attempt: any argument that would turn positivity at finite λ into RH. Weil's criterion needs positivity for every λ and every test function, and no step in this archive bridges that gap.

## Reproducibility

Every number in v33 comes from the scripts in the companion archive v33\_code.zip, which should be uploaded with this version. They need Python 3, mpmath, numpy and python-flint (Arb); the full set runs in under an hour on two cores, most of it the N = 80 prime-side case.

| File | What it produces |
| --- | --- |
| make\_points.py | The first 500 zeta zeros (Arb, 55 digits) and the GUE, Poisson and evenly spaced control sets, seed 20261007 |
| checks.py | Corrections 1 and 2 (zero mode of R, Hessian signs, generic Fisher pole) and the explicit-formula check |
| weil\_gram.py | The zero-side matrix M and its two smallest eigenvalues (Correction 3); run as weil\_gram.py zeros500.txt 3 10,20,45 |
| weil\_prime.py | The prime-side Weil matrix W (New result) |
| ccm\_test.py | Eigenvalues of W and the zeros of its minimiser's transform; run as ccm\_test.py 45 3 |
| results/ | The outputs quoted in this document |

detect.py (with zeros3000.txt) produces the detection experiments of the last section; run as detect.py 4 250 2000 150 0.1,0.3,0.49 100 for λ = 4, basis size 250, 2000 zeros, height 150, three values of δ, 100 digits.

margin.py computes the L²-normalised margin and the prime-weight tolerances (run as margin.py 3 45); rigidity.py confirms the tolerances by bisection (rigidity.py 3 25); the compressed prime mass S\_eff comes from s\_eff.py (output in results/s\_eff.txt).

recon.py locates all real zeros of the minimiser transform up to the basis limit and matches them with 120-digit reference zeros (zeros500\_120d.txt); run as recon.py 3 160. Outputs in results/recon\_\*.txt.

factor.py computes G = Ψ/Ξ on a grid (factor.py 3 60); phi\_test.py checks Riemann's Φ against ξ and prints f\_λ/(cΦ) on \[0, L\] (phi\_test.py 3 45).

window.py samples w\_λ(t/L) and fits the three trial shapes (window.py 5 120); jacobian.py computes the response of the reconstructed zeros to each prime weight (jacobian.py 3 45 1e-45).

certify.py builds the λ = 2 block in ball arithmetic and certifies its positivity (certify.py 2 30 200; about 15 minutes). Outputs in results/certify\_l2\_n20.txt and certify\_l2\_n30.txt.

prolate.py computes the prolate eigenvalues λ\_n(c) for c = Ω log 2 (Bouwkamp's tridiagonal Legendre method; even modes fall below 10⁻²⁰ at index 56, 66, 76 for Ω = 80, 100, 120); prolate\_cross.py computes the prime shift operator in the prolate basis and its low/high coupling (results/prolate\_\*.txt). symbol\_B.py is the floating-point exploration of the envelope form B in the Legendre basis; certify\_symbol.py is its Arb certificate (run as certify\_symbol.py 2 80 110 40 28; 9 minutes; results/certify\_symbol\_l2.txt and the cached entries symbol\_cache\_\*.json). certify\_rigidity.py produces the eight indefiniteness certificates of Theorem 3 (run as certify\_rigidity.py 2 24; 3.5 minutes; results/certify\_rigidity\_l2.txt and rigidity\_cert\_l2\_n24.json).

Precision is 2.6N + 80 decimal digits throughout, enough to resolve μ₀ ≈ 10⁻²·²ᴺ. Doubling the number of Gauss–Legendre nodes (192 to 384) for the archimedean integral changes W by less than 10⁻⁸⁰ at N = 10.

## Zenodo record text

The current record's title and keywords still describe the original physics framework ("Cascade Dynamics", "time travel"). The text below replaces them for version v33.

**Title**

Configuration Space Temporality and the FDBC Programme — Research Archive v33: Corrections, Withdrawals and a Non-Circular Weil-Form Test

**Description**

Version 33 is a correction release. A full re-audit of the FDBC/Riemann Hypothesis material, with every number recomputed, withdraws all claims that RH has been proved or reduced to a single open assumption.

Corrections: (1) the Fisher pole identity |I(s)|·|s − ρ|² = 1 holds at a simple zero of any analytic function and carries no information about RH; (2) the rigidity matrix has an exact zero eigenvalue and the transverse Hessian of the log-gas energy is its negative, so the confinement constant is undefined; (3) the exponential dominance of the Weil minimiser is falsified: extended to N = 200, α(N) falls as roughly 1/N, and GUE, Poisson and evenly spaced points give the same values (test P2); (4) Assumption A, Assumption B and the CMP v25–v27 proof chain are withdrawn.

New result: the Weil quadratic form computed from the prime powers up to 9 and the archimedean place alone, in a cosine basis on \[−log 3, log 3\], has a minimiser whose Fourier transform vanishes at the first 17 zeta zeros, the first to 36 digits. This reproduces in a new basis a phenomenon studied by Connes, Consani and Moscovici. It is a numerical observation, not a proof.

A final section, scale-resolved Weil positivity, shows by computation that positivity at scale λ detects an off-line zero only below height 2πλ², and proves, with ball arithmetic, that the truncated Weil form at λ = 2 (primes 2 and 3, test functions supported in \[1/2, 2\]) is positive on the whole window with margin between 1.38 × 10⁻¹³ and 7.6 × 10⁻¹³, and that it becomes indefinite when Λ(2) or log 2 is perturbed at the eleventh digit; the positivity statement was obtained independently and for a larger window by X. Zhu (arXiv:2608.24827, August 2026), which the accompanying paper cites. A last section measures the positivity margin from the prime side (below 10⁻³⁶ at λ = 3), shows that positivity pins the prime weights to tens of digits, and introduces three new objects: pseudo-zero measures, the prime-driven Krein system, and the compressed prime mass. The accompanying paper "Scale-Resolved Weil Positivity" (PDF and source in the paper/ folder) collects these results with the detection law, the range law, the Φ-window structure of the minimiser, and the conditional coercivity theorem. Earlier files remain for the record; those superseded are listed in the v33 document. Code and data: v33\_code.zip (scripts, results and the certificates of the paper's Theorems 1–3).

**Keywords**

Riemann hypothesis · Weil explicit formula · Weil quadratic form · zeta zeros · Bost–Connes system · prolate spheroidal wave functions · high-precision computation · erratum · negative results · configuration space temporality

**Version label:** v33 · **Files to add:** this document as PDF, v33\_code.zip

## Sources

- F. D. Blum et al., [Configuration Space Temporality — Complete Research Archive v32](https://zenodo.org/records/18859602), Zenodo, DOI 10.5281/zenodo.18859602. All 73 files were downloaded and audited for v33.
- A. Connes and C. Consani, [Spectral triples and ζ-cycles](https://ems.press/journals/lem/articles/11033001), L'Enseignement Mathématique 69 (2023), 93–148 ([arXiv:2106.01715](https://arxiv.org/abs/2106.01715v1)).
- A. Connes, C. Consani and H. Moscovici, [Zeta zeros and prolate wave operators](https://arxiv.org/abs/2310.18423v1), arXiv:2310.18423.
- OpenAI, [openai/math](https://github.com/openai/math), family 003: [The Quasi-Riemann Hypothesis: a zero-free half-plane Re s > 7/8](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf) and its [Lean scope note](https://github.com/openai/math/blob/main/lean/docs/003.md).
- J.-B. Bost and A. Connes, Hecke algebras, type III factors and phase transitions with spontaneous symmetry breaking in number theory, Selecta Math. 1 (1995), 411–457 (cited from memory, not re-opened).
- X.-J. Li, The positivity of a sequence of numbers and the Riemann hypothesis, J. Number Theory 65 (1997), 325–333 (cited from memory, not re-opened).
- X. Zhu, [Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau–Widom decay law](https://arxiv.org/abs/2608.24827), arXiv:2608.24827 (v2, 3 September 2026). Opened and read for the paper's Section 7 and 8.
- A. Groskin, [A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form](https://arxiv.org/abs/2607.02828), arXiv:2607.02828 (July 2026); abstract read.
- M. Suzuki, [Weil's quadratic form via the screw function](https://arxiv.org/abs/2606.09096), arXiv:2606.09096 (June 2026); abstract read.
- F. D. Blum, [Scale-Resolved Weil Positivity](https://claude.ai/code/artifact/2c1ec67d-9e12-4244-bc8e-4e68c0f963ca) (the accompanying paper, this session).
