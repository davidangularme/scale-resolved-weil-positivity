# Scale-Resolved Weil Positivity: Reconstruction of Zeta Zeros from Finitely Many Primes, a Detection Law, and a Certified Block

Oct 7, 2026 · @Frédéric David Blum

Frédéric David Blum, Catalyst AI Research, Haifa. Draft for Zenodo / arXiv (math.NT). Computations and drafting assisted by Claude (Anthropic).

## Abstract

Weil's criterion states that the Riemann Hypothesis is equivalent to the positivity of a quadratic form, built from the primes and the archimedean place through the explicit formula, on all test functions. Restricting the test functions to those supported in \[1/λ, λ\] gives a form W\_λ that depends only on the prime powers up to λ²; positivity was known unconditionally for λ ≤ √2, where no prime contributes (Yoshida, Connes–Consani), and since September 2026 for λ ≤ e^{0.8} ≈ 2.23 (Zhu). We study W\_λ for λ between 2 and 5 by high-precision computation and prove one rigorous statement. (i) Detection: a zero at 1/2 + δ + iγ contributes 4 Re Φ(γ − iδ)² to W\_λ(f), and W\_λ becomes indefinite for any δ > 0 when γ lies below the resolution height H(λ) = 2πλ², where the zero density first reaches the Nyquist density of the test functions; above H(λ) detection is exponentially suppressed. (ii) Reconstruction: the minimiser of W\_λ has a Fourier transform whose real zeros coincide with the zeros of ζ below H(λ) to a precision of order exp(−cH(λ)); with the primes up to 25 the first 38 zeros are recovered, the first to 105 digits, and the accuracy decays linearly in height to nothing near 2H(λ). (iii) Structure: the minimiser equals Riemann's function Φ (the theta-function side of ξ) multiplied by a smooth taper on \[−log λ, log λ\], so that its transform is Ξ times a positive envelope up to the horizon. (iv) Rigidity: the positivity margin of W\_λ is 7.6 × 10⁻¹³ at λ = 2 and below 10⁻³⁶, 10⁻⁷⁴ and 10⁻¹³² at λ = 3, 5 and 6, and positivity at λ = 5 pins the weight Λ(2) to 73 digits. (v) Theorems: with ball arithmetic, W\_2(f) ≥ 1.38 × 10⁻¹³ ‖f‖² for every even f supported in \[−log 2, log 2\], which encloses the exact margin between 1.38 × 10⁻¹³ and 7.6 × 10⁻¹³ (an independent certification, at L = log 2, of a theorem of Zhu valid for L ≤ 0.8); and a certified rigidity statement: the same form becomes indefinite when Λ(2) is changed by a relative amount of order 10⁻¹¹ or log 2 is shifted by an amount we certify. Under RH we prove that W\_λ is coercive for every λ, by Beurling's sampling theorem, and that a small margin localises the minimiser's zeros. Nothing here bears on the Riemann Hypothesis beyond the finite scales treated.

## 1. Introduction and statement of results

Let ρ = ½ + i r run over the non-trivial zeros of ζ, and for a real even f supported in \[−L, L\] write Φ\_f(r) = ∫ f(t) e^{irt} dt. Weil's explicit formula expresses the sum over zeros of Φ\_f(ρ)Φ\_f(1 − ρ̄)̄ through the primes and the archimedean place, and Weil's criterion says that RH holds if and only if this quantity is non-negative for every f. If all zeros are on the line it equals Σ\_γ |Φ\_f(γ)|². We call W\_λ the prime side of this identity restricted to f supported in \[−L, L\], L = log λ: it involves only the prime powers n ≤ λ², and RH is equivalent to W\_λ ≥ 0 for all λ. Connes and Consani \[CC21\] proved W\_λ ≥ 0 for λ ≤ √2, the range in which no prime enters, by compressing the scaling action to Sonin's space; to our knowledge no λ > √2 is proved.

This paper is a quantitative study of W\_λ for 2 ≤ λ ≤ 6. All computations are on the prime side; no zero of ζ enters the construction of W\_λ, and the zeros are used only to compare with. Our results are of four kinds: certified theorems at λ = 2 (Theorems 1–3), observations at λ ≤ 6, a conditional theorem under RH (Theorem 4), and a dictionary. Three 2026 preprints on the same object, in particular Zhu's certified positivity for support up to L = 0.8, are discussed in Section 8; they were found after this work was done, and the paper has been revised to cite them and to state precisely what is ours.

**Theorem 1 (certified block positivity).** Let λ = 2. The form W\_2, whose prime part involves only 2 and 3, is positive definite on the subspace span{cos(kt) : 0 ≤ k ≤ 30} of L²\[−log 2, log 2\]. The proof is a ball-arithmetic computation: every matrix entry is enclosed in an interval of radius below 10⁻¹¹⁸ and positivity is certified by Gershgorin's theorem after conjugation by an approximate eigenbasis (Section 7). The same holds for the block k ≤ 20, and at λ = 3 (primes up to 9) for the block k ≤ 45 (Theorem 1′). Theorem 2 (Section 7) completes the λ = 2 case: W\_2 ≥ 1.38 × 10⁻¹³ on all of L²\[−log 2, log 2\], and Theorem 3 certifies the rigidity of the prime data there.

**Proposition 2 (detection).** For a zero ρ = ½ + δ + iγ with δ > 0, the four zeros ρ, ρ̄, 1 − ρ, 1 − ρ̄ contribute exactly 4 Re Φ\_f(γ − iδ)² = 4\[(Re Φ\_f)² − (Im Φ\_f)²\] to the zero side, with Im Φ\_f(γ − iδ) = ∫ f(t) sinh(δt) sin(γt) dt. In particular, to second order in δ, replacing an on-line zero at γ by such a quadruplet changes W\_λ(f) by 2Φ\_f(γ)² − 4δ²Φ\_f′(γ)². Numerically (Section 3), with the true zeros and one zero displaced by δ, W\_λ is indefinite for δ as small as 0.02 when γ < H(λ) = 2πλ², and the indefiniteness decays super-exponentially in γ above H(λ): at λ = 4 and δ = 0.3 it is 10⁻² at γ = 99, 10⁻⁷ at 111, 10⁻¹² at 124, 10⁻³³ at 150 and below 10⁻⁹⁹ at 200. H(λ) is the height at which the zero density (1/2π) log(γ/2π) reaches the Nyquist density L/π of functions of exponential type L.

**Observation 3 (reconstruction, the range law).** Let f\_λ be the minimiser of W\_λ on span{cos(kt) : k ≤ N} and Ψ\_λ its transform. The real zeros of Ψ\_λ coincide with the ordinates of the zeros of ζ: at λ = 5 (primes up to 25, N = 120) all 38 zeros below 120 are recovered, none spurious, the first to 104.7 digits and the 38th to 43.8; at λ = 3 (primes up to 9, N = 160) the accuracy falls linearly from 33.6 digits at the first zero to nothing near 2.2 H(3), at about one digit per four units of height; at λ = 2 (primes 2 and 3) the first zero is obtained to 10.5 digits. The minimiser's zeros are the Gauss-quadrature nodes of the moment problem that W\_λ defines (Section 4).

**Observation 4 (structure).** f\_λ(t) = c Φ(t) w\_λ(t), where Φ is Riemann's function with ξ(½ + ir) = ∫ Φ(u)e^{iru}du and w\_λ is a smooth even taper on \[−L, L\], equal to 1 near the origin and vanishing at the edge; correspondingly Ψ\_λ = Ξ · G\_λ with G\_λ smooth and positive up to H(λ). The windows at λ = 3 and 5 agree to 0.02 up to t/L = 0.6 and differ at the edge (Section 5).

**Observation 5 (margin and rigidity).** The L²-normalised margin m(λ) = inf W\_λ(f)/‖f‖² is 10⁻¹²·¹ at λ = 2 (converged in N), below 10⁻³⁶·⁶ at λ = 3 and below 10⁻⁷⁴·² at λ = 5. Replacing Λ(2) by (1 + ε)Λ(2) makes W\_λ indefinite for |ε| > 10⁻³⁵ at λ = 3 and |ε| > 10⁻⁷³ at λ = 5; these are rigorous upper bounds on the tolerance, since indefiniteness on a subspace implies indefiniteness. The first-order response of the reconstructed zeros to a change in any prime weight is of rank one, carried by the second eigenfunction of W\_λ (Section 6).

**On certification cost (Section 7).** Positivity of W\_λ on frequencies above 2π exp(2S) is automatic from the growth of the archimedean multiplier, S = Σ\_{n≤λ²} Λ(n)/√n, and this is what makes Theorem 2 a finite computation. The norm of the shift operator compressed to \[−L, L\] is only 0.33–0.40 of S, but it cannot replace S in a pointwise envelope (Section 7, Remark 3), and Zhu's barrier theorem shows the cutoff 2πe^{2S} is optimal within such envelopes; the method therefore stops near λ ≈ e³.

What is not claimed. None of this bears on the Riemann Hypothesis beyond the scales treated: positivity at a finite scale sees displaced zeros only below H(λ), RH is verified far above any H(λ) reachable here, and the reconstruction recovers positions of on-line zeros without saying anything about zeros above the reconstruction horizon. Theorem 1 is a finite-dimensional statement; the infinite-dimensional theorem at λ = 2 is set out as a programme, not proved.

## 2. The prime-side Weil form at scale λ

**Conventions.** For a real even h supported in \[−2L, 2L\] let ĥ(r) = ∫ h(x)e^{irx}dx, so that ĥ is real and even on ℝ and entire. The explicit formula in the form we use is

```latex
\sum_{\rho} \hat h\!\left(\tfrac{\rho - 1/2}{i}\right) = \hat h(\tfrac{i}{2}) + \hat h(-\tfrac{i}{2}) - h(0)\log\pi - \gamma_E\,h(0) + \int_0^{\infty} \frac{e^{-t}h(0) - e^{-t/4}h(t/2)}{1-e^{-t}}\,dt - 2\sum_{n \ge 2} \frac{\Lambda(n)}{\sqrt n}\, h(\log n),
```

where the sum on the left runs over the non-trivial zeros ρ and ĥ is evaluated at (ρ − ½)/i, which is the ordinate γ when ρ = ½ + iγ. The archimedean integral is the time-domain form of (1/2π)∫ ĥ(r) Re ψ(¼ + ir/2) dr; we checked the two forms agree to 20 digits, and we checked the whole identity on a Gaussian test function against the sum over the first 500 zeros: both sides equal 0.2693081132190328436967269 to 25 digits.

**The form.** For f real, even, supported in \[−L, L\], put h = f ∗ f̃ (autocorrelation), so that ĥ(r) = |Φ\_f(r)|² for real r and h is supported in \[−2L, 2L\]. Define W\_λ(f) as the right-hand side above with this h. Only n ≤ λ² contribute, and n = λ² itself contributes nothing since h vanishes at the endpoint. If RH holds, W\_λ(f) = Σ\_γ |Φ\_f(γ)|² (each pair ±γ counted); in general it equals Σ\_ρ Φ\_f(ρ′)Φ\_f(1 − ρ̄′)̄ with ρ′ = (ρ − ½)/i.

**Galerkin spaces.** We use the basis g\_k(t) = cos(kt), 0 ≤ k ≤ N, on \[−L, L\]. Its transforms are ĝ\_0(r) = 2 sin(rL)/r and ĝ\_k(r) = sin((k − r)L)/(k − r) + sin((k + r)L)/(k + r), and the autocorrelations h\_{jk} = g\_j ∗ g\_k are elementary trigonometric expressions on \[0, 2L\]. The matrix W\_{jk} = W\_λ(g\_j, g\_k) (polarised) is computed with: the pole term 2a\_j a\_k, a\_j = ∫ cos(jt)e^{t/2}dt in closed form; the archimedean integral by Gauss–Legendre quadrature with 192 nodes on \[0, 4L\] plus the exact tail −log(1 − e^{−4L})h(0), stable to 10⁻⁸⁰ under doubling the nodes; the prime terms exactly. Working precision is 2.6N + 80 decimal digits because the cosine basis on an interval of length 2L < 2π is severely ill-conditioned (the Galerkin eigenvalues in coefficient coordinates reach 10⁻²·²ᴺ); this conditioning affects nothing basis-independent, such as the sign of the form on the subspace, the L²-normalised margin computed with the Gram matrix, or the zeros of the minimiser's transform.

**Validation against the zero side.** Let M\_{jk} = 2Σ\_{γ≤γ₅₀₀} ĝ\_j(γ)ĝ\_k(γ) be the sum over the first 500 zeros. Under RH, W − M is the sum over the remaining zeros, hence positive semidefinite. At λ = 3, N = 10 we find W − M with eigenvalues between 1.6 × 10⁻⁵⁷ and 0.027, and W\_{00} = 0.0897680466576194 against M\_{00} = 0.0848879826779344. This is the only place the zeros enter, and only as a check.

**Known facts used.** RH is verified for |γ| ≤ 3 × 10¹² \[PT21\]. Weil positivity for support in \[2^{−1/2}, 2^{1/2}\] is due to Yoshida \[Yo92\] and is Theorem 1 of \[CC21\]; for support in \[e^{−0.8}, e^{0.8}\] it is Theorem 1.2 of Zhu \[Zh26\]. Beurling's sampling theorem is used in the form of \[Sei04, Thm. 2.7\]. The prolate spheroidal functions and their double orthogonality are from \[SP61\].

## 3. The detection principle

**The identity.** Let f be real and even with support in \[−L, L\] and Φ = Φ\_f. For a real r, Φ(r) is real. For a zero ρ = ½ + δ + iγ with δ ≠ 0 the functional equation and conjugation give the four zeros ρ, ρ̄, 1 − ρ, 1 − ρ̄, at ordinates (ρ − ½)/i ∈ {γ − iδ, −γ − iδ, γ + iδ, −γ + iδ}. Since ĥ(r) = Φ(r)Φ(−r) = Φ(r)² for even f, and Φ(r̄) = Φ(r)̄, their total contribution to the zero side is

```latex
2\,\Phi(\gamma - i\delta)^{2} + 2\,\Phi(\gamma + i\delta)^{2} = 4\,\mathrm{Re}\,\Phi(\gamma - i\delta)^{2} = 4\big[(\mathrm{Re}\,\Phi)^{2} - (\mathrm{Im}\,\Phi)^{2}\big],
```

with Re Φ(γ − iδ) = ∫ f cosh(δt) cos(γt) dt and Im Φ(γ − iδ) = ∫ f sinh(δt) sin(γt) dt. An on-line pair ±γ contributes 2Φ(γ)². Expanding in δ, the quadruplet contributes 4Φ(γ)² − 4δ²Φ′(γ)² + O(δ⁴), so W\_λ(f) < 0 is possible as soon as some f makes Φ vanish at the on-line zeros near γ and at γ itself while keeping Φ′(γ) large.

**The resolution height.** A function of exponential type L can interpolate prescribed values on a set of points only where the points are sparser than its Nyquist density L/π. The zeros have density (1/2π) log(γ/2π), which reaches L/π at

```latex
H(\lambda) = 2\pi e^{2L} = 2\pi \lambda^{2}.
```

Below H(λ) the zeros are sparse for PW\_L and a displaced zero is detectable for every δ > 0; above H(λ) the on-line zeros form a sampling set and detection requires the exponential mechanism (Im Φ)² ≈ (∫₀ᴸ φ sinh δt)² to beat an on-line contribution of order (log γ)∫₀ᴸφ², that is roughly e^{2δL} ≳ 2δ log γ, which gains only like (log log γ)/(2L).

**Experiments.** We take the first 2000–3000 true zeros, replace the pair nearest a chosen height by the quadruplet with parameter δ, and compute the smallest eigenvalue of the resulting zero-side matrix in the cosine basis (N = 160–250, 100–120 digits). A negative value proves indefiniteness on that subspace, hence of the form.

| λ | H(λ) | height of the displaced zero | δ | smallest eigenvalue |
| --- | --- | --- | --- | --- |
| 5 | 157 | 98.8 | 0.02 | −3.6 × 10⁻³ |
| 5 | 157 | 98.8 | 0.10 | −1.0 × 10⁻¹ |
| 5 | 157 | 98.8 | 0.30 | −1.2 |
| 4 | 100.5 | 98.8 | 0.30 | −9.1 × 10⁻³ |
| 4 | 100.5 | 111.0 | 0.30 | −3.1 × 10⁻⁷ |
| 4 | 100.5 | 124.3 | 0.30 | −1.8 × 10⁻¹² |
| 4 | 100.5 | 150.1 | 0.30 | −1.0 × 10⁻³³ |
| 4 | 100.5 | 201.3 | 0.30 | precision floor (10⁻⁹⁹) |
| 3.5 | 77 | 59.3 | 0.30 | −6.1 × 10⁻² |
| 3.5 | 77 | 98.8 | 0.30 | −2.4 × 10⁻¹² |
| 3 | 56.5 | 98.8 | ≤ 0.49 | below 10⁻³⁰ |

The λ = 5 values scale as δ² and change in the fourth digit when the number of zeros is varied from 800 to 3000 or the basis from 120 to 200 modes. The switch-off above H(λ) is steep but not a wall: at λ = 4 the indefiniteness loses about five orders of magnitude per ten units of height beyond H.

**Consequence.** RH is verified to height 3 × 10¹² \[PT21\], so a finite scale can carry new information about displaced zeros only if H(λ) exceeds that, λ ≳ 7 × 10⁵, with primes up to 5 × 10¹¹ entering. Conversely, positivity at a small scale is a statement about zeros at all heights, but only at resolution 1/L: it can fail only through a window of width about 1/L at some height in which the displaced zeros outweigh the on-line ones. Zero-density estimates control such configurations on average over long ranges, not in individual windows; that is the precise gap between what is known and positivity at any λ > √2.

## 4. Reconstruction of the zeros: the range law

**Procedure.** Let f\_λ be the eigenvector of the smallest eigenvalue of W\_λ on span{cos(kt) : k ≤ N} and Ψ\_λ = Φ\_{f\_λ}. We locate every real zero of Ψ\_λ on (0, N\] by sign changes on a grid of step 0.05 followed by root refinement at full precision, and match each to the nearest ordinate γ of a zero of ζ (reference values to 120 digits from Arb). A match is spurious if the distance exceeds 0.5.

**Results.**

| λ | primes used | N | H(λ) | zeros matched | spurious | digits at γ₁ = 14.13 | digits at the last matched zero |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 2, 3 | 60 | 25.1 | 9 of 13 below 60 | 3 | 10.5 | 0.7 at γ = 59.3 |
| 3 | up to 9 | 160 | 56.5 | 45 of 58 below 160 | 8 | 33.6 | 0.8 at γ = 158.9 |
| 5 | up to 25 | 120 | 157 | 38 of 38 below 120 | 0 | 104.7 | 43.8 at γ = 118.8 |

At λ = 3 the number of correct digits falls linearly with height: 33.6 at γ₁, 15.8 at γ₁₂ = 56.4 ≈ H, 10.4 at γ₁₇ = 69.5, 4.2 at γ₂₅ = 88.8, 2.4 at γ₂₉ = 98.8, and matches become unreliable above 110. A least-squares line gives digits ≈ 30.5 − 0.25γ, crossing zero at 123 = 2.2 H(3); at λ = 2 the line is 11.9 − 0.23γ, crossing at 51.5 = 2.05 H(2). At λ = 5 the basis limit N = 120 lies below H, so the fall-off near 120 is partly the edge of the basis; the accuracy at the first zero, 104.7 digits, is 0.67 H(5), against 0.59 H(3) and 0.42 H(2).

**The law.** Below H(λ) the reconstruction is exact to a precision of order exp(−cH(λ)) with c between 1 and 1.5 in natural logarithms; above H(λ) the information decays at a fixed rate, about one decimal digit per four units of height, and is exhausted near 2H(λ). The primes in (λ², λ′²\] therefore fix the zeros in (2πλ², 2πλ′²\], range by range, and positivity at scale λ′ is the consistency condition between consecutive ranges.

**Interpretation: Gauss quadrature of the zero measure.** Under RH, W\_λ(f) = Σ\_γ Φ\_f(γ)² is the quadratic form of the discrete measure ν = Σ\_γ δ\_γ on transforms of functions supported in \[−L, L\]. The minimiser on an (N+1)-dimensional subspace is the analogue of the orthogonal polynomial of degree N+1 for ν, and its zeros are the nodes of the corresponding Gauss quadrature. For a discrete measure the nodes converge to the atoms once the degrees of freedom exceed the number of atoms in the range, which is what the range law describes: the count of zeros below H(λ), about 2λ² log λ, equals the Nyquist count of a function supported in \[−L, L\] seen at the zero density. The number of primes used, about λ²/(2 log λ), is smaller by a factor 4(log λ)²: the information is carried by the exact values of log p in the kernel, not by the count of primes.

**Interpretation: Krein extension.** Independently of RH, Krein's extension theorem says W\_λ ≥ 0 if and only if some positive even measure ν on ℝ has Fourier transform equal to the prime kernel of W\_λ on \[−2L, 2L\]. The set E\_λ of such measures is convex and determined by the primes up to λ²; RH says E\_λ ≠ ∅ for all λ, and then Σδ\_γ lies in all of them. The maximum-entropy extension has spectral density proportional to 1/|Ψ\_λ|², so its atoms are the zeros of Ψ\_λ; the reconstructed zeros above are its support, computed without reference to ζ. This is the phenomenon studied by Connes and Consani \[CC23\] and Connes, Consani and Moscovici \[CCM23\] with prolate functions; the cosine Galerkin basis reproduces it and gives the quantitative range law.

## 5. The scale-λ equation: the minimiser is Riemann's Φ windowed

**Riemann's representation.** With Ξ(r) = ξ(½ + ir),

```latex
\Xi(r) = \int_{-\infty}^{\infty} \Phi(u)\,e^{iru}\,du, \qquad \Phi(u) = \sum_{n\ge1}\big(2\pi^{2}n^{4}e^{9u/2} - 3\pi n^{2}e^{5u/2}\big)e^{-\pi n^{2}e^{2u}},
```

Φ is even, positive, and decays doubly exponentially: Φ(0) = 0.447, Φ(0.73) = 5.9 × 10⁻⁴, Φ(1.10) = 1.4 × 10⁻⁹. (We verified the representation numerically: the integral reproduces ξ(½) = 0.4971 and vanishes at γ₁ to 10⁻¹⁸.)

**The window.** Define w\_λ(t) = f\_λ(t)/(cΦ(t)) with c fixed by w\_λ(0) = 1. Sampled at s = t/L:

| s = t/L | 0.30 | 0.42 | 0.54 | 0.60 | 0.66 | 0.72 | 0.78 | 0.84 | 0.90 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| λ = 2 (N = 45) | 0.919 | 0.821 | 0.668 | 0.570 | 0.461 | 0.345 | 0.232 | 0.132 | 0.056 |
| λ = 3 (N = 45) | 0.901 | 0.762 | 0.535 | 0.396 | 0.256 | 0.134 | 0.051 | 0.012 | 0.001 |
| λ = 5 (N = 120) | 0.931 | 0.806 | 0.558 | 0.392 | 0.223 | 0.090 | 0.020 | 0.0015 | 1.2 × 10⁻⁵ |

So f\_λ = cΦ w\_λ with w\_λ a smooth even taper, 1 near the origin and vanishing at the edge of the support. At λ = 3 and 5 the windows agree to within 0.02 up to s = 0.6 and then separate, the edge closing earlier relative to the support as λ grows; λ = 2 is gentler throughout. Among the one-parameter families exp(−βs²), the prolate ground function ψ₀(c, s), and the bump exp(−βs²/(1 − s²)), only the bump fits (β = 0.80, 1.76, 2.59; rms of log w on s ≤ 0.9 of 0.12, 0.13, 0.30). The scaling t/L is therefore not the right one for the edge; a universal profile, if one exists, needs λ = 4, 7 and 10 to be seen. The λ = 5 window computed with N = 60 was not converged (it exceeded 1 in the middle) and is discarded; the margin at λ = 5 is likewise not converged at N = 60 (Section 6).

**The envelope.** G\_λ = Ψ\_λ/Ξ is smooth and positive from r = 0 to the horizon at every scale computed. At λ = 3 (N = 60), log₁₀ G rises from −47.6 at r = 0 to −44.2 at r = 57 with no sign change and fits −47.6 − 0.0013r + 0.00106r² with rms 0.005; the first sign changes appear at r = 61, 62, 65, where the reconstruction ends. At λ = 5 (N = 60, resolved to r ≈ 54) the envelope is nearly flat; the curvature of G is not universal either.

**The equation.** The minimiser satisfies the Euler–Lagrange equation of the form, (K\_λ ∗ f)(t) = m(λ) f(t) on \[−L, L\], with K\_λ the archimedean kernel together with the impulses −Λ(n)n^{−1/2} at ±log n and the pole kernel. Substituting f = Φw:

```latex
\big(K_\lambda * (\Phi\,w_\lambda)\big)(t) = m(\lambda)\,\Phi(t)\,w_\lambda(t), \qquad |t| \le \log\lambda .
```

Read on the prime side this says that the primes up to λ², through K\_λ, reproduce the theta-function side Φ on the window \[−log λ, log λ\], with an error governed by the tail of Φ beyond the window, which is of order exp(−πλ²) = exp(−H/2). The accuracy measured at the first zero, exp(−(1 to 1.5)H), is of that form. Under RH the range law is a statement in sampling theory about how small a function of the form Φw, supported in \[−L, L\], can be on the zero set; unconditionally it describes the maximum-entropy pseudo-zero measure. The equation defines w\_λ without reference to any zero. It does not bear on RH: whether the transform of a windowed Φ has only real zeros says nothing about Ξ itself, the Pólya and de Bruijn multiplier theorems running in the opposite direction.

## 6. Margin, rigidity of the prime weights, and the first-order response

**The margin.** Let m(λ) = inf W\_λ(f)/‖f‖²\_{L²} over the Galerkin space, computed as the smallest generalised eigenvalue of (W, G) with G the Gram matrix of the cosines. Galerkin values are upper bounds for the true infimum.

| λ | N | log₁₀ m(λ) | converged? |
| --- | --- | --- | --- |
| 2 | 20 / 30 / 45 / 60 | −12.10 / −12.11 / −12.12 / −12.12 | yes (H = 25 < N) |
| 3 | 30 / 45 | −31.2 / −36.6 | no (H = 57 > N) |
| 4 | 30 / 60 / 90 | −40.8 / −60.9 / −71.0 | nearly (H = 101 ≈ N); steps −20 then −10 |
| 5 | 30 / 60 | −47.8 / −74.2 | no (H = 157 > N) |
| 6 | 90 / 130 | −108.6 / −132.4 | no (H = 226 > N) |

The margin is super-exponentially small in λ and the true values at λ ≥ 3 are below those shown. The Galerkin value converges only once the basis resolves frequencies up to the resolution height: cos(kt) carries frequency k, so N must exceed H(λ) = 2πλ², which is why λ = 2 (H = 25) is converged at N = 20 while λ = 6 (H = 226) is still falling at N = 130. In the Gauss-quadrature picture of Section 4 the margin is the leakage of the minimiser onto the zeros above the resolution height, Σ\_{γ > H} |Ψ\_λ(γ)|², which is governed by exactly the zeros that scale λ cannot resolve. The same quantity, as a function of L = log λ, is studied by Zhu \[Zh26\], who gives certified upper bounds λ\*(L): 1.05 × 10⁻⁶ at L = 0.5, 2.27 × 10⁻¹⁷ at L = 0.8, 8.08 × 10⁻³⁰ at L = 1.0, 9.98 × 10⁻⁴⁹ at L = 1.2, and 3.19 × 10⁻²⁸³ at L = 2.0, together with the empirical law −ln λ\*(L) ≈ 2π² N(T\*)/ln N(T\*), T\* = 2πe^{2L} = H(λ), the Landau–Widom plunge rate of the time-band limiting operator. Our values interleave consistently with that table: log₁₀ m(2) = −12.12 (L = 0.693) sits between Zhu's −5.98 at L = 0.5 and −16.64 at L = 0.8, and our upper bound −36.6 at λ = 3 (L = 1.099) sits between his −29.1 at L = 1.0 and −48.0 at L = 1.2.

**Rigidity.** Replace Λ(n) by (1 + ε)Λ(n) for one prime power n. To first order the smallest eigenvalue moves by ε T\_n(f₀), where T\_n(f₀) is the n-th prime term of W\_λ evaluated on the minimiser, so the form becomes indefinite at |ε| ≈ ε\_n := m(λ)/|T\_n(f₀)|. Since indefiniteness on a subspace implies indefiniteness of the form, ε\_n is a rigorous upper bound on the admissible relative change of Λ(n); a direct bisection confirms the first-order values.

| λ (N) | n = 2 | n = 3 | n = 5 | n = 7 | largest n present |
| --- | --- | --- | --- | --- | --- |
| 2 (45) | 10⁻¹⁰·⁶ | 10⁻⁷·⁵ | — | — | n = 3 |
| 3 (45) | 10⁻³⁵·³ | 10⁻³³ | 10⁻²⁷ | 10⁻²⁰ | n = 8: 10⁻¹²·³ |
| 5 (60) | 10⁻⁷³·¹ | 10⁻⁷² | 10⁻⁶⁹ | 10⁻⁶⁶ | n = 23: 10⁻¹⁷·⁵ |

(Values for n = 3, 5, 7 at λ = 3 and 5 are read from the N = 30 run, shifted by the change of m(λ) between N = 30 and the larger N; the n = 2 and last-prime values are from the larger run.) The prime most recently entered is the least constrained, because the minimiser's autocorrelation at a lag near 2L is small; the prime 2 is the most constrained. Positivity at λ = 5 determines Λ(2) to 73 digits. The consequence for strategy is quantitative: the next prime enters with weight Λ(n)/√n ≈ 1/λ against a slack of 10⁻⁷⁴, so positivity at the next scale is not positivity at this scale plus a perturbation, and no inductive argument over primes can proceed through the margin.

**First-order response of the reconstructed zeros.** Let S\_{jn} = ∂γ\_j(λ)/∂ε\_n be the displacement of the j-th reconstructed zero per unit relative change of Λ(n), computed by finite differences at ε = 10⁻⁴⁵ and 10⁻⁴⁴ (agreement to 10⁻¹⁰) at λ = 3, N = 45, for the eight zeros below 45. The expected signature of the explicit formula, a factor cos(γ\_j log n), is absent at this order. Instead S is close to rank one: every zero moves in the same direction, with an amplitude A(γ\_j) rising by about two orders of magnitude per zero (4 × 10² at γ₁ to 3 × 10¹⁶ at γ₈ for n = 2) and a prime factor B(n) falling steeply with the lag (1, 10⁻³, 10⁻⁶·⁵, 10⁻¹⁵, 10⁻²² for n = 2, 4, 5, 7, 8 relative to n = 2). First-order perturbation theory explains this: the correction to the minimiser is dominated by the second eigenfunction e₁, whose eigenvalue is nearest, so the response direction is e₁ for every n, and A(γ\_j) ≈ −Ψ₁(γ\_j)/Ψ₀′(γ\_j), largest where the zero is least accurately pinned. Whether this two-mode structure survives in the L²-normalised problem as N → ∞, where the eigenvalue spread is not inflated by the conditioning of the cosine basis, has not been checked.

## 7. Certified positivity

**Theorem 1.** For λ = 2, W\_2 is positive definite on V\_N = span{cos(kt) : 0 ≤ k ≤ N} ⊂ L²\[−log 2, log 2\] for N = 20 and N = 30.

**Proof (computer-assisted, ball arithmetic).** All quantities are computed in Arb through python-flint, so that each is a ball guaranteed to contain the true value. (a) The autocorrelations h\_{jk}(x) for 0 ≤ x ≤ 2L, the pole coefficients a\_j and the prime terms Λ(n)n^{−1/2}h\_{jk}(log n) are closed trigonometric and exponential expressions, evaluated in ball arithmetic. (b) The archimedean integral over \[ε, 4L\], ε = 10⁻¹²⁰, is computed with Arb's rigorous integration (acb\_calc\_integrate) to tolerance 10⁻¹⁸⁰; the integrand is analytic on the path. The contribution of \[0, ε\], where the integrand has a removable singularity, is bounded by 2Mε with M a bound on the derivative of the numerator, M ≤ (5/4)|h(0)| + (2L(j + k) + 1)/2, and added as a ball. The tail beyond 4L is the exact −h(0) log(1 − e^{−4L}). (c) The resulting matrix W has entries of radius below 10⁻¹¹⁸. (d) An approximate orthogonal eigenbasis Q is computed in 200-digit floating point; the ball matrix B = QᵀWQ is formed in ball arithmetic, and Gershgorin's theorem is applied with balls: for every row, the lower bound of the diagonal entry exceeds the sum of the upper bounds of the moduli of the off-diagonal entries. The determinant of Q is enclosed in a ball not containing 0. Hence every matrix in the ball B is positive definite, so QᵀW\_{true}Q and therefore W\_{true} are positive definite. For N = 30 the worst Gershgorin gap is 7.7 × 10⁻⁹⁰ ± 10⁻⁹⁶ (coefficient coordinates; the L² margin is about 10⁻¹²); for N = 20 it is 9.3 × 10⁻⁶⁰. The two inputs not themselves machine-checked are the closed forms in (a), which were validated against the zero side to 10⁻⁵⁷ (Section 2), and the derivative bound in (b), which is elementary. ∎

**Theorem 1′ (λ = 3).** The same certificate holds at λ = 3, where the prime powers 2, 3, 4, 5, 7, 8 enter: W\_3 is positive definite on span{cos(kt) : 0 ≤ k ≤ 45} ⊂ L²\[−log 3, log 3\]. Entry radii are below 2.1 × 10⁻¹¹⁸, the worst Gershgorin gap is 3.8 × 10⁻⁹⁸ ± 3 × 10⁻¹⁰⁴ in coefficient coordinates (L² margin about 10⁻³⁷), and det Q = 1.0000 ± 3 × 10⁻¹⁰. Construction time 28 minutes at 200 digits.

**The tail: why a prolate splitting of W fails, and what works instead.** To pass from V\_N to all of L²\[−L, L\] one needs W\_2 ≥ 0 on the orthogonal complement together with control of the cross terms, against a margin of 10⁻¹². Separating frequencies with prolate spheroidal functions ψ\_n of bandwidth Ω on \[−L, L\] \[SP61\] does not do it by itself: we computed the prime shift operator P in the prolate basis at Ω = 80 (c = ΩL = 55.5, 90 even modes) and found ‖P\_LH‖ ≈ 0.8 for every cut K between 10 and 80, because the modes in the Landau–Widom transition region (n ≈ 30–50, concentration eigenvalues between 10⁻² and 0.99) are neither in-band nor out-of-band and a translation by log 2 or log 3 couples them with norm of order one; only the fully in-band modes n ≤ 12 decouple, with row norms below 10⁻⁹. The remedy is to bound the whole symbol before splitting. In frequency,

```latex
W_\lambda(f) = \frac{1}{2\pi}\int_{\mathbb R} |\hat f(r)|^2\,\sigma_\lambda(r)\,dr + 2\hat f(i/2)^2,\qquad \sigma_\lambda(r) = \operatorname{Re}\psi\!\left(\tfrac14+\tfrac{ir}{2}\right) - \log\pi - 2\sum_{\log n<2L}\frac{\Lambda(n)}{\sqrt n}\cos(r\log n),
```

and since Re ψ(¼ + iy) increases with |y|, σ\_λ(r) ≥ m\_out := Re ψ(¼ + iΩ/2) − log π − 2S for |r| ≥ Ω, where S = ΣΛ(n)/√n over the prime powers present. For λ = 2 (S = 1.1244) this is positive as soon as Ω > 2πe^{2S} = 59.6. Replacing σ by m\_out outside \[−Ω, Ω\] gives a form B ≤ W\_λ in which every term except m\_out‖f‖² factors through the band-limiting operator G\_in = Π M\_{1\[−Ω,Ω\]} Π of Slepian, so that in any orthonormal basis the entries of B − m\_out I are bounded by (max|σ| + m\_out)√(ℓ\_j ℓ\_k), ℓ\_k the in-band energy of the k-th basis function. This is the reduction of Zhu \[Zh26, Thm. 1.1\], which we found independently; it makes the cross terms small because the primes are inside the symbol, not outside it.

**Theorem 2 (positivity of W\_2 on the whole window).** For every real even f ∈ L²\[−log 2, log 2\],

```latex
W_2(f) \;\ge\; 1.38\times 10^{-13}\,\|f\|_{L^2}^2 .
```

Together with the Galerkin value m(2) = 10^{−12.12} this encloses the exact margin: 1.38 × 10⁻¹³ ≤ inf W\_2(f)/‖f‖² ≤ 7.6 × 10⁻¹³.

**Proof (computer-assisted, ball arithmetic).** Take Ω = 80, so m\_out = 0.29532 and max\_{|r|≤Ω}|σ\_2| ≤ 7.622 (the archimedean part is monotone, so its modulus is at most max(|ψ(¼)|, Re ψ(¼ + 40i))). Basis: the orthonormal even Legendre polynomials p\_k(t) = √((2k+1)/2L) P\_k(t/L), k = 0, 2, …, D with D = 110 (K = 56 modes), whose transforms are p̂\_k(r) = √((2k+1)L/2)·2 i^k j\_k(rL) with j\_k the spherical Bessel function, and whose pole coefficients ⟨p\_k, cosh(t/2)⟩ = √((2k+1)L/2)·2 i\_k(L/2) are modified spherical Bessel values. (a) The block B\_PP = \[(1/π)∫₀^Ω p̂\_j p̂\_k (σ − m\_out) dr\] + m\_out I + 2χχᵀ is computed with Arb's rigorous integration (the integrand is analytic on the path; Re ψ is continued as (ψ(¼ + ir/2) + ψ(¼ − ir/2))/2), with entry radii below 1.6 × 10⁻¹⁷; Gershgorin's theorem on QᵀB\_PPQ with Q an approximate eigenbasis gives λ\_min(B\_PP) ≥ 1.3805 × 10⁻¹³ (floating-point value 1.3892 × 10⁻¹³). (b) Tail: by |J\_ν(x)| ≤ (x/2)^ν/Γ(ν + 1) the in-band energy of p\_k is ℓ\_k ≤ (2/π)(ΩL)^{2k+1}/((2k+1)!!)², so for f orthogonal to the block ⟨f, G\_in f⟩ ≤ ε\_Q‖f‖² with ε\_Q ≤ Σ\_{k>D}ℓ\_k ≤ 1.08 × 10⁻⁴²; the in-band energies of the block functions themselves fall from 0.994 (k = 0) through 0.5 (k = 38) and 0.02 (k = 56) to 10⁻²³ (k = 88). (c) Hence B\_QQ ≥ m\_out − (7.622 + m\_out)ε\_Q ≥ 0.2953, ‖B\_PQ‖ ≤ (7.622 + m\_out)√ε\_Q + 2‖χ‖‖Qχ‖ ≤ 1.2 × 10⁻¹⁹ (the Legendre tail of cosh(t/2) beyond degree 110 has norm below 10⁻¹⁹), and the Schur complement condition λ\_min(B\_PP) > ‖B\_PQ‖²/λ\_min(B\_QQ) = 4.8 × 10⁻³⁸ holds with room to spare: B ≥ (λ\_min(B\_PP) − 4.8 × 10⁻³⁸)I on L²\[−L, L\], and W\_2 ≥ B. Running time 9 minutes at 40 digits (1596 rigorous integrals). ∎

**Remarks.** (1) The same statement for L ≤ 0.8 (autocorrelation support 1.6, still only the primes 2 and 3 with 4) is Zhu's Theorem 1.2 \[Zh26\], with the lower bound 8.9 × 10⁻¹⁸, proved by the same reduction with a Legendre block of 200 modes and a verified Cholesky factorisation; Theorem 2 is an independent certification at L = log 2 (support 1.386), with a Gershgorin certificate and the explicit two-sided enclosure of the margin. (2) The prime 5 enters at 2L = log 5 = 1.609, just above Zhu's 0.8, and raises 2πe^{2S} from 119 to 503; our λ = 3 (2L = 2.197, six prime powers, 2πe^{2S} = 3570, margin 10⁻³⁷) would need a block of several thousand Legendre modes. This is Zhu's barrier \[Zh26, Thm. 1.4\]: within pointwise envelopes the cutoff 2πe^{2S} is optimal because the prime comb attains 2S, so the block size is doubly exponential in L and the method stops near L ≈ 3. (3) The compressed prime mass S\_eff of the next paragraph cannot simply replace S in this argument: the envelope must bound the comb pointwise in frequency, and replacing the in-band prime terms by an operator-norm bound destroys the cancellation on the low block that the 10⁻¹² margin lives on.

**Theorem 3 (certified rigidity of the prime data at λ = 2).** Let W\_2^{(ε)} be the form W\_2 with Λ(2) replaced by (1 + ε)Λ(2), and W\_2^{\[η\]} the form with log 2 replaced by log 2 + η in the prime term (weight unchanged); define the analogous forms for the prime 3. Each of the following forms is indefinite on L²\[−log 2, log 2\]:

| perturbed datum | indefinite for | certified value W(f₀) |
| --- | --- | --- |
| Λ(2) → (1 + ε)Λ(2) | every ε ≥ 2.62 × 10⁻¹¹ | −6.03 × 10⁻²⁸ ± 10⁻³⁴ |
| Λ(2) → (1 + ε)Λ(2) | every ε ≤ −1.84 × 10⁻⁷ | −1.81 × 10⁻¹⁴ ± 10⁻¹⁹ |
| Λ(3) → (1 + ε)Λ(3) | every ε ≥ 3.61 × 10⁻⁸ | −8.06 × 10⁻²⁸ ± 10⁻³³ |
| Λ(3) → (1 + ε)Λ(3) | every ε ≤ −1.61 × 10⁻⁵ | −5.95 × 10⁻²⁴ ± 10⁻²⁹ |
| log 2 → log 2 + η | η = +2.55 × 10⁻⁸ | −2.46 × 10⁻¹⁴ ± 10⁻¹⁹ |
| log 2 → log 2 + η | η = −2.33 × 10⁻¹² | −6.17 × 10⁻²⁸ ± 10⁻³³ |
| log 3 → log 3 + η | η = +6.59 × 10⁻⁷ | −2.02 × 10⁻²³ ± 10⁻²⁹ |
| log 3 → log 3 + η | η = −1.23 × 10⁻⁹ | −7.85 × 10⁻²⁸ ± 10⁻³³ |

With Theorem 2, the set of ε for which W\_2^{(ε)} is positive is therefore an interval contained in (−1.84 × 10⁻⁷, 2.62 × 10⁻¹¹) for the prime 2 and in (−1.61 × 10⁻⁵, 3.61 × 10⁻⁸) for the prime 3: positivity on the window \[1/2, 2\] alone pins Λ(2) to eleven digits from above and log 2 to twelve digits from below.

**Proof.** For each row an explicit even trigonometric polynomial f₀ = Σ\_{k≤24} a\_k cos(kt) is exhibited (the minimiser of the floating-point Galerkin matrix of the perturbed form at the stated parameter, which lies 5 % beyond the bisected boundary), and the value of the perturbed form on f₀ is computed in ball arithmetic from the rigorous matrix of Theorem 1 (N = 24, 142 digits, entry radii below 7 × 10⁻¹¹⁹) and the closed-form prime matrices at the perturbed weight or position; the ball lies strictly below zero. Since f₀ ∈ L²\[−L, L\], the perturbed form is indefinite on the whole space. For the weight perturbations the dependence on ε is affine with T\_n(f₀) < 0 (resp. > 0), so the conclusion extends to every ε beyond the certified one. ∎

**Remark.** The asymmetry is itself information: lowering log 2 by 2 × 10⁻¹² destroys positivity while raising it by 2 × 10⁻⁸ does not, and lowering Λ(2) is tolerated 10⁴ times more than raising it. The first-order prediction ε\* = m/|T\_n(f₀)| of Section 6 (10^{−10.6} for the prime 2) matches the certified upper boundary 2.5 × 10⁻¹¹. Zhu's upper bounds on λ\*(L) make the analogous statements at larger L sharper still (the margin at L = 2 is 3 × 10⁻²⁸³), but they remain statements about the truncated form: it is the window, not the primes, that is rigid.

**Certification cost and the compressed prime mass.** On frequencies above 2πe^{2S} positivity is automatic, so the block size scales like 4L e^{2S}. The naive S = Σ\_{n≤λ²}Λ(n)/√n (with log n < 2L) can be replaced by S\_eff(λ), the norm of the operator Σ\_n (Λ(n)/√n)(S\_{log n} + S\_{−log n})/2 compressed to L²\[−L, L\], which we computed in the cosine Galerkin basis with the Gram matrix (a lower bound on the true norm, converging from below):

| λ | S | S\_eff | S\_eff/S | block size with S | block size with S\_eff |
| --- | --- | --- | --- | --- | --- |
| 2 | 1.12 | 0.45 | 0.40 | 26 | 7 |
| 3 | 3.17 | 1.14 | 0.36 | 2.5 × 10³ | 43 |
| 4 | 4.97 | 1.82 | 0.37 | 1.2 × 10⁵ | 210 |
| 5 | 7.16 | 2.48 | 0.35 | 1.1 × 10⁷ | 910 |
| 10 | 16.9 | 5.54 | 0.33 | 4 × 10¹⁵ | 6 × 10⁵ |

S\_eff measures how much of the prime mass the window actually feels, and it is the right constant for operator-norm arguments (the per-prime criterion of Section 8 and the frame identity of Section 9); it is not the right constant for the envelope certificate of Theorem 2, which must bound the prime comb pointwise in frequency, and there Zhu's barrier \[Zh26, Thm. 1.4\] shows 2S is optimal. So the block size is exponential in λ/log λ and the method stops near λ ≈ e³; the reason is structural: the archimedean multiplier grows only logarithmically in frequency while the prime mass grows linearly in λ, and short Dirichlet polynomials ΣΛ(n)n^{−1/2−ir} attain their maximum at some frequencies below e^{cλ}. Any certification that does not use the cancellation in that polynomial frequency by frequency pays e^{2S\_eff}, and that cancellation is itself a statement about zeros.

## 8. Discussion

**What is new here.** The explicit detection identity 4 Re Φ(γ − iδ)² and the resolution height H(λ) = 2πλ² with its numerical confirmation; the quantitative range law for reconstruction from finitely many primes, including the 105-digit value of γ₁ from the primes up to 25; the identification f\_λ = cΦw\_λ of the minimiser with Riemann's Φ windowed, and the equation it satisfies; the measured margins and the rigidity of the prime weights; the rank-one first-order response; a ball-arithmetic certificate of W\_2 ≥ 1.38 × 10⁻¹³ on the whole window with a two-sided enclosure of the margin (Theorem 2, independent of and consistent with Zhu \[Zh26\]); the first certified rigidity statement for the prime data (Theorem 3); the measurement of the compressed prime mass; the coercivity and localisation theorems under RH (Theorem 4); and the quantitative finding that a prolate splitting of the form itself fails through the transition modes, so that the primes must be placed inside the symbol.

**What is known and closely related.** The reconstruction of low zeros from the Weil form on an interval is the phenomenon of Connes and Consani's ζ-cycles \[CC23\], where the first thirty-one zeros are obtained with prolate functions, and of the prolate wave operators of Connes, Consani and Moscovici \[CCM23\]; our contribution there is the cosine basis, the explicit prime-side construction with validation, and the measured range law. The positivity of the archimedean form for support in \[2^{−1/2}, 2^{1/2}\] is \[Yo92, CC21\], and its extension to \[e^{−0.8}, e^{0.8}\] with primes is \[Zh26\]; the Galerkin truncations themselves are those of Connes–van Suijlekom and \[CCM23\], studied as such by Groskin \[Gr26\]. That the Weil minimiser's zeros are Gauss nodes is a direct reading of the moment-problem structure and may be folklore among specialists; we have not found it written. Krein's extension theorem and the maximum-entropy extension are classical. The identification with Φ is, to our knowledge, not in the literature in this form, though it is the kind of statement experts would expect once stated.

**What is not claimed.** No zero-free region, no statement about zeros above the heights treated, no reduction of RH. Theorems 1 and 1′ are finite-dimensional; Theorem 2 is a statement about one window, \[1/2, 2\], inside the range already covered by Zhu \[Zh26\]; Theorem 3 is a statement about the truncated form, not about the primes. The range law, the window structure and the margin values are numerical observations at three scales; their forms (slope 0.25 digits per unit height, horizon 2H, accuracy exp(−cH)) are empirical fits. The rank-one response may be partly an artefact of the cosine basis conditioning. The parts of the detection analysis above H(λ) are heuristic.

**Relation to the de Bruijn–Newman flow.** A window on Φ in the variable u is close in spirit to the de Bruijn family H\_t = ∫ e^{tu²}Φ(u)e^{iru}du at negative t, for which Rodgers and Tao \[RT20\] show that non-real zeros must appear (Λ ≥ 0). The scale-λ minimiser keeps its zeros real up to the horizon and we have no information beyond; whether the tapered Φ of Section 5 has non-real zeros above the horizon, and at what height, is a natural question connecting the two pictures.

**Relation to the 2026 literature.** Three preprints, found after the computations of this paper were complete, work on the same object. Zhu \[Zh26\] studies λ\*(L) = inf W(f)/‖f‖² on \[−L, L\] from both sides: a finite reduction through a pointwise envelope of the Weil symbol (the argument of our Theorem 2, found independently), certified positivity for L ≤ 0.8 with the lower bound 8.9 × 10⁻¹⁸, certified upper bounds up to L = 2, the same resolution height T\* = 2πe^{2L} = H(λ) by the same density argument, the Landau–Widom law for the decay of λ\*(L), and a barrier theorem: the envelope method needs a frequency cutoff above 2πe^{A\_L} with A\_L = 2ΣΛ(n)/√n, doubly exponential in L, so it cannot be pushed past L ≈ 3.2. Groskin \[Gr26\] gives an exact dictionary between Galerkin coefficient vectors of the Connes–van Suijlekom and Connes–Consani–Moscovici truncations and band-limited Guinand–Weil test functions whose zero sums equal the quadratic values, and a two-sided certification rule for the archimedean cutoff. Suzuki \[Su26\] places the results of Yoshida, Bombieri and Connes–Consani in one operator framework through the screw function and de Branges spaces. Against this background the contributions specific to the present paper are: the reconstruction of the zeros from the prime-side minimiser and its range law (Section 4), the detection identity for off-line zeros (Section 3), the factorisation f\_λ = Φ w\_λ of the minimiser and the scale equation (Section 5), the rigidity of the prime data with its certified version, Theorem 3 (Sections 6 and 7), and the frame identity (Section 9). Our block certificates (Theorems 1 and 1′) and the margin values are consistent with, and for L = log 2 subsumed by, Zhu's results; Theorem 2 is an independent certification of his Theorem 1.2 at L = log 2 with a different implementation.

**The honest frontier.** Every criterion we formulated, detection, rigidity, a per-prime reflection-coefficient inequality in the Krein system of the prime kernel, reduces to the one quantity m(λ), the positivity margin at scale λ, which is the leakage of Φw\_λ onto the unresolved zeros. RH is the statement m(λ) > 0 for all λ. Seen from the zeros, m(λ) is a sampling constant of the zero set for functions of exponential type L, which fails only through a cluster of displaced zeros in a window of width 1/L at some height, a configuration that zero-density estimates control on average but not in windows. Seen from the primes, it is the regularity of the Krein system driven by the von Mangoldt impulses. Connes and Consani supply the positive structure for the archimedean part alone; Zhu's envelope, and our Theorem 2, buy the first two primes by brute force on one window and prove nothing beyond it; the same structure with the shifts included at every scale is the open problem, and nothing in this paper supplies it.

## 9. The Riemann-wave frame: one identity and a dictionary

The structure f\_λ = Φw of Section 5 turns the Weil form into the frame operator of an explicit system of functions, and that places the whole scale-λ theory inside the theory of exponential systems on an interval, where the relevant theorems exist.

**The identity.** For w ∈ L²\[−L, L\] put f = Φw. Since the transform of Φ is Ξ and Ξ vanishes at the ordinates γ\_k, the zero side of the explicit formula for f reads, when all zeros are on the line,

```latex
\sum_{k} \big|\langle w,\ \Phi\, e^{i\gamma_k t}\rangle_{L^2[-L,L]}\big|^{2} \;=\; \langle w,\ \Phi\, K_\lambda\, \Phi\, w\rangle ,
```

where K\_λ is the prime kernel of W\_λ (archimedean kernel, pole kernel, impulses −Λ(n)n^{−1/2} at ±log n for n ≤ λ²) and the right side is explicit in the primes. The left side is the frame operator of the Riemann waves e\_k(t) = Φ(t)e^{iγ\_k t} on the interval. If a zero is off the line, ρ′ = γ − iδ, its term becomes ⟨w, Φe^{δt}e^{iγt}⟩⟨w, Φe^{−δt}e^{iγt}⟩̄, no longer a square, which is Proposition 2 in this language.

**The dictionary.** With the identity, each object of this paper is a classical object of the theory of exponential systems in Paley–Wiener spaces:

| This paper | Exponential systems on \[−L, L\] | Theorem available |
| --- | --- | --- |
| resolution height H(λ) = 2πλ² | the height where the zero density reaches the Nyquist density L/π | Landau's necessary density conditions for sampling and interpolation |
| exact reconstruction below H | the zeros below H are an interpolating sequence for PW\_L (up to separation) | Beurling, Seip: interpolation when the upper density is below L/π |
| no detection above H, reconstruction fading to 2H | above H the zeros are a sampling set; the minimiser can only be small there, not zero | Beurling–Malliavin, Landau: sampling when the lower density exceeds L/π |
| margin m(λ) | the lower frame bound of the Riemann waves, weighted by Φ | sampling constants of Paley–Wiener spaces |
| window w\_λ | the vector of the frame operator with the smallest eigenvalue; orthogonal to the Riemann waves below H | extremal problems of Landau–Pollak–Slepian type |
| detection of δ > 0 below H | a complex frequency in a system that is not yet a frame | — |

**What the tool buys and what it does not.** Under RH, every empirical law of this paper becomes a statement about the zero set as a sequence of frequencies, and the theorems of Landau, Beurling, Seip and Slepian are the natural route to proving them: that the minimiser vanishes at the zeros below (1 − ε)H(λ), that its accuracy decays at a rate fixed by the excess density above H, and that log m(λ) is asymptotic to an explicit functional of the zero counting function. The first of these is Theorem 4 below; the other two are theorems to be written, not conjectures to be hoped for, and they would turn Observations 3 to 5 into results. Unconditionally the identity is exactly the explicit formula on the test functions Φw, and the positivity of the frame operator is RH; the tool organises the programme and makes its finite-scale laws provable, and it does not supply the positive structure that the shifts obstruct.

**Theorem 4 (under RH: coercivity and localisation).** Assume RH. (i) For every λ > 1 there is c(λ) > 0 such that W\_λ(f) ≥ c(λ)‖f‖²\_{L²} for every real even f ∈ L²\[−L, L\]; in particular m(λ) > 0. (ii) If ‖f‖ = 1 and W\_λ(f) = m, then |f̂(γ)| ≤ √m at every zero ½ + iγ. If moreover |f̂′(γ)| ≥ s > 0 and |f̂″| ≤ M on \[γ − δ, γ + δ\] with δ = 2√m/s and 2Mδ ≤ s, then f̂ has a real zero in (γ − δ, γ + δ).

**Proof.** (i) Under RH, W\_λ(f) = Σ\_γ f̂(γ)² with f̂ in the Paley–Wiener space PW\_L (entire of exponential type L, ‖f̂‖²\_{L²(ℝ)} = 2π‖f‖²). The zeros have density (1/2π) log(t/2π) per unit height, which exceeds the Nyquist density L/π of PW\_L above H(λ) = 2πλ². Under RH the number of zeros in \[t, t + h\] is at most C(h log t + log t/log log t) (Littlewood), so from the zeros above any height one can extract a δ-separated subsequence Γ whose lower uniform density D⁻(Γ) exceeds L/π once δ < 1/(4CL). By Beurling's sampling theorem for uniformly discrete sets with D⁻(Γ) > L/π \[Beu89; Sei04, Thm. 2.7\], Γ is a set of sampling for PW\_L: ‖F‖² ≤ C\_Γ Σ\_{γ∈Γ} |F(γ)|² for all F ∈ PW\_L. Hence W\_λ(f) ≥ Σ\_{γ∈Γ} f̂(γ)² ≥ 2π‖f‖²/C\_Γ. (ii) The first claim is Σ\_γ f̂(γ)² = m. For the second, f̂ is real on ℝ because f is real and even, and Taylor's formula gives f̂(γ ± δ) = f̂(γ) ± δ f̂′(γ) + R with |R| ≤ Mδ²/2 ≤ sδ/4 = √m/2. With δ|f̂′(γ)| = 2√m and |f̂(γ)| ≤ √m, the values f̂(γ + δ) and f̂(γ − δ) have opposite signs, so f̂ vanishes between them. ∎

**What Theorem 4 does and does not give.** Part (i) is the qualitative half of the range law and of Observation 5: the form is coercive on the window, and the margin is positive for every λ, because the zeros above H(λ) oversample PW\_L. Part (ii) turns a small margin into a localisation of the minimiser's zeros: at λ = 2, m ≈ 7.6 × 10⁻¹³ gives δ ≈ 1.7 × 10⁻⁶/s at every zero whose slope s is not small, which is the worst case over all zeros; the 10 digits observed at γ₁ reflect that f̂ is far smaller than √m there. The exponential rate above H (a quarter digit per unit height, horizon 2H) is not covered: it is a quantitative inequality of Remez or Logvinenko–Sereda type on PW\_L for a set of density below Nyquist, and stating it precisely is the open problem the law leaves. Unconditionally, (i) is exactly what Zhu's certified lower bounds \[Zh26\] establish for L ≤ 0.8, and Theorem 2 for L = log 2.

## Reproducibility and data

All computations are in the companion archive v33\_code.zip (Python 3; mpmath 1.4, numpy 2.5, scipy 1.18, python-flint 0.9 with Arb). Reference zeros: zeros500\_120d.txt (first 500 ordinates to 120 digits, from Arb's acb\_zeta\_zeros), zeros3000.txt (first 3000 to 35 digits).

| Script | Produces | Section |
| --- | --- | --- |
| weil\_prime.py | the prime-side matrix W\_λ in the cosine basis | 2 |
| checks.py | explicit-formula check on a Gaussian; W − M test | 2 |
| detect.py | zero-side matrix with one displaced zero; smallest eigenvalue | 3 |
| recon.py | zeros of the minimiser transform vs true zeros | 4 |
| phi\_test.py, window.py | Riemann's Φ check; the window w\_λ(t/L) and its fits | 5 |
| factor.py | the envelope G\_λ = Ψ\_λ/Ξ | 5 |
| margin.py, rigidity.py | L² margin m(λ); prime-weight tolerances (first order and bisection) | 6 |
| jacobian.py | first-order response of the reconstructed zeros | 6 |
| certify.py | ball-arithmetic construction of W\_2 and the Gershgorin certificate | 7 |
| s\_eff.py | compressed prime mass S\_eff(λ) | 7 |
| prolate.py, prolate\_cross.py | prolate eigenvalues λ\_n(c) at Ω = 80, 100, 120; the prime shift operator in the prolate basis and its low/high coupling | 7 |
| symbol\_B.py, certify\_symbol.py | the envelope form B in the Legendre basis (exploratory, then Arb certificate of Theorem 2) | 7 |
| certify\_rigidity.py | ball-arithmetic indefiniteness certificates of Theorem 3 (rigidity\_cert\_l2\_n24.json) | 7 |

Running times on two cores: W\_λ builds in seconds to minutes for N ≤ 60 and about 50 minutes for N = 160 at λ = 3; the N = 30 block certificate at λ = 2 takes 15 minutes, the Theorem 2 certificate 9 minutes (1596 rigorous integrals at 40 digits) and the Theorem 3 certificates 3.5 minutes. Working precision is 2.6N + 80 digits throughout, except the certificate (200 digits) and the detection runs (100–120 digits). Every number quoted in the text is in the results/ directory of the archive.

## References

- &#91;CC21\] A. Connes, C. Consani, Weil positivity and trace formula, the archimedean place, Selecta Math. (N.S.) 27 (2021), no. 4, 77. [arXiv:2006.13771](https://arxiv.org/abs/2006.13771v1).
- &#91;CC23\] A. Connes, C. Consani, Spectral triples and ζ-cycles, L'Enseignement Math. 69 (2023), 93–148. [arXiv:2106.01715](https://arxiv.org/abs/2106.01715v1).
- &#91;CCM23\] A. Connes, C. Consani, H. Moscovici, Zeta zeros and prolate wave operators. [arXiv:2310.18423](https://arxiv.org/abs/2310.18423v1).
- &#91;PT21\] D. J. Platt, T. S. Trudgian, The Riemann hypothesis is true up to 3·10¹², Bull. Lond. Math. Soc. 53 (2021), 792–797.
- &#91;RT20\] B. Rodgers, T. Tao, The De Bruijn–Newman constant is non-negative, Forum Math. Pi 8 (2020), e6. [arXiv:1801.05914](https://arxiv.org/abs/1801.05914v2).
- &#91;SP61\] D. Slepian, H. O. Pollak, Prolate spheroidal wave functions, Fourier analysis and uncertainty I, Bell Syst. Tech. J. 40 (1961), 43–63.
- &#91;Os12\] A. Osipov, Certain upper bounds on the eigenvalues associated with prolate spheroidal wave functions. [arXiv:1206.4541](https://arxiv.org/abs/1206.4541v1).
- &#91;BK14\] A. Bonami, A. Karoui, Uniform bounds of prolate spheroidal wave functions and eigenvalues decay, C. R. Math. Acad. Sci. Paris 352 (2014). [article](https://comptes-rendus.academie-sciences.fr/mathematique/articles/10.1016/j.crma.2014.01.004/).
- &#91;Weil52\] A. Weil, Sur les «formules explicites» de la théorie des nombres premiers, Comm. Sém. Math. Univ. Lund (1952), 252–265.
- &#91;Krein40\] M. G. Krein, Sur le problème du prolongement des fonctions hermitiennes positives et continues, C. R. (Doklady) Acad. Sci. URSS 26 (1940), 17–22.
- &#91;Bl26\] F. D. Blum, Configuration Space Temporality: Complete Research Archive, Zenodo, DOI 10.5281/zenodo.18859602, version 33 (corrections and withdrawals).
- &#91;OAI26\] OpenAI, openai/math, family 003, The quasi-Riemann hypothesis (preprints dated 30 September and 5 October 2026), cited for context only; unverified at the time of writing. [repository](https://github.com/openai/math).
- &#91;Zh26\] X. Zhu, Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau–Widom decay law, [arXiv:2608.24827](https://arxiv.org/abs/2608.24827) (v2, 3 September 2026).
- &#91;Gr26\] A. Groskin, A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form, [arXiv:2607.02828](https://arxiv.org/abs/2607.02828) (July 2026).
- &#91;Su26\] M. Suzuki, Weil's quadratic form via the screw function, [arXiv:2606.09096](https://arxiv.org/abs/2606.09096) (June 2026).
- &#91;Yo92\] H. Yoshida, On Hermitian forms attached to zeta functions, in: Zeta functions in geometry (Tokyo 1990), Adv. Stud. Pure Math. 21 (1992), 281–325.
- &#91;Beu89\] A. Beurling, The collected works of Arne Beurling, Vol. 2, Harmonic Analysis (L. Carleson et al., eds.), Birkhäuser, 1989, pp. 341–365 (balayage of Fourier–Stieltjes transforms; interpolation and sampling).
- &#91;Sei04\] K. Seip, Interpolation and sampling in spaces of analytic functions, Univ. Lecture Ser. 33, Amer. Math. Soc., 2004.

The last two references are cited from memory for bibliographic details not re-opened here (\[Weil52\], \[Krein40\]) and should be checked before submission.
