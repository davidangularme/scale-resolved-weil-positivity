# Scale-Resolved Weil Positivity

Frédéric David Blum (Catalyst AI Research, Haifa) — October 2026.
Companion repository of the preprint of the same title, Zenodo record [10.5281/zenodo.23229271](https://doi.org/10.5281/zenodo.23229271). The earlier *Configuration Space Temporality* archive is a separate record, [10.5281/zenodo.18859602](https://doi.org/10.5281/zenodo.18859602), whose version 33 ([10.5281/zenodo.23232290](https://doi.org/10.5281/zenodo.23232290), 8 October 2026) carries the corrections and the Resolution Law; the document in `archive-v33/` is the full audit of that archive and is included here because both papers cite it.

## Contents

| Folder | What it holds |
| --- | --- |
| `paper/` | *Scale-Resolved Weil Positivity: Reconstruction of Zeta Zeros from Finitely Many Primes, a Detection Law, and a Certified Block* (PDF + Markdown source) |
| `archive-v33/` | *CST/FDBC Archive v33 — Corrections and a Non-Circular Weil-Form Test* (PDF + Markdown): supplementary note auditing the earlier CST/FDBC archive, with the withdrawn claims and the first computations behind the paper |
| `code/` | All scripts (Python, mpmath, python-flint/Arb) and the result files quoted in the paper, including the ball-arithmetic certificates of Theorems 1–3 |

## Main results of the paper

* **Theorem 2.** For every real even f supported in [−log 2, log 2]: W₂(f) ≥ 1.38 × 10⁻¹³ ‖f‖². Together with the Galerkin value this encloses the exact margin: 1.38 × 10⁻¹³ ≤ inf W₂(f)/‖f‖² ≤ 7.6 × 10⁻¹³. (Independent certification, at L = log 2, of a theorem of X. Zhu, [arXiv:2608.24827](https://arxiv.org/abs/2608.24827), valid for L ≤ 0.8.)
* **Theorem 3.** Certified rigidity: the same form becomes indefinite when Λ(2) is raised by a relative 2.6 × 10⁻¹¹ or lowered by 1.8 × 10⁻⁷, or when log 2 is lowered by 2.3 × 10⁻¹² or raised by 2.5 × 10⁻⁸ (eight ball-arithmetic certificates).
* **Theorem 4 (under RH).** W_λ is coercive for every λ (Beurling's sampling theorem), and a small margin localises the zeros of the minimiser's transform.
* **Observations at λ = 2…6.** The detection identity 4 Re Φ(γ − iδ)² and the resolution height H(λ) = 2πλ²; reconstruction of the zeta zeros below H(λ) from the primes ≤ λ² (γ₁ to 105 digits from the primes up to 25); the minimiser equals Riemann's Φ times a smooth window; margins and rigidity of the prime data; the frame identity.

Nothing here bears on the Riemann Hypothesis beyond the finite scales treated.

## Reproducing the certificates

```
pip install mpmath python-flint numpy scipy
cd code
python3 certify_symbol.py 2 80 110 40 28     # Theorem 2, ~9 min on one core
python3 certify_rigidity.py 2 24             # Theorem 3, ~3.5 min
python3 certify.py 2 30                      # Theorem 1 (block), ~15 min at 200 digits
```
See `code/README.txt` for the other scripts and the results they produce.

## Licence

Code: MIT. Text (paper and archive document): CC BY 4.0.

## Citation

See `CITATION.cff`. Computations and drafting were assisted by Claude (Anthropic).
