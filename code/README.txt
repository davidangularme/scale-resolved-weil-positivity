v33 code bundle - Scale-Resolved Weil Positivity (F. D. Blum, October 2026)
Requirements: python3, mpmath, python-flint (Arb), numpy, scipy.
All scripts are prime-side: no zero of zeta enters the construction of W_lambda; zeros*.txt are used only for comparison.

Scripts (section of the paper):
  weil_prime.py   prime-side Weil matrix in the cosine basis (2)
  weil_gram.py, checks.py, make_points.py   zero-side matrix, validation, control point sets (2, 3)
  detect.py       displaced-zero detection experiments (3)
  recon.py        zeros of the minimiser transform vs true zeros, range law (4)
  phi_test.py, window.py, factor.py   Phi-window structure of the minimiser (5)
  margin.py, rigidity.py, jacobian.py   L2 margin, prime tolerances, first-order response (6)
  certify.py      ball-arithmetic block certificate (Theorems 1, 1') (7)
  prolate.py, prolate_cross.py   prolate eigenvalues; prime shift operator in the prolate basis (7)
  symbol_B.py, certify_symbol.py   envelope form B in the Legendre basis; Arb certificate of Theorem 2 (7)
                  run: python3 certify_symbol.py 2 80 110 40 28   (9 min)
  certify_rigidity.py   eight indefiniteness certificates of Theorem 3 (7)
                  run: python3 certify_rigidity.py 2 24           (3.5 min)
  s_eff.py        compressed prime mass (7)
Results in results/. Precision 2.6N+80 digits unless stated.
Related work found after completion: X. Zhu, arXiv:2608.24827 (same reduction, certified positivity for L <= 0.8).
