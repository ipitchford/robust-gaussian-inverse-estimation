# Quantitative robust estimation for Gaussian inverse problems

Anonymous · 1.0.0-candidate · 18 October 2026

This manuscript proves the proposed spectral estimator's asymptotic optimality among linear estimators for Gaussian inverse problems under samplewise norm-bounded attacks. It gives a uniform finite-dimensional additive risk-regret bound, an exact anisotropic-jitter implementation, a sharp isotropic-jitter boundary, and a certified two-band separation.

The conclusion concerns normalized additive risk, not risk ratios or convergence of estimator matrices.

## Package map

- `paper.pdf`, `paper.tex`, `paper.md`: complete manuscript.
- `CLAIMS.json`: machine-readable claim and evidence map.
- `code/replay.py`: clean replay orchestrator.
- `code/inverse_risk_checks.py`: finite numerical checks using exact Gaussian norm expectations and scalar optimization.
- `code/analytic_isotropic_gap.py`: exact rational certificate for the two-band gap.
- `figures/inverse_jitter_illustration.png`: explanatory figure.
- `AI_INDEX.md`, `ASSURANCE.md`, `STATUS.md`: reading and assurance boundaries.

Run `python code/replay.py --output replay-output.json` from the repository root. Original prose, data and figures are CC0-1.0; original code is MIT licensed.

This is an unrefereed AI-assisted research candidate. Finite replay is not proof-assistant formalisation, independent reproduction or human peer review.
