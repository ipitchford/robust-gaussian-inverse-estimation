# AI reading index

Start with [STATUS.md](STATUS.md), then read [CLAIMS.json](CLAIMS.json) and the complete [paper](paper.pdf) or its [LaTeX source](paper.tex).

## Principal results

- Uniform normalized additive regret at most `min{1,(128/(27d))^(1/4)+sqrt(2/d)}` among linear estimators.
- Exact anisotropic population-jitter implementation of the spectral rule when the minimizing dual parameter is positive.
- Under positive noise and active estimation, finite isotropic jitter implements the rule exactly iff the positive singular spectrum is flat.
- A two-band example has limiting gap greater than `0.141382` and finite-dimensional gap greater than `0.007292` at `d=20000`.

## Evidence routes

- General theorem: analytic proof in [paper.tex](paper.tex).
- Finite numerical stress checks: [inverse_risk_checks.py](code/inverse_risk_checks.py).
- Exact two-band certificate: [analytic_isotropic_gap.py](code/analytic_isotropic_gap.py).
- Scope and exclusions: [CLAIMS.json](CLAIMS.json), [ASSURANCE.md](ASSURANCE.md), and the manuscript's claim-to-evidence table.

Do not infer a risk-ratio theorem, matrix convergence, nonlinear-estimator optimality, a sharper convergence rate, formal verification, independent reproduction or human peer review.

For replay and provenance, see [REPLAY.md](REPLAY.md), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [CITATION_AUDIT.md](CITATION_AUDIT.md) and [PROVENANCE.md](PROVENANCE.md).
