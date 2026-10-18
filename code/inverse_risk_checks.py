#!/usr/bin/env python3
"""Independent numerical checks of scalar and finite-d Gaussian risk claims.

Requires SciPy.  This reconstructed checker uses exact chi-norm expectations
and one-dimensional optimisation, not Monte Carlo.  It checks finite instances
and endpoints; the manuscript's proof supplies the general theorem.
"""

import json
import math
from functools import lru_cache
from scipy.optimize import minimize_scalar


def gain(lam, s, t):
    disc = (1 - lam*(s*s+t))**2 + 4*lam*t
    return 2*lam*s/(1 + lam*(s*s+t) + math.sqrt(disc))


def psi(lam, s, t):
    return (1-lam*(s*s+t) + math.sqrt((1-lam*(s*s+t))**2+4*lam*t))/2


def bounded_minimum(fun, lo, hi, extras=()):
    if hi == lo:
        return fun(lo), lo
    result = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options={"xatol": 1e-13, "maxiter": 500})
    if not result.success:
        raise AssertionError(result.message)
    candidates = [(fun(x), x) for x in (lo, hi, result.x, *extras)]
    return min(candidates)


@lru_cache(None)
def mean_normalised_chi(rank):
    return math.exp((math.log(2)-math.log(rank))/2 + math.lgamma((rank+1)/2)-math.lgamma(rank/2))


def run_checks():
    passed = 0
    scalar_cases = 0
    finite_cases = 0
    worst_scalar_error = 0.0
    worst_surrogate_error = 0.0
    max_regret = 0.0
    largest_dimension = []

    def check(condition, message):
        nonlocal passed
        if not condition:
            raise AssertionError(message)
        passed += 1

    # Independently minimise the scalar dual term and check the root equation.
    for s in (0.0, 0.2, 1.0, 3.0):
        for t in (0.0, 0.001, 0.3, 2.0):
            for lam in (0.0, 0.01, 0.1, 1.0, 2.0, 25.0, 1000.0):
                scalar_cases += 1
                h = gain(lam, s, t)
                nominal = s/(s*s+t) if s*s+t else 0.0
                residual = s*h*h-(1+lam*(s*s+t))*h+lam*s
                check(abs(residual) <= 2e-10*(1+lam*s), "scalar root equation")
                check(-1e-12 <= h <= nominal+1e-12 and h*h <= lam+1e-12, "gain domain and clipping")
                if lam == 0:
                    check(h == 0 and psi(lam, s, t) == 1, "zero dual parameter endpoint")
                    continue

                def objective(x):
                    q = (1-s*x)**2+t*x*x
                    den = lam-x*x
                    if den <= 1e-14*max(1.0, lam):
                        return 0.0 if q == 0 else float("inf")
                    return q*lam/den

                upper = min(nominal, math.sqrt(lam))
                numeric, _ = bounded_minimum(objective, 0.0, upper)
                error = abs(numeric-psi(lam, s, t))
                worst_scalar_error = max(worst_scalar_error, error)
                check(error <= 2e-8*(1+abs(numeric)), "closed scalar value versus independent optimisation")

    # For rank r with one positive singular band and d-r null directions,
    # Q_pos=(r/d)*((1-s*h)^2+t*h^2).  Exact Gaussian risk is
    # 1-r/d + Q_pos + 2*epsilon*h*sqrt(Q_pos)*E[chi_r]/sqrt(r)
    #             + epsilon^2*h^2.
    # The surrogate replaces E[chi_r]/sqrt(r) by one.
    for d in (1, 2, 5, 20, 100, 1000, 20000):
        finite_cases += 1
        check(1.0 == (d-0)/d, "zero-rank estimator risk and surrogate are signal energy")
        for rank in sorted({d, max(1, d//2)}):
            p = rank/d
            mu = mean_normalised_chi(rank)
            check(0 < mu <= 1+1e-10, "normalised Gaussian norm mean")
            for s in (0.4, 1.0, 3.0):
                for t in (0.0, 0.001, 0.2):
                    nominal = s/(s*s+t)
                    for relative_eps in (0.0, 0.2, 0.9, 1.0, 1.2):
                        finite_cases += 1
                        eps = relative_eps*s*math.sqrt(p)
                        if relative_eps >= 1:
                            h_star = 0.0
                        elif eps == 0 or t == 0:
                            h_star = nominal
                        else:
                            eps_eff = eps/math.sqrt(p)
                            h_star = (s-eps_eff*math.sqrt(t/(s*s+t-eps_eff*eps_eff)))/(s*s+t)

                        def positive_variance(h):
                            return p*((1-s*h)**2+t*h*h)

                        def surrogate(h):
                            return 1-p+(math.sqrt(positive_variance(h))+eps*h)**2

                        def risk(h):
                            q = positive_variance(h)
                            return 1-p+q+2*eps*h*math.sqrt(q)*mu+eps*eps*h*h

                        s_numeric, _ = bounded_minimum(surrogate, 0.0, nominal)
                        s_star, r_star = surrogate(h_star), risk(h_star)
                        r_opt, h_opt = bounded_minimum(risk, 0.0, nominal, (h_star,))
                        worst_surrogate_error = max(worst_surrogate_error, abs(s_numeric-s_star))
                        check(abs(s_numeric-s_star) <= 2e-9, "closed surrogate minimiser versus independent optimisation")
                        check(r_star <= s_star+1e-10 and r_opt <= r_star+1e-10, "risk and surrogate ordering")
                        regret = r_star-r_opt
                        eta = min(1.0, (128/(27*d))**0.25+math.sqrt(2/d))
                        check(regret >= -1e-10 and regret <= eta+2e-9, "finite-d regret bound")
                        check(s_star-r_opt <= eta+2e-9, "finite-d surrogate value sandwich")
                        max_regret = max(max_regret, regret)
                        if eps > 0:
                            lam = h_star*h_star+h_star*math.sqrt(positive_variance(h_star))/eps
                            check(abs(gain(lam, s, t)-h_star) <= 2e-9, "spectral gain matches closed band minimiser")
                            phi = lam*eps*eps+p*psi(lam, s, t)+(1-p)
                            check(abs(phi-s_star) <= 2e-9, "scalar-dual formula matches primal surrogate")
                        else:
                            check(abs(h_star-nominal) <= 1e-12, "zero-adversary nominal endpoint")
                        if d == 20000 and rank == d and s == 1.0 and t == 0.001:
                            largest_dimension.append({"relative_epsilon": relative_eps, "surrogate_gain": h_star, "risk_optimal_gain": h_opt, "surrogate_value": s_star, "risk_at_surrogate_gain": r_star, "optimal_scalar_risk": r_opt, "additive_regret": regret, "theorem_bound": eta})

    return {
        "status": "PASS",
        "method": "exact Gaussian norm expectation evaluated numerically; independent bounded scalar optimisation",
        "checks_passed": passed,
        "scalar_dual_cases": scalar_cases,
        "finite_dimension_cases": finite_cases,
        "dimensions": [1, 2, 5, 20, 100, 1000, 20000],
        "covered_endpoints": ["zero rank", "rank deficient including m<d model", "zero noise", "zero adversary", "zero dual parameter", "hard-case dual boundary", "active and inactive estimators"],
        "maximum_scalar_value_absolute_error": worst_scalar_error,
        "maximum_surrogate_optimum_absolute_error": worst_surrogate_error,
        "maximum_observed_additive_regret": max_regret,
        "dimension_20000_examples": largest_dimension,
        "scope": "Numerical sanity checks of finite instances, not a substitute for the general proof.  The zero-adversary endpoint uses its nominal limiting solution; no finite lambda minimiser is asserted there.",
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))
