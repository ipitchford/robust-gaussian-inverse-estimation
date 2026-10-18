#!/usr/bin/env python3
"""Exact, standalone algebraic certificate for the two-band isotropic gap.

Run with Python 3; only the standard library is required.  No grid, stored
directions, numerical optimiser, or floating-point comparison is used.
The finite-dimensional conclusion uses the surrogate comparison proved in
the accompanying manuscript; this script checks its exact numerical inputs.
"""

from fractions import Fraction as F
import json


def add(a, b):
    """Polynomial addition; coefficients are in increasing degree order."""
    return [
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ]


def scale(a, k):
    return [k * x for x in a]


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def square(a):
    return multiply(a, a)


def evaluate(a, x):
    value = 0
    for coefficient in reversed(a):
        value = value * x + coefficient
    return value


def certificate():
    checks = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    s, p = [F(2, 5), F(1)], [F(1, 5), F(4, 5)]
    t, eps2 = F(1, 1000), F(1, 5)
    # sqrt(p_j)=sqrt_p_scaled[j]/sqrt(5), e_j=attack[j]/sqrt(5).
    sqrt_p_scaled, attack = [F(1), F(2)], [F(4, 5), F(3, 5)]
    check("band weights sum to one", sum(p) == 1)
    check("fixed attack has squared norm epsilon squared", sum(x*x for x in attack)/5 == eps2)
    check("scaled square roots of band weights", [x*x/5 for x in sqrt_p_scaled] == p)
    expected_numerators = [[F(8, 25), F(1)], [F(3, 5), F(2)]]
    for j in range(2):
        # h=s/(s^2+c): numerator of (1-sh)^2+t*h^2 is c^2+t*s^2.
        residual = add([s[j]**2, F(1)], [-s[j]**2])
        check(f"band {j+1}: residual variance numerator", add(square(residual), [t*s[j]**2]) == [t*s[j]**2, 0, 1])
        check(f"band {j+1}: dropping noise gives a lower bound", t*s[j]**2 > 0)
        check(f"band {j+1}: fixed attack numerator", [s[j]*attack[j], sqrt_p_scaled[j]] == expected_numerators[j])

    # y=25c gives L=1/5*((y+8)^2/(y+4)^2+(2y+15)^2/(y+25)^2).
    y_numerators, y_denominators = [[8, 1], [15, 2]], [[4, 1], [25, 1]]
    for j in range(2):
        check(f"band {j+1}: change of variable", [25*expected_numerators[j][0], expected_numerators[j][1]] == y_numerators[j] and [25*s[j]**2, 1] == y_denominators[j])
    denominator = multiply(square(y_denominators[0]), square(y_denominators[1]))
    numerator = add(multiply(square(y_numerators[0]), square(y_denominators[1])), multiply(square(y_numerators[1]), square(y_denominators[0])))
    polynomial = add(scale(numerator, 2000), scale(denominator, -5303))
    expected_polynomial = [34170000, 1162600, -1004423, 8426, 4697]
    check("L-5303/10000 has the claimed numerator", polynomial == expected_polynomial)
    quadratic = [288616, 92972, 4697]
    centred = add(add(multiply(square([-9, 1]), quadratic), scale([-9, 1], -1173044)), [234708])
    check("centred polynomial identity", centred == polynomial)
    check("quadratic factor positive and increasing on nonnegative inputs", all(x > 0 for x in quadratic))
    check("quadratic factor at nine", evaluate(quadratic, 9) == 1505821)
    check("positive constant for the interval zero to nine", 234708 > 0)
    discriminant = 1173044**2 - 4*1505821*234708
    check("claimed discriminant", discriminant == -37680715136)
    check("quadratic lower bound is strictly positive", 1505821 > 0 and discriminant < 0)
    # For y<=9 each term in the centred identity is nonnegative, with a
    # positive constant.  For y>=9 set z=y-9 and use the positive quadratic
    # 1505821*z^2-1173044*z+234708.  These establish L>5303/10000 everywhere.

    lam = F(104196875, 10**8)
    h = [F(416267115648, 10**12), F(982718129443, 10**12)]
    upper_value = lam * eps2
    for j in range(2):
        den = lam - h[j]**2
        check(f"upper witness band {j+1}: strict dual domain", den > 0)
        q = (1 - s[j]*h[j])**2 + t*h[j]**2
        upper_value += p[j]*q*lam/den
    upper_threshold = F(388918, 10**6)
    upper_margin = upper_threshold - upper_value
    expected_margin = F(3168557554921612509053869041613787853379731, 25868592782876951709628705995150319490966201600000)
    check("exact upper certificate margin", upper_margin == expected_margin)
    check("surrogate witness strictly below 0.388918", upper_margin > 0)
    limiting_gap = F(5303, 10000) - upper_threshold
    check("strict limiting gap threshold", limiting_gap == F(141382, 10**6))

    d = 20000
    fourth_root_upper, square_root_value = F(12409, 100000), F(1, 100)
    fourth_power_margin = fourth_root_upper**4 - F(128, 27*d)
    check("exact fourth-power margin", fourth_power_margin == F(192972343173947, 2700000000000000000000))
    check("strict fourth-root upper bound", fourth_power_margin > 0)
    check("exact square-root term", square_root_value > 0 and square_root_value**2 == F(2, d))
    eta_upper = fourth_root_upper + square_root_value
    check("finite-dimensional error threshold", eta_upper == F(13409, 100000))
    finite_gap = limiting_gap - eta_upper
    check("strict finite-dimensional gap threshold", finite_gap == F(7292, 10**6))
    check("integral band multiplicities", all((d*x).denominator == 1 for x in p))

    return {
        "status": "PASS",
        "method": "standard-library exact rational and polynomial arithmetic",
        "checks_passed": len(checks),
        "checks": checks,
        "normalisation": {"signal_energy": "1", "singular_values": list(map(str, s)), "band_weights": list(map(str, p)), "noise_ratio": str(t), "epsilon_squared": str(eps2)},
        "polynomial_coefficients_increasing_degree": polynomial,
        "discriminant": discriminant,
        "isotropic_surrogate_lower_bound": "0.5303",
        "unrestricted_surrogate_strict_upper_bound": "0.388918",
        "exact_upper_margin": str(upper_margin),
        "limiting_additive_gap": "> 0.141382",
        "dimension": d,
        "exact_fourth_power_margin": str(fourth_power_margin),
        "eta_strict_upper_bound": "0.13409",
        "finite_dimensional_additive_gap": "> 0.007292",
        "finite_dimensional_application": "Uses the manuscript's proved bound I_d-T_d >= Delta-eta_d at an isotropic minimiser with R<=1, Q+epsilon^2 B^2<=1 and q_i<=1/d.",
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
