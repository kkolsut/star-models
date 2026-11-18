"""Reference values of the generalised Fermi-Dirac integrals for the tests of
crates/stellar-core/src/physics/fermi_dirac.rs.

F_k(eta, theta) = int_0^inf x^k sqrt(1 + theta x / 2) / (exp(x - eta) + 1) dx

computed with mpmath at 40 significant digits (tanh-sinh quadrature with break
points around x = eta). Run:  python3 tools/reference/fermi_dirac_reference.py > \
    crates/stellar-core/src/physics/fermi_dirac_reference.in
"""
import mpmath as mp

mp.mp.dps = 40


def fd_with(k, eta, theta, pts):
    # mp.quad works to an absolute tolerance: integrate F / scale with a scale
    # of the size of F (e^eta for eta < 0, eta^(k+1) above), and x = u^2,
    # which removes the x^(1/2) behaviour at the origin.
    scale = mp.exp(eta) if eta < 0 else (1 + eta) ** (k + 1)
    if eta < 0:
        f = lambda u: 2 * u ** (2 * k + 1) * mp.sqrt(1 + theta * u * u / 2) * mp.exp(-u * u) / (1 + mp.exp(eta - u * u))
    else:
        f = lambda u: 2 * u ** (2 * k + 1) * mp.sqrt(1 + theta * u * u / 2) / (mp.exp(u * u - eta) + 1) / scale
    upts = [mp.sqrt(v) if v != mp.inf else mp.inf for v in pts]
    return mp.quad(f, upts, maxdegree=12) * (scale if eta < 0 else scale)


def fd(k, eta, theta):
    """Two independent sets of break points; they must agree to 1e-16."""
    eta = mp.mpf(eta)
    theta = mp.mpf(theta)
    e = max(eta, 0)
    a = [0, 1, 5, 20] + ([eta / 2, eta - 10, eta, eta + 10] if eta > 25 else []) + [e + 60, mp.inf]
    b = [0, mp.mpf('0.5'), 3, 12] + ([eta * mp.mpf('0.7'), eta - 5, eta + 3, eta + 30] if eta > 25 else []) + [e + 90, mp.inf]
    a = sorted(set(a), key=lambda v: float(v) if v != mp.inf else 1e300)
    b = sorted(set(b), key=lambda v: float(v) if v != mp.inf else 1e300)
    va, vb = fd_with(k, eta, theta, a), fd_with(k, eta, theta, b)
    assert abs(va / vb - 1) < mp.mpf('1e-16'), (k, eta, theta, va, vb)
    return va


etas = [-100, -60, -30, -10, -3, -1, -0.3, 0, 0.3, 1, 2, 5, 10, 20, 50, 100, 300, 1000, 3000, 10000]
thetas = [1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0]
rows = []
for eta in etas:
    for th in thetas:
        xb = float(mp.log(1 + mp.exp(eta)))
        if th * xb > 10:
            continue
        vals = [fd(mp.mpf(k) / 2, eta, th) for k in (1, 3, 5)]
        rows.append("    ({:.17e}, {:.17e}, [{}]),".format(eta, th, ", ".join(mp.nstr(v, 20, min_fixed=0, max_fixed=0) for v in vals)))
print("[")
print("\n".join(rows))
print("]")
