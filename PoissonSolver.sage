from sage.all import *
from itertools import combinations

# Optional helper from util, if present
try:
    load("util.sage")
except Exception:
    pass


def cartan_interior(alpha, Q):
    """
    Cartan-style right interior product i_Q alpha for a decomposable multivector Q.
    Sage's multivector interior_product returns k! times this for a k-vector, so we normalize.
    """
    try:
        k = Q.degree()
    except AttributeError:
        raise TypeError("Q must be a multivector field with a degree() method")
    res = Q.interior_product(alpha)
    if k <= 1:
        return res
    return (1 / factorial(k)) * res


def _scalar_expr(sf, chart):
    """Return a symbolic expression for a scalar field on a given chart."""
    if hasattr(sf, "expr"):
        try:
            return sf.expr(chart)
        except Exception:
            try:
                return sf.expr()
            except Exception:
                pass
    return SR(sf)


def solve_poisson(omega, f, frame=None):
    """
    Solve i_Y omega = f for a multivector field Y.

    omega -- k-form
    f     -- l-form with l <= k
    frame -- (optional) frame to use for coordinates; defaults to omega.domain().default_frame()

    Returns: a multivector of degree k-l that satisfies the equation if a (least-squares/pseudo)
    solution exists. If the system is inconsistent, raises a RuntimeError.
    """
    M = omega.domain()
    if frame is None:
        frame = M.default_frame()
    chart = M.default_chart()

    k = omega.degree()
    try:
        l = f.degree()
    except Exception:
        try:
            f = omega.domain().scalar_field(f)
            l = 0
        except Exception as exc:
            raise TypeError("f must be a differential form or scalar field") from exc
    if l > k:
        raise ValueError("Need l <= k to solve i_Y omega = f")
    r = k - l
    if r == 0:
        # Trivial case: need omega == f
        diff = omega - f
        try:
            if diff.is_zero():
                return M.scalar_field(1)
        except Exception:
            pass
        raise RuntimeError("No scalar solution: omega and f differ")

    try:
        n = len(frame)
    except Exception:
        n = M.dimension()
    if r > n:
        raise RuntimeError("Degree mismatch: cannot build an r-vector with r > dimension")

    # Build basis of r-vectors and l-form evaluation slots
    r_combos = list(combinations(range(n), r))
    l_combos = list(combinations(range(n), l))
    basis_r = []

    # Build contraction matrix A (rows: l-form components, cols: r-vector basis elements)
    A_entries = []
    alpha_list = []
    for rc in r_combos:
        mv = frame[rc[0]]
        for idx in rc[1:]:
            mv = mv.wedge(frame[idx])
        basis_r.append(mv)
        alpha_list.append(cartan_interior(omega, mv))

    for lc in l_combos:
        for alpha in alpha_list:
            if l == 0:
                sf = alpha
            else:
                sf = alpha(*[frame[i] for i in lc])
            A_entries.append(_scalar_expr(sf, chart))

    A = matrix(SR, len(l_combos), len(r_combos), A_entries)

    # Build RHS b from f
    b_entries = []
    for lc in l_combos:
        if l == 0:
            sf = f
        else:
            sf = f(*[frame[i] for i in lc])
        b_entries.append(_scalar_expr(sf, chart))
    b = vector(SR, b_entries)

    # Solve the linear system.
    sol = None
    pivot_cols = A.pivots()
    if pivot_cols:
        Ap = A.matrix_from_columns(pivot_cols)
        try:
            ATA = Ap.transpose() * Ap
            rhs = Ap.transpose() * b
            coords = ATA.solve_right(rhs)
            full = [0] * A.ncols()
            for idx, c in zip(pivot_cols, coords):
                full[idx] = c
            sol = vector(SR, full)
        except Exception:
            sol = None
    if sol is None:
        try:
            sol = A.solve_right(b)
        except Exception:
            pass
    if sol is None:
        try:
            sol = A.pseudoinverse() * b
        except Exception as exc:
            raise RuntimeError("No solution found or system inconsistent") from exc

    # Residual check (raise if clearly inconsistent)
    residual = A * sol - b
    if not all(r == 0 for r in residual):
        raise RuntimeError("System inconsistent; residual not zero")

    # Assemble multivector Y
    Y = basis_r[0].parent().zero()
    for coeff, mv in zip(sol, basis_r):
        try:
            zero_coeff = coeff.is_zero()
        except Exception:
            zero_coeff = False
        if not zero_coeff:
            Y += coeff * mv
    return Y


def solve_hamilton(omega, f, frame=None):
    """
    Solve i_Y omega = d f for a multivector Y by wrapping solve_poisson.

    f can be a scalar field (0-form) or a general differential form.
    """
    if hasattr(f, "exterior_derivative"):
        df = f.exterior_derivative()
    else:
        try:
            M = omega.domain()
            sf = M.scalar_field(f)
            df = sf.exterior_derivative()
        except Exception as exc:
            raise TypeError("f must be a differential form or scalar field with exterior_derivative") from exc
    return solve_poisson(omega, df, frame=frame)


def is_poisson(omega, f, frame=None):
    """
    Attempt to solve i_Y omega = f and i_X omega = d f.
    Returns (True, Y, X) if both exist, else (False, None, None).
    """
    try:
        Y = solve_poisson(omega, f, frame=frame)
    except Exception:
        return False, None, None
    try:
        X = solve_hamilton(omega, f, frame=frame)
    except Exception:
        return False, None, None
    return True, Y, X
