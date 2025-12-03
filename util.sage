from sage.all import *

def simplify_form_full(alpha):
    """
    Return a copy of the differential form with all scalar components simplified
    via ``simplify_full``.
    """
    beta = alpha.copy()
    beta.apply_map(lambda c: c.simplify_full())
    return beta

__all__ = ["simplify_form_full"]
