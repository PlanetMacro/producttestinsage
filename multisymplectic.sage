from sage.all import *
load("util.sage")
load("PoissonSolver.sage")

# Extended multiphase space P with coordinates (x1,x2,x3,x4,q,p1,p2,p3,p4,p)
P = Manifold(10, 'P', structure='smooth')
X.<q,p1,p2,p3,p4,p,x1,x2,x3,x4> = P.chart()

# Default frame and coframe
e = P.default_frame()
dq, dp1, dp2, dp3, dp4, dp, dx1, dx2, dx3, dx4 = e.coframe()

# Canonical multisymplectic (n+1)-form from math-ph/0202043:
# omega = dq ^ dp^mu ^ d^n x_mu - dp ^ d^n x   with mu = 1..4
d4x = dx1.wedge(dx2).wedge(dx3).wedge(dx4)
d3x_mu = [
    dx2.wedge(dx3).wedge(dx4),           # mu = 1
    -dx1.wedge(dx3).wedge(dx4),          # mu = 2
    dx1.wedge(dx2).wedge(dx4),           # mu = 3
    -dx1.wedge(dx2).wedge(dx3),          # mu = 4
]
dp_mu = [dp1, dp2, dp3, dp4]
omega = sum(dq.wedge(dp_mu[mu]).wedge(d3x_mu[mu]) for mu in range(4)) - dp.wedge(d4x)

print("Multisymplectic 5-form omega:")
omega = simplify_form_full(omega)
disp = str(omega.display())
print(disp)

# Poisson solver example: solve i_Y omega = q dx1∧dx2∧dx3
target = q * dx1.wedge(dx2).wedge(dx3)
Y_poisson = solve_poisson(omega, target, frame=e)
print("Poisson solver Y for i_Y omega = q dx1^dx2^dx3:")
print(Y_poisson.display())
check = cartan_interior(omega, Y_poisson)
print("Check i_Y omega (should match target):")
print(check.display())
