from sage.all import *

load("util.sage")
load("PoissonSolver.sage")

# most simple symplectic space (x1, x2)
M = Manifold(2, "M", structure="smooth")
X.<x1, x2> = M.chart()

# Default frame and coframe
e = M.default_frame()
dx1, dx2 = e.coframe()

# Symplectic form omega = dx1 ^ dx2
omega = dx1.wedge(dx2)

print("Symplectic 2-form omega:")
omega = simplify_form_full(omega)
disp = str(omega.display())
print(disp)

# Vector field X = d/dx1
X1 = e[0]

# Contraction i_X omega using first index position
iX_omega = omega.contract(0, X1)
iX_omega = simplify_form_full(iX_omega)
print("Contraction i_X omega with X = d/dx1 (first-slot convention):")
disp_contr = str(iX_omega.display())
print(disp_contr)

# Bivector X1 ^ X2 and contraction on first slot successively
X2 = e[1]
iX1X2_omega = omega.contract(0, X1).contract(0, X2)
print("Contraction i_{X1^X2} omega with X1=d/dx1, X2=d/dx2 (first-slot convention):")
print(iX1X2_omega.display())

# Built-in contraction with a multivector via interior_product (multivector ⟂ form)
X1wX2 = X1.wedge(X2)
iX1_builtin = X1.interior_product(omega)
print("Built-in interior_product i_{X1} omega (multivector method):")
print(iX1_builtin.display())
iX1wX2_builtin = X1wX2.interior_product(omega)
print("Built-in interior_product i_{X1^X2} omega (multivector method, Sage's antisymmetrization, k! factor):")
print(iX1wX2_builtin.display())

# Normalized Cartan right interior product for multivectors
iX1wX2_cartan = cartan_interior(omega, X1wX2)
print("Cartan right interior i_{X1^X2} omega (normalized, matches i_{X2} i_{X1}):")
print(iX1wX2_cartan.display())

# Poisson solver example: solve i_Y omega = dx2
Y_poisson = solve_poisson(omega, dx2, frame=e)
print("Poisson solver Y for i_Y omega = dx2:")
print(Y_poisson.display())
check_poisson = cartan_interior(omega, Y_poisson)
print("Check i_Y omega with solver:")
print(check_poisson.display())

# Hamilton solver example: f = x1^2, solve i_X omega = d f
f_scalar = x1**2
X_ham = solve_hamilton(omega, f_scalar, frame=e)
print("Hamilton solver X for i_X omega = d(x1^2):")
print(X_ham.display())
print("Check i_X omega (should be d(x1^2) = 2*x1*dx1):")
print(cartan_interior(omega, X_ham).display())
