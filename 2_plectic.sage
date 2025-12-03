from sage.all import *

load("util.sage")
load("PoissonSolver.sage")

# 2-plectic form on a 3D manifold with coordinates (x1, x2, x3)
M = Manifold(3, "M", structure="smooth")
X.<x1, x2, x3> = M.chart()

# Default frame and coframe
e = M.default_frame()
dx1, dx2, dx3 = e.coframe()

# 2-plectic form omega = dx1 ^ dx2
omega = dx1.wedge(dx2)

print("2-plectic 2-form omega on a 3-manifold:")
omega = simplify_form_full(omega)
disp = str(omega.display())
print(disp)

# Vector field X1 = d/dx1
X1 = e[0]

# Contraction i_X omega using first index position
iX_omega = omega.contract(0, X1)
iX_omega = simplify_form_full(iX_omega)
print("Contraction i_X omega with X = d/dx1 (first-slot convention):")
disp_contr = str(iX_omega.display())
print(disp_contr)

# Bivector X1 ^ X2 and contraction on first slot successively
X2 = e[1]
X3 = e[2]
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
f1 = dx2
Y_poisson = solve_poisson(omega, f1, frame=e)
print("Poisson solver Y for i_Y omega = dx2:")
print(Y_poisson.display())
check_poisson = cartan_interior(omega, Y_poisson)
print("Check i_Y omega with solver (should be dx2):")
print(check_poisson.display())

# Hamilton solver example: abstract scalar field f(x1, x2) and solve i_X omega = d f
f_func = function("f")(x1, x2)
f_scalar = M.scalar_field(f_func, name="f")
X_ham = solve_hamilton(omega, f_scalar, frame=e)
print("Hamilton solver X for i_X omega = d f (f depends on x1, x2):")
print(X_ham.display())
print("Check i_X omega (should be d f):")
print(cartan_interior(omega, X_ham).display())

# is_poisson helper: check existence on Poisson examples
ok_f1, Y_f1, X_f1 = is_poisson(omega, f1, frame=e)
print("is_poisson for f1 = dx2 (expect True):", ok_f1)
if ok_f1:
    print("Y (i_Y omega = f1):", Y_f1.display())
    print("X (i_X omega = d f1):", X_f1.display())

ok_f2, Y_f2, X_f2 = is_poisson(omega, f_scalar, frame=e)
print("is_poisson for f2 = f (abstract 0-form, expect True):", ok_f2)
if ok_f2:
    print("Y (i_Y omega = f2):", Y_f2.display())
    print("X (i_X omega = d f2):", X_f2.display())

# Non-Poisson examples
g1 = dx3
ok_g1, Y_g1, X_g1 = is_poisson(omega, g1, frame=e)
print("is_poisson for g1 = dx3 (expect False):", ok_g1)

g2 = dx1.wedge(dx3)
ok_g2, Y_g2, X_g2 = is_poisson(omega, g2, frame=e)
print("is_poisson for g2 = dx1^dx3 (expect False):", ok_g2)
