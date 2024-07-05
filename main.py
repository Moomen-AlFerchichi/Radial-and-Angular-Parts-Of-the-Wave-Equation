import numpy as np
import matplotlib.pyplot as plt
from math import factorial
#----------------------------------------------------------------------------
# Radial Part of psi :
#----------------------------------------------------------------------------
# Define The associated Laguerre polynomial Using Reccursion Method
def Laguerre_Associated(k, x, a, mem={}):
    # Define index 0 case
    if k == 0:
        return 1
    # Define index 1 case
    elif k == 1:
        return -x + a + 1
    # Import result if stored in memory
    if (k, x, a) in mem:
        return mem[(k, x, a)]
    # Calculating results based on https://en.wikipedia.org/wiki/Laguerre_polynomials#Recurrence_relations
    result = ((2 * k - 1 + a - x) * Laguerre_Associated(k - 1, x, a, mem) - (k - 1 + a) * Laguerre_Associated(k - 2, x, a, mem)) / k
    # Save the result for a specific k,x and a
    # Important for not exceeding the maximum reccursion depth
    mem[(k, x, a)] = result
    return result

# Define the radial function
def fonction_radial(n, l, r, a0):
    # Decomposing the formula elements
    p = 2 * r / (n * a0)
    norm_const = np.sqrt((2 / (n * a0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    # Returning result for a given n,l,a0 parameters and r
    return norm_const*r**l * np.exp(-p / 2) * Laguerre_Associated(n - l - 1, p, 2 * l + 1)

#----------------------------------------------------------------------------
# Angular part of psi : 
#----------------------------------------------------------------------------
def Polynome_legendre(n, x, mem={}):
    if n == 0:
        return 1
    elif n == 1:
        return x
    if (n, x) in mem:
        return mem[(n, x)]
    result = ((2*n-1)*x*Polynome_legendre(n-1, x, mem) - (n-1)*Polynome_legendre(n-2, x, mem)) / n
    mem[(n, x)] = result
    return result

# Same principle as the generelized_Laguerre function
def associated_legendre(m, l, x, mem={}):
    if m < 0 or m > l:
        return 0
    if m == 0:
        return Polynome_legendre(l, x, mem)
    if (l, m, x) in mem:
        return mem[(l, m, x)]
    if l == m:
        result = (-1)**m * (1 - x**2)**(m/2) * double_factorial(2*m-1)
    else:
        if l-1 not in mem:
            mem[(l-1, m, x)] = associated_legendre(m, l-1, x, mem)
        if l-2 not in mem:
            mem[(l-2, m, x)] = associated_legendre(m, l-2, x, mem)
        # Calculating results based on https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Recurrence_formula
        result = ((2 * l - 1)*x*associated_legendre(m, l-1, x, mem) - (l + m - 1)*associated_legendre(m, l - 2, x, mem)) / (l - m)
    mem[(l, m, x)] = result
    return result

# Double factorial function using reccursion
def double_factorial(n):
    if n == 0 or n == -1:
        return 1
    else:
        return n * double_factorial(n-2)

# Spherical Harmonic
def fonction_angulaire(m, l, th, phi):
    # Normalization constant
    norm_const = (-1)**m * np.sqrt((2*l+1) * factorial(l-abs(m)) / (4 * np.pi * factorial(l+abs(m))))
    # Value of the associated_Legendre must be vectorized for proper ploting
    return np.vectorize(associated_legendre)(m, l, np.cos(th)) * norm_const * np.real(np.exp(1.j * m * phi))

# Simulation Parameters
print("Make sure that n!=0 and n<l<m | Soiez sure que n!=0 et n<l<=m")
print("Errors can accumulate at higher parameters | des erreurs peut s'accumuler au parametres élevés")
n = int(input("n (principal quantum number)|(nombre quantique principal): "))
l = int(input("l (azimuthal quantum number)|(nombre quantique azimutal) : "))
m = int(input("magnetic quantum number|(nombre quantique magnétique) : "))
a0 = 5.291772107e-11 

# Angular Function
r_values = np.linspace(0.0, 1e-9, 100) 
phi = np.linspace(0, 2 * np.pi, 100)
theta = np.linspace(0, np.pi, 100)
phi, theta = np.meshgrid(phi, theta)
angular_function = np.abs(fonction_angulaire(m, l, theta, phi))

#projection of coordiantes in Cartesian x,y and z 
x = angular_function * np.sin(theta) * np.cos(phi)
y = angular_function * np.sin(theta) * np.sin(phi)
z = angular_function * np.cos(theta)

# Create the figure
plt.style.use("dark_background")
fig = plt.figure(figsize=(14, 6))

# Plot 2D Radial Function
ax1 = fig.add_subplot(121)
R_nl = np.array([fonction_radial(n, l, r, a0) for r in r_values])
plt.plot(r_values, R_nl)
ax1.set_title(f'2D Radial Function $R_{{{n}}}^{{{l}}}(r)$')
ax1.set_xlabel('r')
ax1.set_ylabel(f'$R_{{{n}}}^{{{l}}}(r)$')


# Plot 3D Angular Function
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(x, y, z, cmap='magma')
ax2.set_title(f'3D Angular Function $Y_{{{l}}}^{{{m}}}(\\theta, \\phi)$\n-colors have no scientific meaning-')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

# Show the plots
plt.tight_layout()
plt.show()