def associated_legendre(l, m, x, mem={}):
    # Immediate result based on parameter m 
    if m < 0 or m > l:
        return 0
    
    if m == 0:
        return Polynome_legendre(l, x, mem)
    
    # Check memory for saved results
    if (l, m, x) in mem:
        return mem[(l, m, x)]
    
    # Using the recurrence relation https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Recurrence_formula
    if l == m:
        result = (-1)**m * (1 - x**2)**(m/2) * double_factorial(2*m-1)
    else:
        if l-1 not in mem:
            mem[(l-1, m, x)] = associated_legendre(l-1, m, x, mem)
        if l-2 not in mem:
            mem[(l-2, m, x)] = associated_legendre(l-2, m, x, mem)
        
	# Calculate result
        result = ((2*l - 1)*x*associated_legendre(l-1, m, x, mem) - (l + m - 1)*associated_legendre(l-2, m, x, mem)) / (l - m)

    # Save current result in memory

    mem[(l, m, x)] = result
    
    # Function output
    return result

def double_factorial(n):
    if n == 0 or n == -1:
        return 1
    else:
        return n * double_factorial(n-2)


#import matplotlib and numpy for visualization
#valeurs_x = np.linspace(-1.0, 1.0, 1000)
#plt.plot(valeurs_x, [associated_legendre(2,2,x) for x in valeurs_x], label=f'P{0}(x)')
#plt.plot(valeurs_x, [associated_legendre(4,2,x) for x in valeurs_x], label=f'P{1}(x)')
#plt.plot(valeurs_x, [associated_legendre(6,2,x) for x in valeurs_x], label=f'P{2}(x)')
#plt.plot(valeurs_x, [associated_legendre(3,2,x) for x in valeurs_x], label=f'P{3}(x)')
#plt.plot(valeurs_x, [associated_legendre(5,2,x) for x in valeurs_x], label=f'P{4}(x)')
#plt.legend()
#plt.show()
