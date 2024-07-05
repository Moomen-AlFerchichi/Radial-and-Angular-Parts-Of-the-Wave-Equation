# Associated Laguerre polynomial
def Polynomial_recursive_L_k_alpha(k, x, a, mem={}):
    # case k=0 and k=1
    if k == 0:
        return 1
    elif k == 1:
        return -x + a + 1
    if (k, x) in mem:
        return mem[(k, x)]
    # reccursion formula @https://en.wikipedia.org/wiki/Laguerre_polynomials#Recurrence_relations substitute n with n-1

    result = ((2 * k - 1 + a - x) * Polynomial_recursive_L_k_alpha(k - 1, x, a) - (k - 1 + a) * Polynomial_recursive_L_k_alpha(k - 2, x, a)) / k

    #update memory	

    mem[(k,x)]=result
    
    #function output
    return result
#import matplotlib and numpy for visualization
#x_values = np.linspace(-3.0, 13.0, 100)
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 0) for x in  x_values],label=f'L_{1}(x)/alpha{0}')
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 1) for x in  x_values],label=f'L_{2}(x)/alpha{1}')
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 2) for x in  x_values],label=f'L_{3}(x)/alpha{2}')
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 3) for x in  x_values],label=f'L_{4}(x)/alpha{3}')
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 4) for x in  x_values],label=f'L_{5}(x)/alpha{4}')
#plt.plot(x_values,[Polynomial_recursive_L_k_alpha(3, x, 5) for x in  x_values],label=f'L_{6}(x)/alpha{5}')
#plt.grid(True)
#plt.legend()
#plt.xlim((-3.0, 13.0))
#plt.ylim((-5.0,10.0))
#plt.show()
