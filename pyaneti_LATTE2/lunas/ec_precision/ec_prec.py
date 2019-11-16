import numpy as np

se = 0.05
K  = 10. # m/s
so = 5.  # m/s

N = - ( np.log10(se) - 0.48 ) / 0.89
N = so / K * 10**(N)
N = N*N

print N
