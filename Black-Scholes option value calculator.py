
import numpy as np
import math

mu = 0
T = 1
S0 = 35
sigma = 0.4
X = 35
r = 0.01
q = 0.03
def norm_dist_cdf(x): 
    return 0.5 * (1 + math.erf(x / np.sqrt(2))) #erf(x) is an integral of the standard normal distribution
  
d1 = (np.log(S0 / X) + (r - q + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
d2 = d1 - sigma * np.sqrt(T)

put_option_price = X * np.exp(-r * T) * norm_dist_cdf(-d2) - S0* np.exp(-q * T) * norm_dist_cdf(-d1)
print(put_option_price)

