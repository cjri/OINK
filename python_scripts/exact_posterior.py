
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import exp

def read_data(file_path, type_first=int):
    data = []
    with open(file_path, 'r') as f:
        try:
            while True:
                s = f.readline().strip()
                if not s :
                    break
                s = s.split(' ')
                data.append((type_first(s[0]), float(s[1])))
        except IOError:
            pass
    return data

prob_no_end = np.array(read_data("julia_R0_acceptance_rates.dat", type_first=float))[:,1]
prob_no_end = prob_no_end/np.sum(prob_no_end)


d = 0.1

R0_vals = np.linspace(0.1, 4, 40)
accept = []

for R0 in R0_vals:

    p0 = fsolve(lambda p: p - (1-d)*exp(-R0)*exp(p*R0), 0)
    p1 = d*p0 / (1 - (1-d)*R0*exp(-R0)*exp(R0*p0))

    p1_hat = p1*R0*exp(-R0)*exp(p0*R0)

    accept.append(p1_hat)


accept = np.array(accept)
accept = accept/np.sum(accept)
plt.figure()
plt.plot(R0_vals, accept)
plt.plot(R0_vals, prob_no_end)
plt.show()
