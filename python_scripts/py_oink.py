
import numpy as np
from random import random
from scipy.stats import poisson, weibull_min
import matplotlib.pyplot as plt
import ray

ray.init()

"""
gsl weibull
p(x) dx = {b \over a^b} x^{b-1}  \exp(-(x/a)^b) dx
double gsl_ran_weibull(const gsl_rng *r, double a, double b)
 
scipy weibull_min 

weibull_min.pdf(x, c, loc, scale)
"""


def simulate():
    T_sim = 10000
    a1 = 7.4026
    b1 = 1.7375
    a2 = 1.0314
    b2 = 1.0025
    time_to_detection = 18

    p_detect = 0.1
    max_detected = 1

    max_cases = 2000
    
    ## Add first case
    def rng_symptom():
    # gsl_ran_weibull(rgen,p.incubation_b,p.incubation_a)
        return weibull_min.rvs(a1, scale=b1)
       
    def rng_infect():
    # gsl_ran_weibull(rgen,p.infection_b,p.infection_a)
        return weibull_min.rvs(a2, scale=b2)

    @ray.remote
    def sample(R0=3, N=1000):
        result = []
        for i in range(N):
            n_detected = 0
            n_cases = 0
            first_detection_time = None

            population = [ (0.0, 0.0 + rng_symptom(), False) ]
            old_population = []
            p = population.pop()
            time = p[0]
            old_population.append(p)
            n_cases += 1

            fail = False
            while time <= T_sim and n_cases<=max_cases:
                n_infect = poisson.rvs(R0)
                for j in range(n_infect):
                    n_cases += 1
                    t_infect = p[1] + rng_infect()
                    t_symptom = t_infect + rng_symptom()
                    detected = random() < p_detect
                    population.append((t_infect, t_symptom, detected))
                    if detected:
                        t_detect = t_symptom + time_to_detection
                        n_detected += 1
                        if first_detection_time is None or t_detect < first_detection_time:
                             first_detection_time = t_detect
                        if n_detected > max_detected:
                             fail=True
                             break

                if not population or fail:
                    break
                population.sort()
                p = population.pop(0)
                old_population.append(p)
                time = p[0]
            result.append((fail, first_detection_time, time, n_cases, n_detected))
        return result

    R0_vals = np.linspace(0.1, 4.0, 40)
    results = [sample.remote(R0=R0, N=1000) for R0 in R0_vals ] 
    results = ray.get(results)
    

        
    plt.figure()
    plt.plot(R0_vals, [sum([(not s[0] and s[4]==1) for s in r]) for r in results])
    
    results = sum(results, [])
        #p = [ (r[1], r[3]) for r in results if not r[0]]
    p = [ r[3] for r in results if not r[0] and r[4]==1]
    if p:
        #x, y = zip(*p)
        plt.figure()
        plt.title(f'n_cases {np.mean(p)}')
        #plt.scatter(x, y)
        plt.hist(p, bins=np.arange(50))
    p = [ r[1] for r in results if not r[0] and r[4]==1]
    if p:
        #x, y = zip(*p)
        plt.figure()
        plt.title(f'n_days {np.mean(p)}')
        #plt.scatter(x, y)
        plt.hist(p, bins=np.arange(50))
    plt.show()

    print(results)
    
    
simulate()
    
    	             
