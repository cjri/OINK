
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import weightedstats

def read_data(file_path):
    data = []
    with open(file_path, 'r') as f:
        try:
            while True:
                s = f.readline().strip()
                if not s :
                    break
                s = s.split(' ')
                data.append((int(s[0]), float(s[1])))
        except IOError:
            pass
    return data

def load_dir_data(path):
    accept = {}
    origin = {}
    sizes = {}
    for i in range(1, 41):
        file_path = path / f"Acceptance_rate{i}.dat"
        accept[i] = read_data(file_path)

    for s in path.glob("Origin_stats_*.dat"):
        i = int(s.name[len("Origin_stats_"):-4])
        origin[i] = read_data(s)

    for s in path.glob("Sizes_time_*.dat"):
        i = int(s.name[len("Sizes_time_"):-4])
        sizes[i] = read_data(s)

    file_path = path / "P_Outbreak_End.dat"
    outbreak_end = read_data(file_path)

    return {'accept':accept, 'origin':origin, 'sizes':sizes, 'outbreak_end':outbreak_end}

def load_all_data(path):
    data = {}
    for subdir in sorted(path.glob('*/')):
        data[subdir.name] = load_dir_data(subdir)
    return data
        
data = load_all_data(Path('output/'))


## Figure 1a


cols = ['r', 'y', 'g', 'b'][::-1]

## Probability that the outbreak is over (no currently infected individuals at observation time)

def calc_first_day(days, probs):
    try:
        return next(d for d, p in zip(days, probs) if p>=0.95)
    except StopIteration:
        return None

# P_Outbreak_End.dat -  OutputProbabilityEnded
# Historical now-casting that there are no infected individuals
plt.figure()
for k, c in zip(data, cols):
    
    oe = data[k]['outbreak_end']
    days, probs = zip(*oe)
    day_first = calc_first_day(days, probs)
    print(k, day_first)
    
    plt.plot(days, probs, label=k, c=c)
    plt.xlim(0, 30)
plt.legend()

## Figure 1b

# Origin_stats_*.dat" - OutputOutbreakTimeStatistics
#  Retrospective (nowcasting well after the event) of outbreak origin times 
plt.figure()
for k, c in zip(data, cols):
    o = data[k]['origin']
    o = o[max(o)]
    days, probs = zip(*o)
    days = np.array(days)
    probs = np.array(probs)
    i = np.argmax(probs)
    m = np.sum(days*probs)
    print(k, 'peak: ', i, 'mean :', m)
    plt.plot(days, probs, label=k, c=c)
    plt.xlim(50, 10)
plt.legend()

def calc_median(days, probs):
    ## Approx
    try:
        return next(d for d, p in zip(days, np.cumsum(probs)) if p>=0.5)
    except StopIteration:
        return None

def other_median(days, probs):
    days = np.array(days, dtype=float)
    probs = np.array(probs, dtype=float)
    u = np.cumsum(probs)
    print('cs', u)
    q = np.searchsorted(u, 0.5 * u[-1])
    return np.where(u[q]/u[-1] == 0.5, 0.5 * (days[q] + days[q+1]), days[q])

"""
test_sizes = [2, 3, 4]
test_probs = [0.25, 0.25, 0.5]

print('median', calc_median(test_sizes, test_probs))
"""

## Figure 1c
## outbreak origin times - "Origin_stats_*.dat OutputOutbreakTimeStatistics
## Historical nowcasting of outbreak origin times (14 days after first detection)

plt.figure()
for k, c in zip(data, cols):
    o = data[k]['sizes']
    o = o[14]
#    print(o)
    sizes, probs = zip(*o)

    print(k, calc_median(sizes, probs))
    print(k, other_median(sizes, probs))
    
    plt.plot(sizes, probs, label=k, c=c)
    plt.xlim(0, 150)
plt.legend()


## Figure 1d
## "Acceptance_rate{i}.dat" - OutputAcceptanceRates
## The acceptance rate P(accepted at t0 | R0), taken at the last timepoint.
## This is then normalized to give a retrospective estimate of the posterior distribution of R0


plt.figure()
for k,c in zip(data, cols):
    a = data[k]['accept']
    o = [ (i, b[-1][1]) for i, b in a.items() ]
    r0val, probs = zip(*o)
    probs = np.array(probs)
    probs/=probs.sum()
    idx = np.argmax(probs)
    print(k, r0val[idx]*0.1, np.sum(probs[0:22]), np.sum(probs[0:21]))
    plt.plot(np.array(r0val)*0.1, probs, label=k, c=c)
plt.legend()
          
    
plt.show()

