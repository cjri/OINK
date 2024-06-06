
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys

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


if len(sys.argv)>1:
    path = sys.argv[1]
else:
    path = 'output/'
  
data = load_all_data(Path(path))

print(list(data))

## Figure 1a


cols = ['r', 'y', 'g', 'b'][::-1]

## Probability that the outbreak is over (no currently infected individuals at observation time)

def calc_first_day(days, probs):
    try:
        return next(d for d, p in zip(days, probs) if p>=0.95)
    except StopIteration:
        return None


## Figure S2a


## Figure 1d
## "Acceptance_rate{i}.dat" - OutputAcceptanceRates
## The acceptance rate P(accepted at t0 | R0), taken at the last timepoint.
## This is then normalized to give a retrospective estimate of the posterior distribution of R0


plt.figure()
d = data['0.04']
    
a = d['accept']
print(list(a[1]))

for u in range(len(a[1])):
    
    o = [ (i, b[u][1]) for i, b in a.items() ]
    r0val, probs = zip(*o)
    probs = np.array(probs)
    probs/=probs.sum()
    idx = np.argmax(probs)
    plt.plot(np.array(r0val)*0.1, probs)

plt.show()    

    
## Figure S2b

# Origin_stats_*.dat" - OutputOutbreakTimeStatistics
#  Nowcasting of outbreak origin times 
plt.figure()

oo = data['0.04']['origin']
print(list(oo))
for u in sorted(oo):
    print(u)
    o = oo[u]
    days, probs = zip(*o)
    days = np.array(days)
    probs = np.array(probs)
    i = np.argmax(probs)
    m = np.sum(days*probs)
    plt.plot(days, probs)
    plt.xlim(50, 10)
#plt.legend()
          
    
plt.show()

