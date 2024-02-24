
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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

print(list(data))

## Figure 1a


cols = ['r', 'y', 'g', 'b'][::-1]

## Probability that the outbreak is over (no currently infected individuals at observation time)

def calc_first_day(days, probs):
    try:
        return next(d for d, p in zip(days, probs) if p>=0.95)
    except StopIteration:
        return None


def calc_median(days, probs):
    ## Approx
    try:
        return next(d for d, p in zip(days, np.cumsum(probs)) if p>=0.5)
    except StopIteration:
        return None


    
## Figure 1c
## outbreak origin times - "Origin_stats_*.dat OutputOutbreakTimeStatistics
## Historical nowcasting of outbreak origin times (14 days after first detection)

plt.figure()
for k, c in zip(data, cols):
    o = data[k]['sizes']


    sel_keys = sorted([v for v in o.keys() if 10<=v<=30])

    medians = []

    for v in sel_keys:    
        oo = o[v]
        sizes, probs = zip(*oo)
        sizes = np.array(sizes)
        probs = np.array(probs)
        print(sizes[0], probs[0])
        medians.append(calc_median(sizes, probs))

    plt.plot(sel_keys, medians, label=k, c=c)
plt.legend()
plt.show()

