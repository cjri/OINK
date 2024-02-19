
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
    for subdir in path.glob('*/'):
        data[subdir.name] = load_dir_data(subdir)
    return data
        
data = load_all_data(Path('output/'))

print(list(data))

## Figure 1a

plt.figure()
for k in data:
    oe = data[k]['outbreak_end']
    print(oe)
    days, probs = zip(*oe)
    plt.plot(days, probs, label=k)
    plt.xlim(0, 30)
plt.legend()

## Figure 1b

plt.figure()
for k in data:
    o = data[k]['origin']
    print(max(o))

    o = o[max(o)]
    days, probs = zip(*o)
    plt.plot(days, probs, label=k)
    plt.xlim(50, 10)
plt.legend()

plt.show()

## Figure 1c
plt.figure()
for k in data:
    o = data[k]['sizes']
    o = o[14]
#    print(o)
    sizes, probs = zip(*o)
    print(np.sum(o[150:]))
    plt.plot(sizes, probs, label=k)
#    plt.xlim(1, 150)
plt.legend()

## Figure 1d
plt.figure()
for k in data:
    a = data[k]['accept']
    o = [ (i, b[-1][1]) for i, b in a.items() ]
    r0val, probs = zip(*o)
    probs = np.array(probs)
    probs/=probs.sum()
    plt.plot(np.array(r0val)*0.1, probs, label=k)
plt.legend()
          
    
plt.show()

