from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# Compare Oink.jl, Oink_timepoint.jl and C++ simulation output

# Compare specific timepoints

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

### Figure 1d

sim_probs = {}

for tp in [1, 2, 28, 90]:

    for k in ['0.1']:
        a = data[k]['accept']
        o = []

        for i,b in a.items():
            o.append((i, next(c[1] for c in b if c[0]==tp)))    
        #print(o)

        r0val, probs = zip(*o)
        probs = np.array(probs)
        probs/=probs.sum()
        sim_probs[tp] = probs

sim_times = {}
for tp in [1, 2, 28, 90]:

    for k in ['0.1']:
        o = data[k]['origin'][tp]
          
        #print(o)

        times, probs = zip(*o)
        probs = np.array(probs)
        probs/=probs.sum()
        sim_times[tp] = probs

path = "output_julia/"
accept = {}
for fn in Path(path).glob('*_R0.dat'):
    i = int(fn.name[:-7])
    accept[i] = np.array(read_data(path + f"{i}_R0.dat", type_first=float))[:,1]

accept = dict((k,accept[k]) for k in sorted(accept))
prob_julia = dict((k, accept[k]/np.sum(accept[k])) for k in accept)

path = "output_julia/"
times = {}
for fn in Path(path).glob('*_time.dat'):
    i = int(fn.name[:-9])
    times[i] = np.array(read_data(path + f"{i}_time.dat", type_first=float))[:,1]

times_julia = dict((k,times[k]) for k in sorted(times))
print(list(times_julia))

prob_no_end = np.array(read_data(path + "julia_R0_acceptance_rates.dat", type_first=float))[:,1]
prob_no_end = prob_no_end/np.sum(prob_no_end)

print(prob_no_end)

print(list(prob_julia))

plt.figure()
plt.plot(sim_probs[90])
plt.plot(prob_julia[28])
plt.plot(prob_no_end)

for tp in [1, 2, 28]:
    plt.figure()
    plt.plot(sim_times[tp])
    plt.plot(times_julia[tp])

times_no_end = np.array(read_data(path + "julia_output_time.dat", type_first=float))[:,1]
times_no_end = times_no_end/np.sum(times_no_end)

plt.figure()
plt.plot(sim_times[90])
plt.plot(times_julia[28])
plt.plot(times_no_end)


plt.show()
