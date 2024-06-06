
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
                data.append((float(s[0]), float(s[1])))
        except IOError:
            pass
    return data


path = "output_julia/"
accept = {}
for fn in Path(path).glob('*_R0.dat'):
    i = int(fn.name[:-7])
    accept[i] = read_data(fn)

accept = dict((k,accept[k]) for k in sorted(accept))
    

plt.figure()
for i in accept:
    r, a = zip(*accept[i])
    a = np.array(a)
    a = a/np.sum(a)
    plt.plot(r, a, label=str(i))
plt.legend()
plt.show()


accept = {}
for fn in Path(path).glob('*_time.dat'):
    i = int(fn.name[:-9])
    accept[i] = read_data(fn)

accept = dict((k,accept[k]) for k in sorted(accept))
    

plt.figure()
for i in accept:
    r, a = zip(*accept[i])
    a = np.array(a)
    a = a/np.sum(a)
    plt.plot(r, a, label=str(i))
plt.legend()
plt.show()
