
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

    file_path = path / "P_Outbreak_End.dat"
    outbreak_end = read_data(file_path)

    return outbreak_end

        

if len(sys.argv)>2:
    path = Path(sys.argv[1])
    path2 = Path(sys.argv[2])
else:
    path = Path('output/')
    path2 = Path('output2/')
  

end4 = load_dir_data(path / '0.04')
end2 = load_dir_data(path2/ '0.04')

## Figure 1a


cols = ['r', 'y', 'g', 'b'][::-1]


## Figure 1d
## "Acceptance_rate{i}.dat" - OutputAcceptanceRates
## The acceptance rate P(accepted at t0 | R0), taken at the last timepoint.
## This is then normalized to give a retrospective estimate of the posterior distribution of R0


plt.figure()
    

days, probs = zip(*end4)

plt.plot(days, probs, label='4')
plt.xlim(0, 30)


days, probs = zip(*end2)

plt.plot(days, probs, label='2')
plt.xlim(0, 30)

plt.show()
