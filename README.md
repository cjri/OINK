# OINK

Code supporting Epidemiological inference at the threshold of data availability: an influenza A(H1N2)v spillover event in the United Kingdom 

John A. Fozard Emma C. Thomson Christopher J. R. Illingworth 

https://www.biorxiv.org/content/10.1101/2024.03.11.584378v1

## Dataset

A single detected case: see Data/Detections.dat. Note that the code has only been tested for a single detected case, and has not been optimized for multiple detected cases.

## Using the code

Requires a c++ compiler (e.g. gcc9 onwards), the GSL library, and Mathematica for plotting outputs.

On ubuntu
```
sudo apt-get install g++ libgsl-dev
make PARALLEL=1
bash scripts/run_simulations_slurm.sh
bash scripts/run_simulations_extra_1.sh
bash scripts/run_simulations_extra_2.sh
bash scripts/run_simulations_extra_3.sh
```
(Note that OMP_NUM_THREADS in these scripts may need to be changed.)

Full instructions can be found in `Instructions.txt`



