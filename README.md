# OINK

Code supporting Epidemiological inference at the threshold of data availability: an influenza A(H1N2)v spillover event in the United Kingdom 

John A. Fozard Emma C. Thomson Christopher J. R. Illingworth 

https://www.biorxiv.org/content/10.1101/2024.03.11.584378v1

## Dataset

A single detected case: see Data/Detections.dat. Note that the code has only been tested for a single detected case, and has not been optimized for multiple detected cases.

## Using the code

Requires a c++ compiler (e.g. gcc9 onwards), the GSL library, and Mathematica for plotting outputs.

On Debian/Ubuntu (e.g. on a GCP E2 instance)
```
sudo apt-get install git
git clone https://github.com/cjri/OINK
sudo apt-get install g++ libgsl-dev make
cd OINK
make PARALLEL=1
bash scripts/run_simulations_slurm.sh
bash scripts/run_simulations_extra_1.sh
bash scripts/run_simulations_extra_2.sh
bash scripts/run_simulations_extra_3.sh
```
(Note that OMP_NUM_THREADS in these scripts may need to be changed.)

Plots shown in the manuscript are generated using the Mathematica notebooks in `Analysis/`.
Without Mathematica, similar plots may be generated using the code found in `python_scripts/'

Python environment requires Python 3 (tested on 3.10.9) `numpy`, `scipy`, `matplotlib` and `weightedstats`; install these via pip or conda.

For Figure 1, after running the simulations (`run_simulations_slurm.sh`)
```
python python_scripts/plot_data_single.py output/
```
For the supplementary figures (with slightly different formatting)

```
python python_scripts/plot_S1.py
python python_scripts/plot_S2.py
python python_scripts/plot_S3.py
```
for the remaining supplementary figures
```
python python_scripts/plot_data_single.py output_variable_detection_delay/
python python_scripts/plot_data_single.py output_enhanced_detection/
python python_scripts/plot_data_single.py output_variable_infection_length/
```

Full instructions can be found in `Instructions.txt`
Further description of the simulation algorithm can be found in `Algorithm.txt`
