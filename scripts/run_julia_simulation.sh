
export JULIA_NUM_THREADS=4
rm -r output_julia/
mkdir output_julia/

while IFS= read -r t; do
#for t in 1 2 28 90; do
    #echo ${t}
    time julia julia_alternative_simulations/Oink_timepoint.jl -n 100000 -o output_julia/${t} -e ${t}
done < Data/Time_points.dat

time julia julia_alternative_simulations/Oink.jl -n 100000 -o output_julia/
wait
