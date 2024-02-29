
export JULIA_NUM_THREADS=4
rm -r output_julia/
mkdir output_julia/

while IFS= read -r t; do
#for t in 1 2 28 90; do
    time julia Oink_timepoint_orig.jl -n 100000 -o output_julia/${t} -e ${t}
done < Data/Time_points.dat

time julia Oink.jl -n 100000
wait
