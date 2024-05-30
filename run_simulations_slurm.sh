

rm -r output/
mkdir output/

NS=1000000

export OMP_NUM_THREADS=16

for pd in 0.1 0.08 0.06 0.04; do
    p=output/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --detect_after_first ${pd} --first_detect 0.0 &
done
wait
