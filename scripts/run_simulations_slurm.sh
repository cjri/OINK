

rm -r output/
mkdir output/

NS=1000000

export OMP_NUM_THREADS=16

for pd in 0.1 0.08 0.06 0.04; do
    p=output/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd}  --first_detect 0.0 &
done
wait

rm -r output2/
mkdir output2/

export OMP_NUM_THREADS=64

pd=0.04
p=output2/${pd}
mkdir ${p}
./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --first_detect 0.0 --max_R0 2 
