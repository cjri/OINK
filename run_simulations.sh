

rm -r output/
mkdir output/

for pd in 0.1 0.08 0.06 0.04; do
    p=output/${pd}
    mkdir ${p}
    ./oink --replicas 1000000 --output_prefix ${p}/ --detect ${pd} --first_detect 0.0 &
done
wait

rm -r output2/
mkdir output2/

for pd in 0.04; do
    p=output2/${pd}
    mkdir ${p}
    ./oink --replicas 1000000 --output_prefix ${p}/ --detect ${pd} --first_detect 0.0 --max_R0 2 &
done
wait
