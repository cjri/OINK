

rm -r output/
mkdir output/

for pd in 0.1 0.08 0.06 0.04; do
    p=output/${pd}
    mkdir ${p}
    ./oink --replicas 100000 --output_prefix ${p}/ --detect ${pd} &
done
wait
