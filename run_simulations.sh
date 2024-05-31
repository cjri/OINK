

rm -r output/
mkdir output/

NS=1000

for pd in 0.1 0.08 0.06 0.04; do
    p=output/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --detect_after_first ${pd} --first_detect 0.0 &
done
wait

rm -r output2/
mkdir output2/

pd=0.04
p=output2/${pd}
mkdir ${p}
./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --first_detect 0.0 --detect_after_first ${pd} --max_R0 2 

rm -r output3/
mkdir output3/

pd=0.1
for pd2 in 0.1 0.2; do
    p=output3/${pd}-${pd2}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --first_detect 0.0 --detect_after_first ${pd} --max_R0 2 
done


wait
