

NS=1000000

rm -r output_variable_infection_length/
mkdir output_variable_infection_length/

for pd in 0.1 0.08 0.06 0.04; do
    p=output_variable_infection_length/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd}  --first_detect 0.0 --infection_length_min 5 --infection_length_max 9 &
done
wait

