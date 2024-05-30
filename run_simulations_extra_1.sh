

rm -r output_variable_symptom_onset/
mkdir output_variable_symptom_onset/

NS=1000000

## Variation between symptom onset and detection

for pd in 0.1 0.08 0.06 0.04; do
    p=output_variable_symptom_onset/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd}  --first_detect 0.0 --detect_min 14 --detect_max 22 &
done
wait


