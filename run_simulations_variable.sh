

rm -r output_variable_symptom_onset/
mkdir output_variable_symptom_onset/

NS=10000

## Variation between symptom onset and detection

for pd in 0.1 0.08 0.06 0.04; do
    p=output_variable_symptom_onset/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd}  --first_detect 0.0 --detect_min 14 --detect_max 22 &
done
wait

## Variation in time infected

rm -r output_variable_infected/
mkdir output_variable_infected/

for pd in 0.1 0.08 0.06 0.04; do
    p=output_variable_infected/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd}  --first_detect 0.0 --infection_length_min 5 --infection_length_max 9 &
done
wait

rm -r output_enhanced_detection/
mkdir output_enhanced_detection/

## Enhanced surveillance after first case

for pd in 0.1 0.08 0.06 0.04; do
    p=output_enhanced_detection/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --detect_enhanced 0.2 --first_detect 0.0 &
done
wait

