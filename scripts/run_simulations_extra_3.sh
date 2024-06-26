

rm -r output_enhanced_detection/
mkdir output_enhanced_detection/

NS=1000000

## Enhanced surveillance after first case

for pd in 0.1 0.08 0.06 0.04; do
    p=output_enhanced_detection/${pd}
    mkdir ${p}
    ./oink --replicas ${NS} --output_prefix ${p}/ --detect ${pd} --detect_enhanced 0.2 --first_detect 0.0 &
done
wait

