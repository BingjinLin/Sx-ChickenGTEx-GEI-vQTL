#!/bin/bash

script="/soft/zhh/vQTL/Coloc/New_bash/All_coloc.vQTL.tGWAS.sh"
max_jobs=8
running=0

while IFS= read -r cmd; do

    [[ -z "$cmd" ]] && continue
    [[ "$cmd" =~ ^[[:space:]]*# ]] && continue

    echo "Submitting: $cmd"

    eval "$cmd" &

    ((running++))

    echo "Currently submitted: $running"

    if [ $? -ne 0 ]; then
        echo "WARNING: $cmd failed!"
    else
        echo "SUCCESS: $cmd completed."
    fi

    if [ "$running" -ge "$max_jobs" ]; then

        wait -n

        ((running--))
        echo "One job finished."
        echo "Currently running: $running"
    fi
done < "$script"
wait

echo "================================"
echo "All scripts completed!"
echo "================================"
