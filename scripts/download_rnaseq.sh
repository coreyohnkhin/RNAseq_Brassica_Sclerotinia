while read run; do
    software/sratoolkit.3.2.1-ubuntu64/bin/prefetch $run -O data
    software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump -S -O data/ data/$run
    rm -rf data/$run
done < $1
