process_file() {
    local input_file="$1"
    local base_name=$(basename "$input_file" | cut -d '.' -f 1)
    local output_dir="data/$base_name"

    mkdir -p "$output_dir"

    echo "Downloading SRA runs from $input_file to directory: $output_dir"

    while read run || [ -n "$run" ]; do
        if [ -z "$run" ]; then
            continue
        fi
        echo "Processing run: $run"
        software/sratoolkit.3.2.1-ubuntu64/bin/prefetch "$run" -O "$output_dir"
        software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump -S -O "$output_dir" "$output_dir/$run"
        rm -rf "$output_dir/$run"
    done < "$input_file"

    echo "All downloads from $input_file completed in: $output_dir"
}

if [ $# -eq 0 ]; then
    echo "Error: No input files provided"
    echo "Usage: $0 input_file1.txt [input_file2.txt] [input_file3.txt] ..."
    exit 1
fi

for input_file in "$@"; do
    if [ ! -f "$input_file" ]; then
        echo "Error: File $input_file does not exist"
        continue
    fi
    process_file "$input_file"
done
