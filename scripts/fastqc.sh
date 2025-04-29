#!/bin/bash

mkdir data/fastqc_output
fastqc data/*fastq.gz -o data/fastqc_output
