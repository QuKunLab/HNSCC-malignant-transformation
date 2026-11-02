#!/usr/bin/env python
# coding: utf-8

import os
kits = {"DZOE2023121932-b1": ["N1","N2","N3","T1","T2","T3"], 
        "DZOE2024011405-b1": ["N1","T1"]}

## fastq quality control
def args_input(kit, sample):
    in_prefix = f"./data/WES/{kit}/Raw_data/{sample}/{sample}."
    out_prefix = f"./fastp/{kit}/{sample}."
    args_base = ["./tools/fastp",
                 "-i", in_prefix+"R1.fastq.gz", "-I", in_prefix+"R2.fastq.gz",
                 "-o", out_prefix+"R1.fastq.gz", "-O", out_prefix+"R2.fastq.gz",
                 "-h", out_prefix+"fastp.html", "-j", out_prefix+"fastp.json"]
    return args_base
args_qc = ["-w","8", "-f","9", "-F","9", "-l","36", "--detect_adapter_for_pe"]

for kit,samples in kits.items():
    for sample in samples:
        print(f"Fastp sample {kit}-{sample}")
        args = args_input(kit, sample)+args_qc
        os.system(" ".join(args))

##bwa  mapping
args_base = ['bwa', 'mem', '-t','18']
def create_args(kit, sample):
    tag = ['-R', f"'@RG\\tID:{kit}-{sample}\\tPL:illumina\\tSM:OSCC'"]
    ref = "./fa/hg38clean.fa"
    prefix = f"./fastp/{kit}/{sample}."
    out = f"./bwa/{kit}-{sample}.mappings.sorted.bam"
    args = tag + [ref, 
                  prefix+"R1.fastq.gz", prefix+"R2.fastq.gz",
                  "|", "samtools view -@ 18 -Sb -",
                  "|", "samtools sort -@ 18 -O BAM -o", out, "-"]
    return args

for kit,samples in kits.items():
    for sample in samples:
        print(f"BWA sample {kit}-{sample}")
        args = args_base + create_args(kit, sample)
        print(" ".join(args))
        os.system(" ".join(args))

