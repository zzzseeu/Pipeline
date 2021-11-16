import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import pysam

from utils import *

GENOME = '/home/xinzhou/GENOME/database'


@log
def fastqc(outdir, thread, R1_read, R2_read):
    fqdir = f'{outdir}/01.fastqc'
    check_dir(fqdir)
    if R2_read != 'None':
        string = f'{R2_read}'
    else:
        string = ''
    cmd = (
        f'fastqc -o {fqdir} -f fastq '
        f'-t {thread} '
        f'{R1_read} {string} '
    )
    fastqc.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@log
def multiqc(inputdir, outdir):
    cmd = (
        f'multiqc --title Multiqc_report '
        f'-d -dd 1 --module fastqc '
        f'--filename multiqc_report.html '
        f'-o {outdir} '
        f'{inputdir} '
    )
    multiqc.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@log
def bowtie2(outdir, sample, thread, R1_read, R2_read, species):
    bt2dir = f'{outdir}/02.mapping/{sample}'
    check_dir(bt2dir)
    if R2_read != 'None':
        string = f'-1 {R1_read} -2 {R2_read} '
    else:
        string = f'-U {R1_read} '
    cmd = (
        f'bowtie2 -p {thread} --trim5 8 --local '
        f'-x {GENOME}/{species}/bowtie2_index/{species} '
        f'{string}| samtools view -bS - | samtools sort  -O bam -o {bt2dir}/{sample}.bam; '
        f'samtools index {bt2dir}/{sample}.bam; '
        f'samtools flagstat {bt2dir}/{sample}.bam > {bt2dir}/{sample}.alignment.stats '
    )
    bowtie2.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
    # process_bam(f'{bt2dir}/{sample}.sam')


@log
def get_mapped_reads(stat_file):
    alignment_stat_file = open(stat_file, 'r')
    for line in alignment_stat_file.readlines():
        if 'in total (QC-passed reads + QC-failed reads)' in line:
            all_reads = int(re.findall(r'\d+\s', line)[0].strip())
        if 'mapped (' in line:
            mapped_reads = int(re.findall(r'\d+\s', line)[0].strip())
    return all_reads, mapped_reads


@log
def calculate_mapping_rate(outdir):
    stats = glob(f'{outdir}/02.mapping/*/*.alignment.stats')
    alignment_stat = open(f'{outdir}/02.mapping/alignment_summary.txt', 'w')
    for i in stats:
        sample = i.split('/')[-1].split('.')[0]
        all_reads, mapped_reads = get_mapped_reads(i)
        mapped_rates = round(mapped_reads/all_reads, 4)
        mapped_rates = float(mapped_rates * 100)
        alignment_stat.write(f'{sample}\t{all_reads}\t{mapped_reads}\t{mapped_rates}%\n')
    alignment_stat.close()


@log
def normalize(outdir, sample, bam, extend_size, species):
    prefix = f'{outdir}/02.mapping/{sample}'
    check_dir(prefix)
    stat_file = f'{outdir}/02.mapping/{sample}/{sample}.alignment.stats'
    all_reads, mapped_reads = get_mapped_reads(stat_file)
    scale = 1000000 / mapped_reads
    scale = format(scale, '.3f')
    cmd = (
        f'macs2 pileup -f BAM --extsize {extend_size} -i {bam} -o {prefix}/{sample}.bdg; '
        f'macs2 bdgopt -i {prefix}/{sample}.bdg -m multiply -p {scale} -o {prefix}/{sample}_temp_normalized.bdg; '
        f'sed -n \'2,$p\' {prefix}/{sample}_temp_normalized.bdg > {prefix}/{sample}.bdg; '
        f'rm {prefix}/{sample}_temp_normalized.bdg; '
        f'bedSort {prefix}/{sample}.bdg {prefix}/{sample}.bdg; '
        f'bedGraphToBigWig {prefix}/{sample}.bdg {GENOME}/{species}/reference/genome.size {prefix}/{sample}.bw '
    )
    normalize.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@log
def process(outdir, sample, thread, R1_read, R2_read, 
            species, extend_size):
    fastqc(outdir, thread, R1_read, R2_read)
    bowtie2(outdir, sample, thread, R1_read, R2_read, species)
    bam = f'{outdir}/02.mapping/{sample}/{sample}.bam'
    normalize(outdir, sample, bam, extend_size, species)


def run_process(args):
    temp_dict = parse_mapfile(args.mapfile)
    samples_dict = get_fq(temp_dict, args.paired_end)
    samples = list(samples_dict.keys())
    R1_read = [samples_dict[i]['R1'] for i in samples]
    if args.paired_end:
        R2_read = [samples_dict[i]['R2'] for i in samples]
    else:
        R2_read = ['None'] * len(samples)
    species = [args.species] * len(samples)
    outdirs = [args.outdir] * len(samples)
    threads = [args.thread] * len(samples)
    extend_size = [args.extend_size] * len(samples)

    # define cores
    if len(samples) < 6:
        cores = len(samples)
    else:
        cores = 5
    # run
    run_pools = []
    with ProcessPoolExecutor(cores) as pool:
        for i in pool.map(process, outdirs, samples, threads, R1_read, R2_read, species, 
                        extend_size):
            run_pools.append(i)
    multiqc(args.outdir, f'{args.outdir}/01.fastqc/multiqc')
    calculate_mapping_rate(args.outdir)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', 
                        help='Mapfile, which includes 3 columns. 1st column is sample name, 2nd column is fastq file path', 
                        required=True)
    parser.add_argument('--species', 
                        help='Species', required=True)
    parser.add_argument('--outdir', 
                        help='Output directory, which should includes fastqc, bowtie2, multiqc and bamCoverage results.', 
                        default='./')
    parser.add_argument('--paired_end', 
                        help='Paired end data or not.', 
                        action='store_true')
    parser.add_argument('--thread', 
                        help='Threads for each step to run.', 
                        default=3)
    parser.add_argument('--extend_size', 
                        help='The extension size in bps. Each alignment read will become a EXTSIZE of fragment, then be piled up. Check description for -B for detail.', 
                        default=200)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    run_process(args)


if __name__ == '__main__':
    main()


