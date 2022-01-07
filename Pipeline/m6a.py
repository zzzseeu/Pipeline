import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import glob
import pysam

from utils import *


@log
def fastqc(outdir, sample, thread, R1_read, R2_read):
    fqdir = f'{outdir}/01.fastqc/{sample}'
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
def bowtie2(outdir, sample, thread, species, genomer_dir, R1_read, R2_read):
    bt2dir = f'{outdir}/02.mapping/{sample}'
    check_dir(bt2dir)
    if R2_read != 'None':
        string = f'-1 {R1_read} -2 {R2_read} '
    else:
        string = f'-U {R1_read} '
    cmd = (
        f'bowtie2 -p {thread} --trim5 8 --local '
        f'-x {genomer_dir}/{species}/bowtie2_index/{species} '
        f'{string} -S {bt2dir}/{sample}.sam >& {bt2dir}/{sample}.alignment.stats '
    )
    bowtie2.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
    process_sam(f'{bt2dir}/{sample}.sam')


@log
def get_mapping_statement(sample, stat_file, paired_end=False):
    alignment_stat_file = open(stat_file, 'r')
    if not paired_end:
        for line in alignment_stat_file.readlines():
            if ' reads; of these:' in line:
                attrs = line.strip(' ').split(' ')
                all_reads = int(attrs[0])
            if ' were unpaired; of these:' in line:
                attrs = line.strip(' ').split(' ')
                unpaired_reads = int(attrs[0])
            if ' aligned 0 times' in line:
                attrs = line.strip(' ').split(' ')
                unmapped_reads = int(attrs[0])
            if ' aligned exactly 1 time' in line:
                attrs = line.strip(' ').split(' ')
                unique_mapped_reads = int(attrs[0])
            if ' aligned >1 times' in line:
                attrs = line.strip(' ').split(' ')
                multi_mapped_reads = int(attrs[0])
            if ' overall alignment rate' in line:
                attrs = line.strip(' ').split(' ')
                alignment_rate = attrs[0]
        # df = pd.DataFrame([[all_reads, unpaired_reads, unmapped_reads, unique_mapped_reads, multi_mapped_reads, alignment_rate]], columns=['All_reads', 'Unpaired_reads', 'Unmapped_reads', 'Unique_mapped_reads', 'Multi_mapped_reads', 'Alignment_rate'])
        s = f'{sample}\t{all_reads}\t{unpaired_reads}\t{unmapped_reads}\t{unique_mapped_reads}\t{multi_mapped_reads}\t{alignment_rate}\n'
    else:
        for line in alignment_stat_file.readlines():
            if ' reads; of these:' in line:
                attrs = line.strip(' ').split(' ')
                all_reads = int(attrs[0]) * 2
            if ' were paired; of these:' in line:
                attrs = line.strip(' ').split(' ')
                paired_reads = int(attrs[0]) * 2
            if ' aligned concordantly exactly 1 time' in line:
                attrs = line.strip(' ').split(' ')
                unique_mapped_paired_reads = int(attrs[0]) * 2
            if ' aligned concordantly >1 times' in line:
                attrs = line.strip(' ').split(' ')
                multi_mapped_paired_reads = int(attrs[0]) * 2
            if ' aligned discordantly 1 time' in line:
                attrs = line.strip(' ').split(' ')
                paired_discordantly_mapped_reads = int(attrs[0]) * 2
            if ' aligned exactly 1 time' in line:
                attrs = line.strip(' ').split(' ')
                unique_mapped_mates = int(attrs[0])
            if ' aligned >1 times' in line:
                attrs = line.strip(' ').split(' ')
                multi_mapped_mates = int(attrs[0])
            if ' overall alignment rate' in line:
                attrs = line.strip(' ').split(' ')
                alignment_rate = attrs[0]
        mapped_reads = unique_mapped_paired_reads + multi_mapped_paired_reads + paired_discordantly_mapped_reads + unique_mapped_mates + multi_mapped_mates
        unmapped_reads = all_reads - mapped_reads
        unique_mapped_reads = unique_mapped_paired_reads + unique_mapped_mates + paired_discordantly_mapped_reads
        multi_mapped_reads = multi_mapped_paired_reads + multi_mapped_mates
        # df = pd.DataFrame([[all_reads, paired_reads, unmapped_reads, unique_mapped_reads, multi_mapped_reads, alignment_rate]], columns=['All_reads', 'Paired_reads', 'Unmapped_reads', 'Unique_mapped_reads', 'Multi_mapped_reads', 'Alignment_rate'])
        s = f'{sample}\t{all_reads}\t{paired_reads}\t{unmapped_reads}\t{unique_mapped_reads}\t{multi_mapped_reads}\t{alignment_rate}\n'
    res = open(f'{os.path.dirname(stat_file)}/mapping_statment.txt', 'w')
    res.write(s)


@log
def mapping_summary(outdir, paired_end=False):
    stat_files = glob.glob(f'{outdir}/02.mapping/*/mapping_statment.txt')
    s = ' '.join(stat_files)
    cmd = (
        f'cat {s} > {outdir}/02.mapping/MappingSummary.txt'
    )
    os.system(cmd)
    if paired_end:
        c = ['sample', 'All_reads', 'Paired_reads', 'Unmapped_reads', 'Unique_mapped_reads', 'Multi_mapped_reads', 'Alignment_rate']
    else:
        c = ['sample', 'All_reads', 'Unpaired_reads', 'Unmapped_reads', 'Unique_mapped_reads', 'Multi_mapped_reads', 'Alignment_rate']
    df = pd.read_csv(f'{outdir}/02.mapping/MappingSummary.txt', sep='\t', header=None)
    df.columns = c
    df.to_csv(f'{outdir}/02.mapping/MappingSummary.txt', sep='\t', index=False)


@log
def process(outdir, sample, thread, R1_read, R2_read, genome_dir, species, paired_end):
    fastqc(outdir, sample, thread, R1_read, R2_read)
    bowtie2(outdir, sample, thread, species, genome_dir, R1_read, R2_read)
    stat_file = f'{outdir}/02.mapping/{sample}/{sample}.alignment.stats'
    get_mapping_statement(sample, stat_file, paired_end)


def run_process(args):
    samples_dict = parse_mapfile(args.mapfile)
    samples_dict = get_fq(samples_dict, args.paired_end)
    samples = [s for s in samples_dict]
    outdirs = [args.outdir] * len(samples)
    threads = [args.thread] * len(samples)
    R1_reads = [samples_dict[s]['R1'] for s in samples_dict]
    if args.paired_end:
        R2_reads = [samples_dict[s]['R2'] for s in samples_dict]
    else:
        R2_reads = ['None'] * len(samples)
    species = [args.species] * len(samples)
    paired_end = [args.paired_end] * len(samples)
    genome_dirs = [args.genome_dir] * len(samples)
    if len(samples) > 5:
        p = 5
    else:
        p = len(samples)
    process_pools = []
    with ProcessPoolExecutor(p) as pool:
        for i in pool.map(process, outdirs, samples, threads, R1_reads, R2_reads, genome_dirs, species, paired_end):
            process_pools.append(i)
    multiqc(f'{args.outdir}/01.fastqc', f'{args.outdir}/01.fastqc')
    mapping_summary(args.outdir, args.paired_end)
    os.system(f'rm {args.outdir}/02.mapping/*/mapping_statment.txt')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', 
                        help='Mapfile, which includes 3 columns. 1st column is sample name, 2nd column is fastq file path', 
                        required=True)
    parser.add_argument('--species', 
                        help='Species', required=True)
    parser.add_argument('--genome_dir', help='Genome index directory', required=True)
    parser.add_argument('--outdir', 
                        help='Output directory, which should includes fastqc, bowtie2, multiqc and bamCoverage results.', 
                        default='./')
    parser.add_argument('--paired_end', 
                        help='Paired end data or not.', 
                        action='store_true')
    parser.add_argument('--thread', 
                        help='Threads for each step to run.', 
                        default=3)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    run_process(args)


if __name__ == '__main__':
    main()
