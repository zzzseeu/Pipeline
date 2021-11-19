import itertools
import logging
import os
import subprocess
import sys
import time
from collections import defaultdict
from datetime import timedelta
from functools import wraps
import glob

import pandas as pd
import pysam


def log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10)  # debug

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper


@log
def parse_mapfile(mapfile):
    """
    Feature
    - Parse mapfile to dict.
    Input
    - Mapfile. Four columns, 1st col is sample name, 2nd col is Fastq data path.
    Output
    - Samples_dict. Key: sample, value: {'path': path/to/fastq}.
    """
    df = pd.read_csv(mapfile, sep='\t', header=None, index_col=0)
    samples = df.index.tolist()
    samples_dict = defaultdict(dict)
    for sample in samples:
        path = df.loc[sample, 1]
        samples_dict[sample]['path'] = path
    return samples_dict


def get_read(library_id, library_path, read='1'):
    read1_list = [f'_{read}', f'_R{read}', f'_R{read}_001']
    fq_list = ['fq', 'fastq']
    suffix_list = ["", ".gz"]
    read_pattern_list = [
        f'{library_path}/*{library_id}{read}.{fq_str}{suffix}' 
        for read in read1_list 
        for fq_str in fq_list 
        for suffix in suffix_list
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = sorted(non_empty for non_empty in fq_list if non_empty)
    fq_list = list(itertools.chain(*fq_list))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            '\n'
            f'Invalid Read{read} path! \n'
            f'library_id: {library_id}\n'
            f'library_path: {library_path}\n'
        )
    return fq_list


def get_fq(samples_dict, paired_end):
    """
    Feature
    - Get fastq path.
    Input
    - Samples_dict. Key: sample, value: {'path': path/to/fastq}.
    Output
    - Samples_dict. Key: sample, value: {'path': path/to/fastq,
                                        'R1': full_path/to/R1 read, 
                                        'R2': full_path/to/R2 read}.
    """
    # define read pattern
    for sample in samples_dict:
        path = samples_dict[sample]['path']
        fq1_list = get_read(sample, path)
        samples_dict[sample]['R1'] = fq1_list[0]
        if paired_end:
            fq2_list = get_read(sample, path, '2')
            samples_dict[sample]['R2'] = fq2_list[0]

    return samples_dict


def check_dir(dir):
    """
    Check directory path. If directory is not exists, will automatically make directory.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)


@log
def process_sam(sam):
    prefix = sam.rstrip('.sam')
    cmd = (
        f'samtools view -S {sam} -b > {prefix}.bam; '
        f'samtools sort -O bam -o {prefix}.sorted.bam {prefix}.bam; '
        f'samtools index {prefix}.sorted.bam; '
        f'rm {sam} '
    )
    process_sam.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
