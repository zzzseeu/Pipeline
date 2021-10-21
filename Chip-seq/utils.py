import pysam
import pandas as pd
import subprocess
import sys
import os
import logging
import time
from datetime import timedelta
from functools import wraps
from collections import defaultdict
from glob import glob


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
    - Mapfile. Four columns, 1st col is sample name, 2nd col is Fastq data path, 3rd col is species.

    Output
    - Samples_dict. Key: sample, value: {'path': path/to/fastq, 'species': species}.

    """
    df = pd.read_csv(mapfile, sep='\t', header=None, index_col=0)
    samples = df.index.tolist()
    samples_dict = defaultdict(dict)
    for sample in samples:
        path = df.loc[sample, 1]
        species = df.loc[sample, 2]
        samples_dict[sample]['path'] = path
        samples_dict[sample]['species'] = species
    return samples_dict


def get_fq(samples_dict, paired_end):
    """
    Feature
    - Get fastq path.

    Input
    - Samples_dict. Key: sample, value: {'path': path/to/fastq, 'species': species}.

    Output
    - Samples_dict. Key: sample, value: {'path': path/to/fastq, 'species': species, 
                                        'R1': full_path/to/R1 read, 'R2': full_path/to/R2 read}.

    """
    for sample in samples_dict:
        path = samples_dict[sample]['path']
        if glob(f'{path}/{sample}*R1.fastq.gz'):
            R1_read = glob(f'{path}/{sample}*R1.fastq.gz')[0]
            samples_dict[sample]['R1'] = R1_read
        else:
            raise Exception('Invalid path to R1 read!')
        if paired_end:
            if glob(f'{path}/{sample}*R2.fastq.gz'):
                R2_read = glob(f'{path}/{sample}*R2.fastq.gz')[0]
                samples_dict[sample]['R2'] = R2_read
            else:
                raise Exception('Invalid path to R2 read!')

    return samples_dict


def check_dir(dir):
    """Check
    """
    if not os.path.exists(dir):
        os.makedirs(dir)


@log
def process_bam(sam):
    prefix = sam.rstrip('.sam')
    cmd = (
        f'samtools view -S {sam} -b > {prefix}.bam; '
        f'samtools sort -O bam -o {prefix}.sorted.bam {prefix}.bam; '
        f'samtools index {prefix}.sorted.bam; '
        f'rm {sam} '
    )
    process_bam.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
