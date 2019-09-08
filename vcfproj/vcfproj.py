from future.utils import lmap
from past.builtins import xrange
import pandas as pd

from collections import namedtuple
import re

import gzip
import bz2
import urllib
import numpy as np
from ncls import NCLS

"""
The VCF reaser is structured as KB's g2gtools
can be accessed here
https://github.com/churchill-lab/g2gtools/tree/master/g2gtools

Input: GTF, VCF
Output: VCF-like file
"""

VCF_FIELDS = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'samples']
VCFRecord = namedtuple('VCFRecord', VCF_FIELDS)

GT_DATA_FIELDS = ['ref', 'left', 'right', 'gt', 'fi', 'phase', 'gt_left', 'gt_right', 'is_snp']
GTData = namedtuple('GTData', GT_DATA_FIELDS)

GENOTYPE_UNPHASED = '/'
GENOTYPE_PHASED = '|'

REGEX_ALT = re.compile("(^[A|C|G|T]+)")

def parse_vcf_line(line):
    """
    Parse a line in the VCF file.
    :param line: a line from the VCF file
    :type line: str
    :return: :class:`.vcf.VCFRecord`
    """

    if isinstance(line, str):
        if line.startswith('#'):
            return None

        elem = line.strip().split('\t')
    elif isinstance(line, list):
        elem = line

    try:
        quality = int(elem[5])
    except ValueError:
        try:
            quality = float(elem[5])
        except ValueError:
            quality = None

    filter_field = None

    if elem[6] != '.':
        filter_field = elem[6].split(';')

    info = elem[7]

    try:
        fmt = elem[8]
    except IndexError:
        fmt = None
    else:
        if fmt == '.':
            fmt = None

    return [elem[0], int(elem[1]), None if elem[2] == '.' else elem[2], elem[3],
                     elem[4].split(','), quality, filter_field, info, fmt, elem[9:]]

def open_resource(resource, mode='rb'):
    """
    Open different types of files and return the handle.
    :param resource: a file located locally or on the internet.  
     Gzip'd and zip'd files are handled.
    :param mode: standard file open modes
    :return: the resource (file) handle
    """
    if not resource:
        return None

    resource = str(resource)

    if resource.endswith(('.gz', '.Z', '.z')):
        return gzip.open(resource, mode)
    elif resource.endswith(('.bz', '.bz2', '.bzip2')):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(('http://', 'https://', 'ftp://')):
        return urllib.urlopen(resource)
    else:
        return open(resource, 'r')


def tiny_vcf_reader(resource):
    reader = open_resource(resource, 'r')

    # skip header
    df = pd.DataFrame(columns=VCF_FIELDS)
    current_line = next(reader)x

    while current_line.startswith('##'):
        current_line = next(reader)

    if current_line.startswith('#'):
        elems = current_line.strip().split('\t')
        samples = elems[9:]
        samples = dict(zip(samples, (x for x in xrange(len(samples)))))
    else:
        print('Improperly formatted VCF file')

    records = []
    count = 0
    for line in reader:
        record = parse_vcf_line(line)
        print(count, end = '\r')
        records += [record]

    return records

def projection(GTF_FILE, VCF_FILE):
    """
    Projects VCF file to transcript coordinates.

    Parameters
    ----------
    GTF_FILE : string containing GTF file name, assumed to be unzipped
    VCF_FILE : string containing VCF file name, can be a gzipped file 

    Returns
    -------
    vcf_txome : pandas dataframe with the following header
    [
        'chrom_x',
        'gene',
        'txome',
        'relative_pos',
        'transcript_length',
        'id',
        'ref',
        'alt',
        'qual',
        'filter',
        'info',
        'format',
        'samples'
    ]

    call it truncated VCF 
    """
    from timeit import default_timer as timer


    start = timer()

    records = tiny_vcf_reader(VCF_FILE)
    chromosomes = set(gtf_df.chrom.values)
    dataframes = []

    for i in chromosomes:
        gtf_df_subset = gtf_df.loc[gtf_df.chrom == i]
        df_subset = df.loc[df.chrom == i]

        start_val = gtf_df.loc[gtf_df.chrom == i].start.values
        end_val = gtf_df.loc[gtf_df.chrom == i].end.values
        indices = gtf_df.loc[gtf_df.chrom == i].index.values

        query_start_val = df.loc[df.chrom == i].pos.values
        query_end_val = df.loc[df.chrom == i].pos.values + 1
        query_indices = df.loc[df.chrom == i].index.values

        ncls = NCLS(np.array(start_val), np.array(end_val), indices)
        result = ncls.all_overlaps_both(query_start_val, query_end_val, query_indices)
        map_df = pd.DataFrame(list(zip(*result)), columns= ['vcf_index' , 'gtf_index'])
        vcf_gtf_subset = pd.merge(
            gtf_df_subset.join(map_df.set_index('gtf_index')),
            df_subset,
            left_on = 'vcf_index',
            left_index = False,
            right_index = True,
        )

        dataframes += [vcf_gtf_subset]
        print('chromosome ',i,' done')


    print('Merging the cromosome...')

    vcf_gtf = pd.concat(dataframes)
    vcf_gtf['relative_pos'] = vcf_gtf['pos'] - vcf_gtf['start']
    vcf_gtf['transcript_length'] = vcf_gtf['end'] - vcf_gtf['start'] + 1


    end = timer()

    print('elapsed ',(end - start)," seconds")

    return(vcf_gtf[
        [
            'chrom_x',
            'gene',
            'txome',
            'relative_pos',
            'transcript_length',
            'id',
            'ref',
            'alt',
            'qual',
            'filter',
            'info',
            'format',
            'samples'
        ]
    ])
