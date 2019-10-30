#!/usr/bin/env python
"""
region_matrix.py

Given a list of N BAM files and a BED file with K NONCONTIGUOUS NONOVERLAPPING
regions, writes the N x K matrix whose (i, j)th element is the AUC of
the i-th sample in the j-th region. Wraps wiggletools, available at
https://github.com/Ensembl/WiggleTools.

Each line of input region BED should be in the following format.
chrom <TAB> start pos (0-based; inclusive) <TAB> end pos (0-based; exclusive)
"""
import multiprocessing
import sys
import subprocess
import os
import signal
import time

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def aucs_from_regions(wiggletools, region_bed, indexed_bam, regions):
    """ Computes AUCs in regions.

        wiggletools: path to wiggletools executable
        region_bed: must have three columns -- chromosome, start coordinate, 
        and end coordinate
        indexed_bam: path to indexed BAM file
        regions: list of regions from region_bed in coordinate-sorted order

        Return value: list of AUCs in same order as region BED
    """
    aucs = [0 for _ in xrange(len(regions))]
    wiggletools_command = (
            '{wiggletools} write_bg - trim {region_bed} {indexed_bam}'.format(
                    wiggletools=wiggletools,
                    region_bed=region_bed,
                    indexed_bam=indexed_bam,
                )
        )
    wiggletools_process = subprocess.Popen(
            wiggletools_command, stdout=subprocess.PIPE, bufsize=-1, shell=True
        )
    # Must allow advancement to first region in while loop below
    current_region_index = -1
    current_region = (None, None, None)
    last_end = None
    for line in wiggletools_process.stdout:
        tokens = line.strip().split('\t')
        chrom, start, end, coverage = (
                tokens[0], int(tokens[1]), int(tokens[2]),
                int(float(tokens[3]))
            )
        if last_end != start:
            '''Must advance to next region of regions list iff we've reached a
            new region in region BED.'''
            while (chrom != current_region[0]
                    or start < current_region[1]
                    or end > current_region[2]):
                current_region_index += 1
                current_region = regions[current_region_index]
        aucs[current_region_index] += (end - start) * coverage
        last_end = end
    wiggletools_process.terminate()
    exit_code = wiggletools_process.wait()
    return (aucs, exit_code)

def star_aucs_from_regions(args):
    """ Converts aucs_from_regions([1,2]) to aucs_from_regions(1,2) 

        args: list of arguments

        Return value: aucs_from_regions(*args)
    """
    return aucs_from_regions(*args)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--regions', type=str, required=True,
            help='BED file with regions in which to compute AUCs; '
                 'sort this beforehand with '
                 '"sort -k1,1 -k2,2n -k3,3n <BED file>"'
        )
    parser.add_argument('--bams', type=str, required=False, nargs='+',
            help='space-separated list of paths to BAM files; all BAMs must '
                 'be indexed'
        )
    parser.add_argument('-m', '--bam-manifest', type=str, required=False,
            help='manifest file listing paths to BAMs and sample names'
        )
    parser.add_argument('--wiggletools', type=str, required=False,
            default='wiggletools',
            help='path to wiggletools executable'
        )
    parser.add_argument('--num-processes', '-p', type=int, required=False,
            default=1,
            help='maximum number of wiggletools processes to run at once'
        )
    args = parser.parse_args()
    if not (os.path.isfile(args.wiggletools) 
            and os.access(args.wiggletools, os.X_OK)):
        raise RuntimeError(
                'Cannot execute wiggletools. Check that the wiggletools '
                'binary exists at the specified path (i.e., the argument of '
                '--wiggletools) and is executable.'
            )
    if args.bams and args.bam_manifest:
        raise RuntimeError(
                'BAMs are specified with either --bams or --bam-manifest, but '
                'not both.'
            )
    elif not (args.bams or args.bam_manifest):
        raise RuntimeError(
                'BAMs must be specified as arguments of --bams or in a '
                'manifest file specified as the argument of --bam-manifest.'
            )
    if args.bams:
        # Use BAMs specified at command line
        bams = [(bam, bam) for bam in args.bams]
    else:
        # Read manifest file
        bams = []
        with open(args.bam_manifest) as bam_stream:
            for line in bam_stream:
                tokens = line.strip().split('\t')
                try:
                    bams.append((tokens[0], tokens[1]))
                except IndexError:
                    # No sample name specified
                    bams.append((tokens[0], tokens[0]))
    args.wiggletools = os.path.abspath(args.wiggletools)
    print >>sys.stderr, 'Reading regions...'
    # Store regions in interval tree to increment in case of overlaps
    regions = []
    with open(args.regions) as region_stream:
        for line in region_stream:
            tokens = line.strip().split('\t')
            chrom, start, end = tokens[0], int(tokens[1]), int(tokens[2])
            regions.append((chrom, start, end))
    regions.sort()
    print >>sys.stderr, 'Done. Computing region matrix...'
    pool = multiprocessing.Pool(args.num_processes, init_worker)
    output_matrix = []
    output_result = pool.map_async(
            star_aucs_from_regions,
            [[args.wiggletools, args.regions, bam[0], regions]
                for bam in bams],
            callback=output_matrix.extend
        )
    while not output_result.ready():
        time.sleep(3)
    pool.close()
    pool.terminate()
    pool.join()
    if any(result[1] for result in output_matrix):
        raise RuntimeError(
                'Nonzero exit code encountered computing output matrix.'
            )
    # Dump final matrix
    print '\t'.join([''] + ['{}:{}-{}'.format(*region) for region in regions])
    for i, bam in enumerate(bams):
        print '\t'.join([bam[1]] + map(str, output_matrix[i][0]))
