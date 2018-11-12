#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
#  Gap2Seq
#  Copyright (C) Leena Salmela, Kristoffer Sahlin, Veli Mäkinen,
#  Alexandru Tomescu, Riku Walve 2017
#
#  Contact: leena.salmela@cs.helsinki.fi
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import os, sys
import subprocess, multiprocessing
import datetime, random

# Helpers for finding executables
script_dir = os.path.dirname(os.path.realpath(__file__))
isexecutable = lambda f: os.path.isfile(f) and os.access(f, os.X_OK)
def find_executable(path_hints, name):
    for path_hint in path_hints:
        path = os.path.abspath(os.path.join(path_hint, name))
        if isexecutable(path):
            return path

    print('%s not found' % name, file=sys.stderr)
    sys.exit(1)

# Helper for checking if a file already exists
def check_file(filename):
    if os.path.isfile(filename):
        print('%s exists' % filename)
        sys.exit(1)

# Find required tools
GAPMERGER = find_executable([script_dir, './', '../'], 'GapMerger')
GAPCUTTER = find_executable([script_dir, './', '../'], 'GapCutter')
GAP2SEQ = find_executable([script_dir, './', '../'], 'Gap2Seq-core')
READFILTER = find_executable([script_dir, './', '../'], 'ReadFilter')

# An object for holding all the data for a library of short reads
class Library:
    def __init__(self, bam, mean_insert_size, std_dev, threshold):
        self.bam, self.mu, self.sd, self.threshold = \
            bam, mean_insert_size, std_dev, threshold

        # Check bam is indexed
        if not os.path.isfile(bam + '.bai'):
            print('%s.bai not found' % bam, file=sys.stderr)
            sys.exit(1)

    def filter_unmapped(self, filename):
        # TODO: Allow ReadFilter to run with just -unmapped-only
        subprocess.check_call([READFILTER,
            '-unmapped-only',
            '-scaffold', '0',
            '-breakpoint', '0',
            '-gap-length', '0',
            '-reads', filename] + self.data(),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def data(self):
        return ['-bam', self.bam,
            '-mean', str(self.mu),
            '-std-dev', str(self.sd)]

class Gap:
    def __init__(self, scaffold, position, gap_length, flank_length,
            left, right, comment, id, insertion):
        self.scaffold, self.position = scaffold, position
        self.left, self.right, self.gap_length = left, right, gap_length
        self.flank_length, self.comment, self.id = flank_length, comment, id
        self.insertion = insertion

    def data(self):
        length = str(self.gap_length)
        if self.insertion:
            length = '0'

        return ['-scaffold', str(self.scaffold),
            '-breakpoint', str(self.position),
            '-flank-length', str(self.flank_length),
            '-gap-length', length]

    def filler_data(self):
        return ['-left', self.left,
            '-right', self.right,
            '-length', str(self.gap_length)]

# Listener for filled gaps, to print progress
successful_gaps, filled_gaps, num_of_gaps = 0, 0, 1
def listener(queue, filename):
    global successful_gaps, filled_gaps, num_of_gaps
    start_time = datetime.datetime.now()
    with open(filename, 'w') as f:
        while 1:
            result = queue.get()
            if result == 'kill':
                break

            success, comment, fill = result

            filled_gaps += 1
            if success:
                successful_gaps += 1

            delta = (datetime.datetime.now() - start_time).total_seconds()
            eta = (delta / filled_gaps) * (num_of_gaps - filled_gaps)
            eta_string = str(datetime.timedelta(seconds=eta))

            print('Progress %.3f%% [%i / %i] %s left\t' % \
                (100*(filled_gaps / num_of_gaps), filled_gaps, num_of_gaps, eta_string),
                end='\r')

            f.write(comment + '\n' + fill + '\n')
            f.flush()
    return successful_gaps

# Runs all the read filtering and gap filling for a single gap
def fill_gap(libraries, gap, k, fuz, solid, derr, max_mem, randseed,
        upper, unique, best, reads=None, queue=None):
    # Cleanup, just to be sure
    reads_base = 'tmp.reads.' + gap.id + '.'
    subprocess.check_call(['rm', '-f', reads_base + '*'])

    # Filter reads
    filtered = False
    if reads == None:
        filtered = True
        reads = []
        with open('tmp.extract.' + gap.id + '.log', 'w') as f:
            filtered_length = 0
            for i, lib in enumerate(libraries):
                reads_file = reads_base + str(i)
                subprocess.check_call([READFILTER,
                    '-reads', reads_file] + \
                    gap.data() + lib.data(), stderr=f, stdout=f)

                # If no reads are filtered, no file exists
                if not os.path.isfile(reads_file):
                    continue

                reads.append(os.path.abspath(reads_file))

                # Count the number of filtered bases
                grep = subprocess.check_output('grep \'^[^>;]\' ' + reads_file + ' | wc -c', shell=True)
                filtered_length += int(grep)

            # Thresholding has to be done here in case of multiple libraries
            threshold = sum([lib.threshold for lib in libraries])
            if (filtered_length / gap.gap_length) < threshold:
                for i, lib in enumerate(libraries):
                    reads_file = 'tmp.reads.' + str(i) + '.unmapped'
                    if os.path.isfile(reads_file):
                        reads.append(os.path.abspath(reads_file))

    # Run Gap2Seq on the gap with the filtered reads
    log, filled = b'', 'tmp.filled.' + gap.id
    with open('tmp.gap2seq.' + gap.id + '.log', 'w') as f:
        wd = os.getcwd()
        wd_new = os.path.join(wd, 'tmp.' + str(random.randrange(1000, 9999)) + '.' + gap.id)
        subprocess.check_call(['mkdir', wd_new])
        os.chdir(wd_new)

        try:
            binary_options = [option for option, on in [('-all-upper', upper), ('-unique', unique), ('-best-only', best)] if on]
            log = subprocess.check_output([GAP2SEQ,
                '-k', str(k),
                '-fuz', str(fuz),
                '-solid', str(solid),
                '-nb-cores', '1',
                '-dist-error', str(derr),
                '-max-mem', str(max_mem),
                '-randseed', str(randseed),
                '-reads', ','.join(reads),
                '-filled', os.path.join(wd, filled)] + binary_options + gap.filler_data(),
                stderr=f)
        except subprocess.CalledProcessError:
            log = b''

        os.chdir(wd)
        subprocess.check_call(['rm', '-rf', wd_new])

    # Parse filled scaffold
    fill = ''
    if log != b'' and os.path.isfile(filled):
        with open(filled, 'r') as f:
            for line in f:
                if line[0] != '>':
                    fill += line.rstrip()

        # Remove logs and temporary/intermediate files
        subprocess.check_call(['rm', '-f',
            'tmp.filled.' + gap.id,
            'tmp.extract.' + gap.id + '.log',
            'tmp.gap2seq.' + gap.id + '.log'])
    else:
        fill = gap.left + 'N' * gap.gap_length + gap.right

    # Cleanup reads
    if filtered:
        subprocess.check_call(['rm', '-f', reads_base + '*'])

    successful = (not 'N' in fill) and (not 'n' in fill)
    if queue != None:
        queue.put((successful, gap.comment, fill))
    else:
        return (successful, gap.comment, fill)

# Run Gap2Seq without read filtering, i.e. just feed all scaffolds to it
# FIXME: This is an ugly way to implement this
def fill_scaffolds(scaffolds, filled, k, fuz, solid, derr, max_mem, randseed,
        upper, unique, best, reads, threads):

    if type(scaffolds) != type(''):
        scaffolds = scaffolds.name

    binary_options = [option for option, on in [('-all-upper', upper), ('-unique', unique), ('-best-only', best)] if on]
    command = [GAP2SEQ,
        '-k', str(k),
        '-fuz', str(fuz),
        '-solid', str(solid),
        '-nb-cores', str(threads),
        '-dist-error', str(derr),
        '-max-mem', str(max_mem),
        '-randseed', str(randseed),
        '-reads', ','.join(reads),
        '-filled', filled,
        '-scaffolds', scaffolds] + binary_options
    # print(' '.join(command))
    subprocess.check_call(command)

# NOTE: Assumes gaps and the bed file are in the same order
def parse_gap(bed, gap, id, insertion=False):
    gap = gap.split('\n')
    comment = gap[0]

    gap = ''.join(gap[1:])

    left = gap[:gap.upper().find('N')]
    right = gap[gap.upper().rfind('N')+1:]
    flank_length = min(len(left), len(right))
    gap_length = len(gap) - len(left) - len(right)

    # Parse gap data from bed file
    gap_data = bed.readline().rstrip().split('\t')
    scaffold = gap_data[0]
    position = int(gap_data[1]) + len(left)

    return Gap(scaffold, position, gap_length, flank_length, left, right, comment, id, insertion)

# Starts multiple gapfilling processes in parallel
def start_fillers(bed, gaps, libraries, queue=None, pool=None, k=31, fuz=10,
        solid=2, derr=500, max_mem=20, reads=None, randseed=0, upper=False,
        unique=False, best=False, insertion=False):
    start_filler = lambda seq, gap_id: fill_gap(libraries, parse_gap(bed, seq,
        str(gap_id), insertion), k, fuz, solid, derr, max_mem, randseed, upper,
        unique, best, reads)

    if pool != None:
        start_filler = lambda seq, gap_id: pool.apply_async(fill_gap,
            args=([libraries, parse_gap(bed, seq, str(gap_id), insertion), k,
                fuz, solid, derr, max_mem, randseed, upper, unique, best,
                reads, queue]))

    gap_id = 0

    jobs = []
    seq = ''
    for gap in gaps:
        if gap[0] == '>' and seq != '':
            jobs.append(start_filler(seq, gap_id))
            gap_id += 1
            seq = ''
        seq += gap

    jobs.append(start_filler(seq, gap_id))

    return jobs

# Run GapCutter, i.e. cut scaffolds into contigs and gaps
def cut_gaps(k, fuz, mask, scaffolds, contigs_file = 'tmp.contigs',
        gap_file = 'tmp.gaps', bed_file = 'tmp.bed'):
    if contigs_file != 'tmp.contigs': check_file(contigs_file)
    if gap_file != 'tmp.gaps': check_file(gap_file)
    if bed_file != 'tmp.bed': check_file(bed_file)

    binary_options = [option for option, on in [('-mask', mask)] if on]
    command = [GAPCUTTER,
        '-k', str(k),
        '-fuz', str(fuz),
        '-scaffolds', scaffolds,
        '-gaps', gap_file,
        '-contigs', contigs_file,
        '-bed', bed_file] + binary_options
    # print(' '.join(command))
    subprocess.check_call(command)

    return open(bed_file, 'r'), open(gap_file, 'r')

# Run GapMerger, i.e. merge contigs and filled gaps back into scaffolds
def merge_gaps(filled, merged, contigs_file='tmp.contigs'):
    command = [GAPMERGER,
        '-scaffolds', merged,
        '-gaps', filled,
        '-contigs', contigs_file]
    # print(' '.join(command))
    subprocess.check_call(command)

# Parse VCF file and extract kmers from reference genome
def cut_vcf(vcf, reference_file, k, fuz, contigs_file = 'tmp.contigs',
        gap_file = 'tmp.gaps', bed_file = 'tmp.bed'):
    # Parse chromosomes from the reference into a dictionary
    # NOTE: This loads the entire reference genome
    reference = {}
    comment, seq = '', ''
    for line in f:
        if line[0] == '>':
            if seq != '':
                reference[comment] = seq
                seq = ''
            comment = line.rstrip()[1:]
        seq += line.rstrip()
    reference[comment] = seq

    # Open a new gap and bed files
    gap, bed = open(gap_file, 'r+'), open(bed_file, 'r+')
    for line in f:
        if line[0] == '#': continue
        fields = line.rstrip().split('\t')

        # Parse VCF fields
        insert_ref = fields[3]
        insert = fields[4][len(insert_ref):]
        comment, start, end = fields[0], int(fields[1]) - 1, \
            (int(fields[1]) - 1) + (len(insert) - 1)

        # Extract k-mers from reference seqeuences
        left_end = min(start+k+fuz, len(reference[comment]))
        right_start = max(end-(k+fuz), 0)
        left = reference[comment][start:left_end]
        right = reference[comment][right_start:end]
        seq = left + 'N' * len(insert) + right

        gap.write('>%s:%i-%i\n%s' % (comment, start, end, seq))
        bed.write('%s\t%i\t%i' % (comment, start-(k+fuz)+1, end+(k+fuz)+1))

    # We will read the files later, so they need to be seeked to beginning
    gap.seek(0)
    bed.seek(0)

    return gap, bed

# Count the number of gaps
def count_gaps(bed):
    global num_of_gaps
    num_of_gaps = 0

    bed.seek(0)
    for line in bed:
        num_of_gaps += 1
    bed.seek(0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Gap2Seq 3.1')

    # wrapper specific arguments
    parser.add_argument('-f', '--filled', required=True, type=str, help="output file for filled scaffolds")
    parser.add_argument('-t', '--threads', type=int, default=1, help="number of threads to use")
    parser.add_argument('--insertion', action='store_true', help=argparse.SUPPRESS)

    # Gap2Seq specific arguments
    parser.add_argument('-k', type=int, default=31, help="k-mer length for DBG  [default 31]")
    parser.add_argument('--fuz', type=int, default=10, help="number of nucleotides to ignore on gap fringes  [default 10]")
    parser.add_argument('--solid', type=int, default=2, help="threshold for solid k-mers for building the DBG [default 2]")
    parser.add_argument('--max-mem', type=float, default=20, help="maximum memory usage of DP table computation in gigabytes (excluding DBG) [default 20]")
    parser.add_argument('--dist-error', type=int, default=500, help="maximum error in gap estimates  [default 500]")
    parser.add_argument('--randseed', type=int, default=0, help="random seed (0 to use current time)  [default 0]")
    parser.add_argument('--all-upper', action='store_true', help="fill all bases in upper case.")
    parser.add_argument('--unique', action='store_true', help="fill only gaps with a unique path of best length")
    parser.add_argument('--best-only', action='store_true', help="consider only paths that have optimal length")
    parser.add_argument('--mask', action='store_true', help=argparse.SUPPRESS)

    # Reads can be given in either a set of fasta-formatted reads
    parser.add_argument('-r', '--reads', type=str, help="short reads")

    # or a set of aligned read libraries
    # Tab-separated list: bam, mean_insert_size, std_dev, threshold
    parser.add_argument('-l', '--libraries', type=str, help="list of aligned read libraries")
    parser.add_argument('-i', '--index', type=int, default=-1, help=argparse.SUPPRESS)

    # One of three options is required for gap data:
    # 1. Cut gaps and bed from some scaffolds
    parser.add_argument('-s', '--scaffolds', type=str, help="scaffolds with gaps")

    # 2. Use pre-cut gaps and bed
    parser.add_argument('-b', '--bed', type=argparse.FileType('r'), help="pre-cut gaps")
    parser.add_argument('-g', '--gaps', type=argparse.FileType('r'), help="")

    # 3. Generate gaps and bed from VCF
    parser.add_argument('-v', '--vcf', type=argparse.FileType('r'), help="variants in reads against reference")
    parser.add_argument('-R', '--reference', type=argparse.FileType('r'), help="")

    args = vars(parser.parse_args())

    # Short read libraries aligned to the scaffolds
    libraries = []
    if args['libraries'] != None:
        with open(args['libraries'], 'r') as f:
            for lib in f:
                arg = lib.split('\t')
                libraries += [Library(arg[0], int(arg[1]), int(arg[2]), int(arg[3]))]

        # Use only 1 library
        if args['index'] != -1:
            libraries = [libraries[args['index']]]

        # Filter unmapped reads
        for i, lib in enumerate(libraries):
            lib.filter_unmapped('tmp.reads.' + str(i) + '.unmapped')
    elif args['reads'] != None:
        args['reads'] = args['reads'].split(',')
    else:
        parser.print_help()
        print('Either [-r/--reads], or [-l/--libraries] is required.')
        exit(1)

    scaffolds_cut = False
    if args['bed'] == None or args['gaps'] == None:
        if args['scaffolds'] != None:
            print('Cutting gaps')
            args['bed'], args['gaps'] = cut_gaps(args['k'], args['fuz'], args['mask'], args['scaffolds'])
            scaffolds_cut = True
            args['final_out'] = args['filled']
            args['filled'] = 'tmp.filled'
        elif args['vcf'] != None:
            print('Parsing VCF')
            args['bed'], args['gaps'] = cut_vcf(args['vcf'], args['reference'], args['k'], args['fuz'])
            args['insertion'] = True
        else:
            parser.print_help()
            print('Either [-s/--scaffolds], [-b/--bed and -g/--gaps], or [-v/--vcf and -R/--reference] are required.')
            exit(1)

    count_gaps(args['bed'])

    if args['libraries'] != None:
        # Gap2Seq divides the max mem evenly between threads, but as we run multiple
        # parallel instances with 1 thread, we need to pre-divide
        max_mem = args['max_mem'] / args['threads']

        queue, pool = None, None
        if args['threads'] > 1:
            manager = multiprocessing.Manager()
            queue = manager.Queue()
            pool = multiprocessing.Pool(args['threads'] + 1)

        # Start listening for filled gaps
        if args['threads'] > 1:
            res = pool.apply_async(listener, (queue, args['filled']))

        print('Starting gapfillers')
        jobs = start_fillers(args['bed'], args['gaps'], libraries, queue=queue,
            pool=pool, k=args['k'], fuz=args['fuz'], solid=args['solid'],
            derr=args['dist_error'], max_mem=max_mem, reads=args['reads'],
            randseed=args['randseed'], upper=args['all_upper'],
            unique=args['unique'], best=args['best_only'],
            insertion=args['insertion'])

        args['bed'].close()
        args['gaps'].close()

        if args['threads'] > 1:
            for job in jobs:
                job.get()
            queue.put('kill')
            successful_gaps = res.get(timeout=1)
            pool.close()
            pool.join()

        # Cleanup unmapped reads
        subprocess.check_call(['rm', '-f', 'tmp.reads.*.unmapped'])
    else:
        fill_scaffolds(args['gaps'], args['filled'], args['k'], args['fuz'],
            args['solid'], args['dist_error'], args['max_mem'],
            args['randseed'], args['all_upper'], args['unique'],
            args['best_only'], args['reads'], args['threads'])

    if scaffolds_cut:
        print('Merging filled gaps and contigs')
        merge_gaps(args['filled'], args['final_out'])

    if args['libraries'] != None:
        print('Filled %i out of %i gaps' % (successful_gaps, num_of_gaps))
