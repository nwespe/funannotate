#!/usr/bin/env python

import sys, multiprocessing, subprocess, os, shutil, argparse, time, inspect
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='prodigal_parallel.py', usage="%(prog)s [options] -i genome.fasta",
    description='''Script runs prodigal in parallel to use multiple processors''',
    epilog="""Modified from augustus-parallel.py written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', default='prodigal.gff3', help='Name of output file')
parser.add_argument('-p','--proteins', default='prodigal.proteins.fa', help='Name of proteins output file')
parser.add_argument('-f','--format', default='gff', help='Output file format')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to run')
parser.add_argument('--debug', action='store_true', help='Keep intermediate files')
parser.add_argument('--logfile', default ='prodigal-parallel.log', help='logfile')
args=parser.parse_args()


def countGFFgenes(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count


def runProdigal(Input):
    if '_part' in Input:
        chr = Input.split('_part')[0]
    else:
        chr = Input
    prod_out = os.path.join(tmpdir, Input+'.prodigal.gff3')
    core_cmd = ['prodigal', '-i', os.path.join(tmpdir, chr+'.fa'), '-a', os.path.join(tmpdir, chr+'.proteins.faa'), '-f', args.format]
    #try using library module
    lib.runSubprocess2(core_cmd, '.', lib.log, prod_out)


log_name = args.logfile
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

#first step is to split input fasta file into individual files in tmp folder
lib.log.debug("Splitting contigs file")
tmpdir = 'prodigal_tmp_'+str(os.getpid())
os.makedirs(tmpdir)
scaffolds = []
global ranges
ranges = {}
with open(args.input, 'rU') as InputFasta:
    for record in SeqIO.parse(InputFasta, 'fasta'):  # removed code to split up large contigs
        name = str(record.id)
        scaffolds.append(name)
        outputfile = os.path.join(tmpdir, name+'.fa')
        with open(outputfile, 'w') as output:
            SeqIO.write(record, output, 'fasta')

#now loop through each scaffold running prodigal
if args.cpus > len(scaffolds):
    num = len(scaffolds)
else:
    num = args.cpus
lib.log.debug("Running Prodigal on %i chunks, using %i CPUs" % (len(scaffolds), num))
lib.runMultiProgress(runProdigal, scaffolds, num)


lib.log.debug("Prodigal prediction is finished, now concatenating results")
with open(args.out, 'w') as output:
    for scaffold in scaffolds:
        scaf_file = os.path.join(tmpdir, scaffold + '.prodigal.gff3')
        scaf_id = scaffold.lstrip('sequence')
        with open(scaf_file) as input:
            for line in input:
                line = line.replace('CDS', 'gene')   # substitute gene for CDS
                line = line.replace('ID=1_', 'ID=' + scaf_id + '_')  # substitute scaffold number for ID=X
                output.write(line)

with open(args.proteins, 'w') as output:
    for scaffold in scaffolds:
        prot_file = os.path.join(tmpdir, scaffold + '.proteins.faa')
        scaf_id = scaffold.lstrip('sequence')
        with open(prot_file) as input:
            output.write(input.read().replace('ID=1_', 'ID=' + scaf_id + '_'))

if not args.debug:
    shutil.rmtree(tmpdir)
lib.log.info('Found {0:,}'.format(countGFFgenes(args.out))+' gene models')
