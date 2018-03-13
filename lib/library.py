from __future__ import division
import os, subprocess, logging, sys, argparse, inspect, csv, time, re, shutil, datetime, glob, platform, multiprocessing, itertools, hashlib, math, types, gzip
from natsort import natsorted
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# get the working directory, so you can move back into DB folder to find the files you need
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
LIB = os.path.join(parentdir, 'lib')
UTIL = os.path.join(parentdir, 'util')
GeneMark2GFF = os.path.join(UTIL, 'genemark_gtf2gff3.pl')

pref_colors=["#CF3C57","#65B23A","#6170DD","#D18738","#D542B5",
"#724A63","#60AABA","#5DB07C","#6C5824","#D74B2B","#6B97D6","#893B2E",
"#B68DB7","#564E91","#ACA13C","#3C6171","#436B33","#D84088",
"#D67A77","#9D55C4","#8B336E","#DA77B9","#D850E5","#B188DF"]

Nogs = {'NOG': 'All organisms (5.0GB)',
'aciNOG': 'Acidobacteria (125.3MB)',
'acidNOG': 'Acidobacteriia (75.4MB)',
'acoNOG': 'Aconoidasida (217.1MB)',
'actNOG': 'Actinobacteria (765.3MB)',
'agaNOG': 'Agaricales (211.1MB)',
'agarNOG': 'Agaricomycetes (236.5MB)',
'apiNOG': 'Apicomplexa (322.7MB)',
'aproNOG': 'Proteobacteria_alpha (638.4MB)',
'aquNOG': 'Aquificae (51.5MB)',
'arNOG': 'Archaea (256.9MB)',
'arcNOG': 'Archaeoglobi (21.8MB)',
'artNOG': 'Arthropoda (725.0MB)',
'arthNOG': 'Arthrodermataceae (111.2MB)',
'ascNOG': 'Ascomycota (1.1GB)',
'aveNOG': 'Aves (186.1MB)',
'bacNOG': 'Bacilli (362.6MB)',
'bactNOG': 'Bacteria (3.3GB)',
'bacteNOG': 'Bacteroidia (199.2MB)',
'basNOG': 'Basidiomycota (356.5MB)',
'bctoNOG': 'Bacteroidetes (508.9MB)',
'biNOG': 'Bilateria (1.7GB)',
'bproNOG': 'Proteobacteria_beta (481.0MB)',
'braNOG': 'Brassicales (275.4MB)',
'carNOG': 'Carnivora (293.5MB)',
'chaNOG': 'Chaetomiaceae (180.9MB)',
'chlNOG': 'Chlorobi (51.3MB)',
'chlaNOG': 'Chlamydiae (39.1MB)',
'chloNOG': 'Chloroflexi (136.8MB)',
'chlorNOG': 'Chloroflexi (75.8MB)',
'chloroNOG': 'Chlorophyta (146.8MB)',
'chorNOG': 'Chordata (1.1GB)',
'chrNOG': 'Chromadorea (392.6MB)',
'cloNOG': 'Clostridia (505.6MB)',
'cocNOG': 'Coccidia (137.4MB)',
'creNOG': 'Crenarchaeota (110.0MB)',
'cryNOG': 'Cryptosporidiidae (105.4MB)',
'cyaNOG': 'Cyanobacteria (254.8MB)',
'cytNOG': 'Cytophagia (164.6MB)',
'debNOG': 'Debaryomycetaceae (145.5MB)',
'defNOG': 'Deferribacteres (41.6MB)',
'dehNOG': 'Dehalococcoidetes (15.0MB)',
'deiNOG': 'Deinococcusthermus (75.4MB)',
'delNOG': 'delta/epsilon (471.4MB)',
'dipNOG': 'Diptera (397.7MB)',
'dotNOG': 'Dothideomycetes (298.2MB)',
'dproNOG': 'Proteobacteria_delta (424.6MB)',
'droNOG': 'Drosophilidae (314.1MB)',
'eproNOG': 'Proteobacteria_epsilon (104.8MB)',
'eryNOG': 'Erysipelotrichi (85.8MB)',
'euNOG': 'Eukaryotes (3.1GB)',
'eurNOG': 'Euryarchaeota (264.7MB)',
'euroNOG': 'Eurotiomycetes (507.2MB)',
'eurotNOG': 'Eurotiales (358.1MB)',
'fiNOG': 'Fishes (641.2MB)',
'firmNOG': 'Firmicutes (728.8MB)',
'flaNOG': 'Flavobacteriia (222.5MB)',
'fuNOG': 'Fungi (1.2GB)',
'fusoNOG': 'Fusobacteria (74.9MB)',
'gproNOG': 'Proteobacteria_gamma (735.0MB)',
'haeNOG': 'Haemosporida (197.1MB)',
'halNOG': 'Halobacteria (106.6MB)',
'homNOG': 'Hominidae (229.9MB)',
'hymNOG': 'Hymenoptera (199.5MB)',
'hypNOG': 'Hypocreales (353.3MB)',
'inNOG': 'Insects (688.9MB)',
'kinNOG': 'Kinetoplastida (259.8MB)',
'lepNOG': 'Lepidoptera (208.0MB)',
'lilNOG': 'Liliopsida (660.0MB)',
'maNOG': 'Mammals (855.5MB)',
'magNOG': 'Magnaporthales (161.3MB)',
'meNOG': 'Animals (1.8GB)',
'metNOG': 'Methanobacteria (38.4MB)',
'methNOG': 'Methanococci (24.5MB)',
'methaNOG': 'Methanomicrobia (99.4MB)',
'necNOG': 'Nectriaceae (200.3MB)',
'negNOG': 'Negativicutes (96.5MB)',
'nemNOG': 'Nematodes (430.0MB)',
'onyNOG': 'Onygenales (282.8MB)',
'opiNOG': 'Opisthokonts (2.8GB)',
'perNOG': 'Peronosporales (154.1MB)',
'plaNOG': 'Planctomycetes (149.3MB)',
'pleNOG': 'Pleosporales (223.4MB)',
'poaNOG': 'Poales (596.3MB)',
'prNOG': 'Primates (448.8MB)',
'proNOG': 'Proteobacteria (1.5GB)',
'rhaNOG': 'Rhabditida (334.5MB)',
'roNOG': 'Rodents (381.4MB)',
'sacNOG': 'Saccharomycetaceae (202.7MB)',
'saccNOG': 'Saccharomycetes (275.9MB)',
'sorNOG': 'Sordariales (296.1MB)',
'sordNOG': 'Sordariomycetes (714.1MB)',
'sphNOG': 'Sphingobacteriia (154.0MB)',
'spiNOG': 'Spirochaetes (121.2MB)',
'spriNOG': 'Supraprimates (635.6MB)',
'strNOG': 'Streptophyta (960.6MB)',
'synNOG': 'Synergistetes (59.5MB)',
'tenNOG': 'Tenericutes (29.9MB)',
'thaNOG': 'Thaumarchaeota (15.3MB)',
'theNOG': 'Thermoplasmata (26.9MB)',
'therNOG': 'Thermotogae (66.5MB)',
'thermNOG': 'Thermococci (31.4MB)',
'treNOG': 'Tremellales (79.9MB)',
'veNOG': 'Vertebrates (1.0GB)',
'verNOG': 'Verrucomicrobia (140.9MB)',
'verrNOG': 'Verrucomicrobiae (73.0MB)',
'virNOG': 'Viridiplantae (1.0GB)'}

COGS = {'J': '(J) Translation, ribosomal structure and biogenesis',
'A': '(A) RNA processing and modification',
'K': '(K) Transcription',
'L': '(L) Replication, recombination and repair',
'B': '(B) Chromatin structure and dynamics',
'D': '(D) Cell cycle control, cell division, chromosome partitioning',
'Y': '(Y) Nuclear structure',
'V': '(V) Defense mechanisms',
'T': '(T) Signal transduction mechanisms',
'M': '(M) Cell wall/membrane/envelope biogenesis',
'N': '(N) Cell motility',
'Z': '(Z) Cytoskeleton',
'W': '(W) Extracellular structures',
'U': '(U) Intracellular trafficking, secretion, and vesicular transport',
'O': '(O) Posttranslational modification, protein turnover, chaperones',
'C': '(C) Energy production and conversion',
'G': '(G) Carbohydrate transport and metabolism',
'E': '(E) Amino acid transport and metabolism',
'F': '(F) Nucleotide transport and metabolism',
'H': '(H) Coenzyme transport and metabolism',
'I': '(I) Lipid transport and metabolism',
'P': '(P) Inorganic ion transport and metabolism',
'Q': '(Q) Secondary metabolites biosynthesis, transport and catabolism',
'R': '(R) General function prediction only',
'S': '(S) Function unknown'}

DBURL = { 'uniprot': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
        'uniprot-release': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt',
        'merops': 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsscan.lib',
        'dbCAN': 'http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt',
        'dbCAN-tsv': 'http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt',
        'dbCAN-log': 'http://csbl.bmb.uga.edu/dbCAN/download/readme.txt',
        'pfam': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.hmm.gz',
        'pfam-tsv': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.clans.tsv.gz',
        'pfam-log': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam.version.gz',
        'outgroups': 'https://osf.io/r9sne/download?version=1',
        'repeats': 'https://osf.io/vp87c/download?version=1',
        'go-obo': 'http://purl.obolibrary.org/obo/go.obo', 
        'mibig': 'http://mibig.secondarymetabolites.org/MIBiG_prot_seqs_1.3.fasta',
        'interpro': 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',
        'gene2product': 'https://raw.githubusercontent.com/nextgenusfs/gene2product/master/ncbi_cleaned_gene_products.txt'}


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'


class gzopen(object):
   """Generic opener that decompresses gzipped files
   if needed. Encapsulates an open file or a GzipFile.
   Use the same way you would use 'open()'.
   """
   def __init__(self, fname):
      f = open(fname)
      # Read magic number (the first 2 bytes) and rewind.
      magic_number = f.read(2)
      f.seek(0)
      # Encapsulated 'self.f' is a file or a GzipFile.
      if magic_number == '\x1f\x8b':
         self.f = gzip.GzipFile(fileobj=f)
      else:
         self.f = f

   # Define '__enter__' and '__exit__' to use in
   # 'with' blocks. Always close the file and the
   # GzipFile if applicable.
   def __enter__(self):
      return self
   def __exit__(self, type, value, traceback):
      try:
         self.f.fileobj.close()
      except AttributeError:
         pass
      finally:
         self.f.close()

   # Reproduce the interface of an open file
   # by encapsulation.
   def __getattr__(self, name):
      return getattr(self.f, name)
   def __iter__(self):
      return iter(self.f)
   def next(self):
      return next(self.f)

def Funzip(input, output, cpus):
    '''
    function to unzip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '--decompress', '-c', '-p', str(cpus), input]
    elif which('bgzip'):
        cmd = ['bgzip', '--decompress', '-c', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '--decompress', '-c', input]
    try:
        runSubprocess2(cmd, '.', log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip(input, output, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '-c', '-p', str(cpus), input]
    elif which('bgzip'):
        cmd = ['bgzip', '-c', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '-c', input]
    try:
        runSubprocess2(cmd, '.', log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip_inplace(input):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    cpus = multiprocessing.cpu_count()
    if which('pigz'):
        cmd = ['pigz', '-f', '-p', str(cpus), input]
    elif which('bgzip'):
        cmd = ['bgzip', '-f', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        runSubprocess(cmd, '.', log)
    except NameError:
        subprocess.call(cmd)


def CheckFASTQandFix(forward, reverse):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from itertools import izip, izip_longest
    # open and check first header, if okay exit, if not fix
    file1 = FastqGeneralIterator(gzopen(forward))
    file2 = FastqGeneralIterator(gzopen(reverse))
    check = True
    for read1, read2 in izip(file1, file2):
        # see if index is valid
        if ' ' in read1[0] and ' ' in read2[0]:
            if read1[0].split(' ')[1].startswith('1') and read2[0].split(' ')[1].startswith('2'):  # std illumina, exit
                break
        elif read1[0].endswith('/1') and read2[0].endswith('/2'):  # also acceptable
            break
        else:  # it is not okay missing paired information
            check = False
            break
    file1.close()
    file2.close()
    if not check:  # now need to fix these reads
        log.info("PE reads do not conform to Trinity naming convention (need either /1 /2 or std illumina), fixing...")     
        # work on forward reads first
        if forward.endswith('.gz'):
            Funzip(forward, forward+'.bak', multiprocessing.cpu_count())
            SafeRemove(forward)
        else:
            os.rename(forward, forward+'.bak')
        # now add ending to reads
        with open(forward+'.fix', 'w') as forwardfix:
            for title, seq, qual in FastqGeneralIterator(open(forward+'.bak')):
                title = title+'/1'
                forwardfix.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        Fzip(forward+'.fix', forward, multiprocessing.cpu_count())      
        SafeRemove(forward+'.bak')
        SafeRemove(forward+'.fix')           
        # now work on reverse reads
        if reverse.endswith('.gz'):
            Funzip(reverse, reverse+'.bak', multiprocessing.cpu_count())
        else:            
            os.rename(reverse, reverse+'.bak')
        with open(reverse+'.fix', 'w') as reversefix:
            for title, seq, qual in FastqGeneralIterator(open(reverse+'.bak')):
                title = title+'/2'
                reversefix.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        # zip back up to original file
        Fzip(reverse+'.fix', reverse, multiprocessing.cpu_count())
        SafeRemove(reverse+'.bak')
        SafeRemove(reverse+'.fix')
    return


def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return


def runSubprocess(cmd, dir, logfile):
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout:
        logfile.debug(stdout)
    if stderr:
        logfile.debug(stderr)


def runSubprocess2(cmd, dir, logfile, output):
    # function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)

def runSubprocess3(cmd, dir, logfile):
    #function where STDOUT pipes to FNULL, capture STDERR in logfile
    FNULL = open(os.devnull, 'w')
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=FNULL, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        logfile.debug(stderr)
        
def runSubprocess4(cmd, dir, logfile):
    #function where STDOUT and STDERR pipes to FNULL
    FNULL = open(os.devnull, 'w')
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=FNULL, stderr=FNULL)
    proc.communicate()

def runSubprocess5(cmd, dir, logfile, input, output):
    #function where STDOUT to file, STDIN as input, STDERR pipes to logfile
    logfile.debug(' '.join(cmd))
    with open(input) as infile:
        with open(output, 'w') as out:
            proc = subprocess.Popen(cmd, cwd=dir, stdin=infile, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)

def runSubprocess6(cmd, dir, logfile, logfile2):
    #function where cmd captured in logfile, but both stdout and stdin piped to additional logfile
    logfile.debug(' '.join(cmd))
    with open(logfile2, 'w') as logout:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stdout:
            logout.write(stdout)
        if stderr:
            logout.write(stderr)

def evmGFFvalidate(input, evmpath, logfile):
    Validator = os.path.join(evmpath, 'EvmUtils', 'gff3_gene_prediction_file_validator.pl')
    cmd = [Validator, input]
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if not stderr:
        return True
    else:
        logfile.debug(stderr)
        False

def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.digest()
    
def sha256_check(file1, file2):
    files = [file1, file2]
    output = [(fname, hashfile(open(fname, 'rb'), hashlib.sha256())) for fname in files]
    if output[0][1] == output[1][1]:
        return True
    else:
        return False

def readBlocks(source, pattern):
    buffer = []
    for line in source:
        if line.startswith(pattern):
            if buffer: yield buffer
            buffer = [ line ]
        else:
            buffer.append( line )
    yield buffer
    
def empty_line_sep(line):
    return line=='\n'

def get_parent_dir(directory):
    return os.path.dirname(directory)

def getSize(filename):
    st = os.stat(filename)
    return st.st_size
    
def checkinputs(filename):
    if not os.path.isfile(filename):
        log.error("%s is not a valid file, exiting" % filename)
        sys.exit(1)
    size = getSize(filename)
    if size < 2: #this is 1 character...
        log.error("%s appears to be empty, exiting" % filename)
        sys.exit(1)

def make_tarfile(output_filename, source_dir):
    import tarfile
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

def multipleReplace(text, wordDict):
    for key in wordDict:
        text = text.replace(key, wordDict[key])
    return text

def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'mafft', 'signalp', 'proteinortho5.pl', 'ete3', 'phyml', 'phobius.pl']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
            else:
                if name == 'signalp':
                    subprocess.Popen([name, '-V'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'dustmasker':
                    subprocess.Popen([name, '-version-full'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'tbl2asn':
                    subprocess.Popen([name, '--help'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'raxmlHPC-PTHREADS':
                    subprocess.Popen([name, '-version'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'ete3':
                    subprocess.Popen([name, 'version'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'phobius.pl':
                    subprocess.Popen([name, '-h'], stdout=devnull, stderr=devnull).communicate()
                else:
                    subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)
        
def checkannotations(input):
    if os.path.isfile(input):
        filesize = getSize(input)
        if int(filesize) < 1:
            return False
        else:
            return True
    else:
        return False

def line_count(fname):
    with open(fname) as f:
        i = -1
        for i, l in enumerate(f):
            pass
    return i + 1

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def getGeneBasename(fastafile):
    bases = []
    with open(fastafile, 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('>'):
                line = line.replace('>', '')
                Base = line.split('_')[0]+'_'
                if not Base in bases:
                    bases.append(Base)
    return bases

def get_version():
    cmd = [os.path.join(parentdir, 'funannotate'), 'version']
    version = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0].rstrip()
    return version

def checkAugustusFunc(base):
    '''
    fucntion to try to test Augustus installation is working, note segmentation fault still results in a pass
    '''
    brakerpass = 0
    buscopass = 0
    version = subprocess.Popen(['augustus', '--version'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
    version = version.split(' is ')[0]
    bam2hints = which(os.path.join(base, 'bin', 'bam2hints'))
    filterBam = which(os.path.join(base, 'bin', 'filterBam'))
    if bam2hints and filterBam:
        brakerpass = 1
    model = os.path.join(parentdir, 'lib', 'EOG092C0B3U.prfl')
    if not os.path.isfile(model):
        log.error("Testing Augustus Error: installation seems wrong, can't file prfl model")
        sys.exit(1)
    profile = '--proteinprofile='+model
    proteinprofile = subprocess.Popen(['augustus', '--species=anidulans', profile, os.path.join(parentdir, 'lib', 'busco_test.fa')], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
    proteinprofile.strip()
    if proteinprofile == '':
        buscopass = 0
    elif not 'augustus: ERROR' in proteinprofile:
        buscopass = 1
    return (version, brakerpass, buscopass)
    
def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

def fmtcols(mylist, cols):
    justify = []
    for i in range(0,cols):
        length = max(map(lambda x: len(x), mylist[i::cols]))
        length += 2
        ljust = map(lambda x: x.ljust(length), mylist[i::cols])
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) / cols
    lines = (' '.join(justify[i::num_lines]) 
             for i in range(0,num_lines))
    return "\n".join(lines)

def list_columns(obj, cols=4, columnwise=True, gap=4):
    """
    Print the given list in evenly-spaced columns.

    Parameters
    ----------
    obj : list
        The list to be printed.
    cols : int
        The number of columns in which the list should be printed.
    columnwise : bool, default=True
        If True, the items in the list will be printed column-wise.
        If False the items in the list will be printed row-wise.
    gap : int
        The number of spaces that should separate the longest column
        item/s from the next column. This is the effective spacing
        between columns based on the maximum len() of the list items.
    """

    sobj = [str(item) for item in obj]
    if cols > len(sobj): cols = len(sobj)
    max_len = max([len(item) for item in sobj])
    if columnwise: cols = int(math.ceil(float(len(sobj)) / float(cols)))
    plist = [sobj[i: i+cols] for i in range(0, len(sobj), cols)]
    if columnwise:
        if not len(plist[-1]) == cols:
            plist[-1].extend(['']*(len(sobj) - len(plist[-1])))
        plist = zip(*plist)
    printer = '\n'.join([
        ''.join([c.ljust(max_len + gap) for c in p])
        for p in plist])
    return printer


def roundup(x):
    return x if x % 100 == 0 else x + 100 - x % 100
    
def maxabs(a, axis=None):
    """Return slice of a, keeping only those values that are furthest away
    from 0 along axis"""
    maxa = a.max(axis=axis)
    mina = a.min(axis=axis)
    p = abs(maxa) > abs(mina) # bool, or indices where +ve values win
    n = abs(mina) > abs(maxa) # bool, or indices where -ve values win
    if axis == None:
        if p: return maxa
        else: return mina
    shape = list(a.shape)
    shape.pop(axis)
    out = np.zeros(shape, dtype=a.dtype)
    out[p] = maxa[p]
    out[n] = mina[n]
    return out
    
def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%I:%M:%S %p]')
    else:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%I:%M:%S %p]')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def renameGFF(input, newname, output):
    with open(output, 'w') as outfile:
        with open(input, 'rU') as infile:
            for line in infile:
                if line.startswith('>'): #remove any fasta sequences
                    continue
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    cols = line.split('\t')
                    #make sure it has correct columns to be GFF
                    if len(cols) == 9:
                        outfile.write('%s\t%s\t%s' % (cols[0], newname, '\t'.join(cols[2:])))               

def countGFFgenes(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count

def countEVMpredictions(input):
    augustus = 0
    genemark = 0
    pasa = 0
    hiq = 0
    other = 0
    total = 0
    with open(input, 'rU') as f:
        for line in f:
            line = line.strip()
            contig, source, feature, start, end, blank, strand, score, info = line.split('\t')
            if feature == 'gene':
                total += 1
                if source == 'Augustus':
                    augustus += 1
                elif source == 'GeneMark':
                    genemark += 1
                elif source == 'pasa_pred':
                    pasa += 1
                elif source == 'other_pred':
                    other += 1
                elif source == 'HiQ':
                    hiq += 1
    return total, augustus, genemark, hiq, pasa, other

def countGMAPtranscripts(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith('###'):
                count += 1
    return count
    
def runMultiProgress(function, inputList, cpus):
    #setup pool
    p = multiprocessing.Pool(cpus)
    #setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    #refresh pbar every 5 seconds
    while True:
        incomplete_count = sum(1 for x in results if not x.ready())
        if incomplete_count == 0:
            break
        sys.stdout.write("     Progress: %.2f%% \r" % (float(tasks - incomplete_count) / tasks * 100))
        sys.stdout.flush()
        time.sleep(1)
    p.close()
    p.join()

def runMultiNoProgress(function, inputList, cpus):
    #setup pool
    p = multiprocessing.Pool(cpus)
    #setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    p.close()
    p.join()

def cleanProteins(inputList, output):
    #expecting a list of protein fasta files for combining/cleaning headers
    #make sure you aren't duplicated sequences names
    #dropping proteins less than 50 amino acids
    seen = set()
    with open(output, 'w') as out:
        for x in inputList:
            with open(x, 'rU') as input:
                for rec in SeqIO.parse(input, 'fasta'):
                    if len(rec.seq) < 50:
                        continue
                    #explicitly check for swissprot and jgi
                    if rec.id.startswith('sp|') or rec.id.startswith('jgi|'):
                        ID = rec.id.split('|')[-1]
                    else:
                        ID = rec.id
                    #now clean up the shit
                    badshit = [':', ';', '/', '\\', '.', ',', '%']
                    for i in badshit:
                        if i in ID:
                            ID = ID.replace(i, '_')
                    if not ID in seen:
                        seen.add(ID)
                    else:
                        #means that ID has already been used, so add a number to it, auto increment
                        counter = 1
                        while ID in seen:
                            oldnum = counter-1
                            ID = ID.replace('_'+str(oldnum), '') + '_'+str(counter)
                            counter += 1
                        seen.add(ID)
                    out.write('>%s\n%s\n' % (ID, rec.seq))
  
def gb2output(input, output1, output2, output3):
    with open(output1, 'w') as proteins:
        with open(output2, 'w') as transcripts:
            with open(output3, 'w') as scaffolds:
                with open(input, 'rU') as gbk:
                    SeqRecords = SeqIO.parse(gbk, 'genbank')
                    for record in SeqRecords:
                        scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                        for f in record.features:
                            if f.type == "CDS":
                                proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0].rstrip('*')))
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))

def sortGFF(input, output, order):
    cmd = ['bedtools', 'sort', '-header', '-faidx', order, '-i', input]
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] == None:
            if stderr[1] != '':
                log.error("Sort GFF failed, unreferenced scaffold present in gene predictions, check logfile")
                sys.exit(1)
 
def checkGenBank(input):
    count = 0
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    count += 1
    if count == 0:
        return False
    else:
        return True
        
def countGenBank(input):
    cds = 0
    trna = 0
    dnas = 0
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            dnas += 1
            for f in record.features:
                if f.type == 'CDS':
                    cds += 1
                elif f.type == 'tRNA':
                    trna += 1
    return dnas, cds, trna
        
def checkFastaHeaders(input, limit):
    length = 0
    names = []
    with open(input, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line.replace('\n', '')
                ID = line.replace('>', '').strip()
                names.append(ID)
                headlen = len(line) - 1 #subtract one character for fasta carrot 
                if headlen > length:
                    length = headlen
    if length > int(limit):
        return (False, names)
    else:
        return (True, names)

def BamHeaderTest(genome, mapping):
    import pybam
    #get list of fasta headers from genome
    genome_headers = []
    with open(genome, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id not in genome_headers:
                genome_headers.append(rec.id)
    #get list of fasta headers from BAM
    bam_headers = []
    with open(mapping, 'rb') as bamin:
        bam = pybam.read(bamin)
        bam_headers = bam.file_chromosomes
    #make sure the bam headers is a list
    if not type(bam_headers) is list:
        log.error("PyBam parsing failed, printing results, funannotate is expecting a list, not this....")
        print bam_headers
        sys.exit(1)
    #now compare lists, basically if BAM headers not in genome headers, then output bad names to logfile and return FALSE
    genome_headers = set(genome_headers)
    diffs = [x for x in bam_headers if x not in genome_headers]
    if len(diffs) > 0:
        log.debug("ERROR: These BAM headers not found in genome FASTA headers\n%s" % ','.join(diffs))
        return False
    else:
        return True

def convertgff2tbl(gff, gffskip, prefix, fasta, proteins, tblout):
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    '''
    function to convert directly from gff to gbl
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['start'])
    GeneCount = 0
    #first get the scaffold names and lengths
    scaffLen = {}
    with open(fasta, 'rU') as seqin:
        for record in SeqIO.parse(seqin, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    Genes = {}
    ParentDB = {}
    orphans = []
    invalid = []
    with open(gff, 'rU') as infile:
        for line in infile:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            ID, Parent, Product, Note = (None,)*4
            if feature == 'gene':
                GeneCount += 1
                ID = attributes.split(';')[0].replace('ID=', '')
                if not ID in Genes:
                    Genes[ID] = {'type': '', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': False, 'proper_stop': False, 'phase': [], 'product': 'hypothetical protein', 'note': '' } 
                else:
                    print("Duplicate Gene IDs found, %s" % ID)
            else: #meaning needs to append to a gene ID as it is mRNA, tRNA, CDS, exon
                info = attributes.split(';')
                for x in info:
                    if x.startswith('ID='):
                        ID = x.replace('ID=', '')
                    elif x.startswith('Parent='):
                        Parent = x.replace('Parent=', '')
                        if '-' in Parent:
                            Parent = Parent.split('-')[0]
                    elif x.startswith('product=') or x.startswith('Product='):
                        Product = x.split('roduct=')[-1]
                    elif x.startswith('note=') or x.startswith('Note='):
                        Note = x.split('ote=')[-1]
                if not Parent:
                    invalid.append(line)
                    continue
                if not Parent in Genes: #GFF file not sorted? try to recapture
                    if Parent in ParentDB:
                        Parent = ParentDB.get(Parent)
                    else:
                        orphans.append(line)                
                if Product:
                    Genes[Parent]['product'] = Product
                if Note:
                    Genes[Parent]['note'] = Note
                if feature == 'mRNA':
                    Genes[Parent]['type'] = 'mRNA'
                    if not ID in ParentDB:
                        ParentDB[ID] = Parent
                elif feature == 'tRNA':
                    Genes[Parent]['type'] = 'tRNA'
                if feature == 'exon':
                    if Parent in Genes:
                        Genes[Parent]['mRNA'].append((start, end))
                elif feature == 'CDS':
                    if Parent in Genes:
                        Genes[Parent]['CDS'].append((start, end))
                        Genes[Parent]['phase'].append(phase)
    if len(invalid) > 0:
        log.error("{:,} GFF features missing Parent/ID annotations and were skipped".format(len(invalid)))
        with open(gffskip, 'w') as invalidgff:
            for x in invalid:
                invalidgff.write('%s\n' % x)
    if len(orphans) > 0:
        log.error("GFF file parsing error, child features found without valid Parent. Check child/parent naming scheme.")
        log.debug('%s' % '\n'.join(orhpans))
        sys.exit(1)
    #now sort dictionary by contig and location
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    SeqRecords = SeqIO.index(fasta, 'fasta')
    counter = 1
    with open(proteins, 'w') as protout:
        for k,v in sortedGenes.items():
            if not prefix:
                locusTag = k
            else:
                locusTag = prefix+'_'+str(counter).zfill(6) 
            if v['strand'] == '+':
                sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0])
                sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0], reverse=True)
                sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0], reverse=True)
            v['mRNA'] = sortedExons
            v['CDS'] = sortedCDS
            renamedGenes[locusTag] = v
            if v['type'] == 'mRNA':
                if len(v['CDS']) == 0: #there is no CDS, drop this "model"
                    log.debug('%s has no CDS, skipping model' % k)
                    continue
                #get the codon_start by getting first CDS phase + 1
                indexStart = [x for x, y in enumerate(v['CDS']) if y[0] == sortedCDS[0][0]]
                codon_start = int(v['phase'][indexStart[0]]) + 1
            else:
                codon_start = ''
            renamedGenes[locusTag]['codon_start'] = codon_start
            if not v['contig'] in scaff2genes:
                scaff2genes[v['contig']] = [locusTag]
            else:
                scaff2genes[v['contig']].append(locusTag)
            #translate to protein sequence, construct cDNA Seq object, translate
            if v['strand'] == '+' and v['type'] == 'mRNA':
                cdsSeq = ''
                for s in v['CDS']:
                    singleCDS = SeqRecords[v['contig']][s[0]-1:s[1]]
                    cdsSeq += (str(singleCDS.seq))
                mySeq = Seq(cdsSeq, IUPAC.ambiguous_dna)
                protSeq = mySeq.translate(cds=False, table=1)
            elif v['strand'] == '-' and v['type'] == 'mRNA':
                cdsSeq = ''
                for s in sorted(v['CDS'], key=lambda tup: tup[0]):
                    singleCDS = SeqRecords[v['contig']][s[0]-1:s[1]]
                    cdsSeq += (str(singleCDS.seq))
                mySeq = Seq(cdsSeq, IUPAC.ambiguous_dna)
                mySeq = mySeq.reverse_complement()
                protSeq = mySeq.translate(cds=False, table=1)
            if v['type'] == 'mRNA':
                if protSeq.endswith('*'):
                    protStop = True
                    protout.write('>%s\n%s\n' % (locusTag,protSeq[:-1]))
                else:
                    protStop = False
                    protout.write('>%s\n%s\n' % (locusTag,protSeq))
                if protSeq.startswith('M'):
                    protStart = True
                else:
                    protStart = False
                renamedGenes[locusTag]['proper_start'] = protStart
                renamedGenes[locusTag]['proper_stop'] = protStop
            counter += 1
             
    #now have scaffolds dict and gene dict, loop through scaff dict printing tbl
    with open(tblout, 'w') as tbl:
        for k,v in natsorted(scaff2genes.items()):
            tbl.write('>Feature %s\n' % k)
            tbl.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
            tbl.write('\t\t\t%s\t%s\n' % ('CFMR', '12345'))
            for genes in v: #now loop through each gene on the scaffold
                geneInfo = renamedGenes.get(genes)
                if geneInfo['type'] == 'mRNA':
                    #check for partial models
                    if geneInfo['proper_start'] and geneInfo['codon_start'] == 1:
                        partialStart = ''
                    else:
                        partialStart = '<'
                    if geneInfo['proper_stop']:
                        partialStop = ''
                    else:
                        partialStop = '>'
                    if geneInfo['strand'] == '+':
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['start'], partialStop, geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[0], partialStop, exon[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[0], exon[1]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[0], partialStop, exon[1]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[0], partialStop, cds[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[0], cds[1]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[0], partialStop, cds[1]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)                                      
                    else:
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['end'], partialStop, geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[1], partialStop, exon[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[1], exon[0]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[1], partialStop, exon[0]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[1], exon[0]))                 
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[1], partialStop, cds[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[1], cds[0]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[1], partialStop, cds[0]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)
                elif geneInfo['type'] == 'tRNA':
                    if geneInfo['strand'] == '+':
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['start'], geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        tbl.write('<%s\t>%s\ttRNA\n' % (geneInfo['start'], geneInfo['end']))
                        if geneInfo['product'] != 'hypothetical protein':
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            tbl.write('\t\t\tnote\t%s\n' % geneInfo['note'])                                    
                    else:
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['end'], geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        tbl.write('<%s\t>%s\ttRNA\n' % (geneInfo['end'], geneInfo['start']))
                        if geneInfo['product'] != 'hypothetical protein':
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            tbl.write('\t\t\tnote\t%s\n' % geneInfo['note'])
    return GeneCount

def GFF2tbl(evm, trnascan, proteins, scaffLen, prefix, Numbering, SeqCenter, SeqRefNum, tblout):
    from collections import OrderedDict
    '''
    function to take EVM protein models and tRNA scan GFF to produce a GBK tbl file as well
    as a new GFF3 file. The function will also rename locus_id if passed.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['start'])
    Proteins = {}
    with open(proteins, 'rU') as prots:
        for record in SeqIO.parse(prots, 'fasta'):
            ID = record.id.replace('evm.model', 'evm.TU')
            start, stop = (True,)*2
            Seq = str(record.seq)
            if not Seq.endswith('*'):
                stop = False
            if not Seq.startswith('M'):
                start = False
            if not ID in Proteins:
                Proteins[ID] = {'start': start, 'stop': stop}    
    Genes = {}
    with open(evm, 'rU') as infile:
        for line in infile:
            if line.startswith('\n'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            if feature == 'gene':
                ID = attributes.split(';')[0].replace('ID=', '')
                if not ID in Genes:
                    Genes[ID] = {'type': 'mRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': Proteins[ID]['start'], 'proper_stop': Proteins[ID]['stop'], 'phase': [], 'product': 'hypothetical protein', 'note': '' }
                else:
                    print("Duplicate Gene IDs found, %s" % ID)
            else: #meaning needs to append to a gene ID
                info = attributes.split(';')
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                Parent = Parent.replace('evm.model', 'evm.TU')
                if feature == 'exon':
                    if Parent in Genes:
                        Genes[Parent]['mRNA'].append((start, end))
                elif feature == 'CDS':
                    if Parent in Genes:
                        Genes[Parent]['CDS'].append((start, end))
                        Genes[Parent]['phase'].append(phase)
    #now load tRNA predictions
    with open(trnascan, 'rU') as trnain:
        for line in trnain:
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            if feature == 'gene':
                ID = attributes.split(';')[0].replace('ID=', '')
                if not ID in Genes:
                    Genes[ID] = {'type':'tRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': False, 'proper_stop': False, 'phase': [], 'product': '', 'note': ''  }
            elif feature == 'tRNA':
                info = attributes.split(';')
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                Product = info[2].replace('product=', '')
                Note = info[3].replace('note=', '')
                Genes[Parent]['product'] = Product
                Genes[Parent]['note'] = Note
            elif feature == 'exon':
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                if Parent in Genes:
                    Genes[Parent]['mRNA'].append((start, end))
                    Genes[Parent]['CDS'].append((start, end))
                    Genes[Parent]['phase'].append('0') 
    #now sort dictionary by contig and location, rename using prefix
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    count = Numbering
    for k,v in sortedGenes.items():
        if prefix:
            locusTag = prefix+'_'+str(count).zfill(6)
        else:
            locusTag = k
        if v['strand'] == '+':
            sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0])
            sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0])
        else:
            sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0], reverse=True)
            sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0], reverse=True)
        #get the codon_start by getting first CDS phase + 1
        indexStart = [x for x, y in enumerate(v['CDS']) if y[0] == sortedCDS[0][0]]
        codon_start = int(v['phase'][indexStart[0]]) + 1
        v['mRNA'] = sortedExons
        v['CDS'] = sortedCDS
        renamedGenes[locusTag] = v
        renamedGenes[locusTag]['codon_start'] = codon_start
        if not v['contig'] in scaff2genes:
            scaff2genes[v['contig']] = [locusTag]
        else:
            scaff2genes[v['contig']].append(locusTag)
        count += 1
    #now have scaffolds dict and gene dict, loop through scaff dict printing tbl
    with open(tblout, 'w') as tbl:
        for k,v in natsorted(scaff2genes.items()):
            tbl.write('>Feature %s\n' % k)
            tbl.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
            tbl.write('\t\t\t%s\t%s\n' % (SeqCenter, SeqRefNum))
            for genes in v: #now loop through each gene on the scaffold
                geneInfo = renamedGenes.get(genes)
                if geneInfo['type'] == 'mRNA':
                    #check for partial models
                    if geneInfo['proper_start'] and geneInfo['codon_start'] == 1:
                        partialStart = ''
                    else:
                        partialStart = '<'
                    if geneInfo['proper_stop']:
                        partialStop = ''
                    else:
                        partialStop = '>'
                    if geneInfo['strand'] == '+':
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['start'], partialStop, geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[0], partialStop, exon[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[0], exon[1]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[0], partialStop, exon[1]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[0], partialStop, cds[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[0], cds[1]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[0], partialStop, cds[1]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)                                      
                    else:
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['end'], partialStop, geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[1], partialStop, exon[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[1], exon[0]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[1], partialStop, exon[0]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[1], exon[0]))                 
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[1], partialStop, cds[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[1], cds[0]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[1], partialStop, cds[0]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)
                elif geneInfo['type'] == 'tRNA':
                    if geneInfo['strand'] == '+':
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['start'], geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0:
                                tbl.write('<%s\t>%s\ttRNA\n' % (exon[0], exon[1]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            tbl.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            tbl.write('\t\t\tnote\t%s\n' % geneInfo['note'])                                    
                    else:
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['end'], geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0:
                                tbl.write('<%s\t>%s\ttRNA\n' % (exon[1], exon[0]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            tbl.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            tbl.write('\t\t\tnote\t%s\n' % geneInfo['note']) 

def getGBKinfo(input):
    accession = None
    organism = None
    strain = None
    isolate = None
    gb_gi = None
    WGS_accession = None
    version = None
    with open(input, 'rU') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            try:
                WGS_accession = 'WGS:'+record.annotations['contig'].split(':')[0].replace('join(', '')[:4]
            except KeyError:
                pass
            try:
                accession = record.annotations['accessions'][0]
            except KeyError:
                pass
            try:
                organism = record.annotations['organism'].replace('Unclassified.', '').rstrip()
            except KeyError:
                pass
            try:
                gb_gi = record.annotations['gi']
            except KeyError:
                pass
            try:
                version = record.annotations['sequence_version']
            except KeyError:
                pass
            for f in record.features:
                if f.type == "source":
                    isolate = f.qualifiers.get("isolate", [None])[0]
                    strain = f.qualifiers.get("strain", [None])[0]
            break
    return organism, strain, isolate, accession, WGS_accession, gb_gi, version

def getGBKLocusTag(input):
    LocusTags = []
    with open(input, 'rU') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            for f in record.features:
                if f.type == 'gene':
                    ID = f.qualifiers['locus_tag'][0]
                    if not ID in LocusTags:
                        LocusTags.append(ID)
    lastTag = natsorted(LocusTags)[-1]
    tag, count = lastTag.split('_')
    justify = len(count)
    return tag, count, justify

def gb2dna(input, output):
    with open(output, 'w') as outfile:
        with open(input, 'rU') as infile:
            for record in SeqIO.parse(infile, 'genbank'):
                outfile.write(">%s\n%s\n" % (record.id, record.seq))

def getID(input, type):
    #function to get ID from genbank record.features
    locusTag = None
    ID = None
    Parent = None
    if type == 'gene':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['gene'][0]
            except KeyError:
                pass
        return locusTag, ID, locusTag
        
    elif type == 'mRNA' or type == 'tRNA':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['transcript_id'][0]
                ID = locusTag
            except KeyError:
                pass
            try:
                Parent = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['transcript_id'][0]
            except KeyError:
                pass 
        return locusTag, ID, Parent
                   
    elif type == 'CDS':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['protein_id'][0]
            except KeyError:
                pass
            try:
                Parent = input.qualifiers['gene'][0]
            except KeyError:
                pass       
        else:
            try:
                ID = input.qualifiers['protein_id'][0]
            except KeyError:
                pass
        return locusTag, ID, Parent
        
def gb2parts(input, tbl, prots, trans, dna):
    from collections import OrderedDict
    '''
    function returns a dictionary of all gene models from a genbank file this function
    can handle multiple transcripts per locus/gene
    '''
    GeneCount = 0
    genes = {}
    transcript_parts = {}
    protein_parts = {}
    scaffolds = {} # contig : {'length': len, 'loci': [list]}
    with open(dna, 'w') as dnaout:
        with open(prots, 'w') as protout:
            with open(trans, 'w') as transout:
                with open(input, 'rU') as filein:
                    for record in SeqIO.parse(filein, 'genbank'):
                        dnaout.write(">%s\n%s\n" % (record.id, record.seq))
                        Contig = record.id
                        if not Contig in scaffolds:
                            scaffolds[Contig] = {'length': len(record.seq), 'loci': []}
                        for f in record.features:
                            #reset for every feature
                            ID,start,end,strand,num_parts,exons,cds,protID,transcriptID,protSeq,product,note,sortedExons,sortedCDS = (None,)*14
                            Fivepartial = ''
                            Threepartial = ''
                            if not f.type in ['gene', 'mRNA', 'tRNA', 'CDS']:
                                continue
                            #get info from features
                            ID, featureID, Parent = getID(f, f.type)
                            start = f.location.nofuzzy_start + 1 
                            end = f.location.nofuzzy_end
                            strand = f.location.strand
                            if strand == 1:
                                str = '+'
                            else:
                                str = '-'
                            num_parts = len(f.location.parts)           
                            if f.type == "gene":
                                GeneCount += 1
                                if unicode(f.location.start).startswith('<'):
                                    Fivepartial = '<'
                                if unicode(f.location.end).startswith('>'):
                                    Threepartial = '>'
                                if not ID in scaffolds[Contig]['loci']:
                                    scaffolds[Contig]['loci'].append(ID)
                                if not ID in genes:
                                    genes[ID] = {'type': None, 'contig': Contig, 'location': (int(start),int(end)), 'strand': str, 'transcripts': [], 'proteins': [], '5partial': Fivepartial, '3partial': Threepartial, 'product': 'hypothetical protein', 'note': None}
                                else:
                                    print("Duplicate Gene ID: %s" % ID)
                                    genes[ID]['location'] = (int(start),int(end))
                                    genes[ID]['strand'] = strand            
                            elif f.type == "mRNA":
                                try:
                                    product = f.qualifiers['product'][0]
                                except KeyError:
                                    product = 'hypothetical protein'
                                exonTuples = []
                                if num_parts < 2: #only single exon
                                    exonTuples.append((int(start),int(end)))
                                else: #more than 1 exon, so loop through
                                    for i in range(0, num_parts):
                                        ex_start = f.location.parts[i].nofuzzy_start + 1
                                        ex_end = f.location.parts[i].nofuzzy_end
                                        exonTuples.append((int(ex_start),int(ex_end)))
                                #now we want to sort the positions I think...
                                if strand == 1:
                                    sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
                                else:
                                    sortedExons = sorted(exonTuples, key=lambda tup: tup[0], reverse=True)
                                #update dictionary
                                if featureID:
                                    transID = featureID
                                else:
                                    transID = ID
                                if not Parent in genes:
                                    genes[Parent] = {'type': f.type, 'contig': Contig, 'location': (int(start),int(end)), 'strand': str, 'transcripts': [transID], 'proteins': [], 'product': product, 'note': None}
                                else:
                                    genes[Parent]['transcripts'].append(transID)
                                    genes[Parent]['type'] = f.type
                                    genes[Parent]['product'] = product
                                if not transID in transcript_parts:
                                    transcript_parts[transID] = {'mRNA': sortedExons, 'strand': str, 'location': (int(start),int(end))}
                                else:
                                    print('Duplicate transcript: %s' % ID)
                                transout.write('>%s %s\n%s\n' % (Parent, transID, f.extract(record.seq)))
                            elif f.type == 'tRNA':
                                try:
                                    product = f.qualifiers['product'][0]
                                except KeyError:
                                    product = 'tRNA-Xxx'
                                try:
                                    note = f.qualifiers['note'][0]
                                except KeyError:
                                    note = None
                                if strand == 1:
                                    sortedExons = (int(start), int(end))
                                else:
                                    sortedExons = (int(end), int(start))
                                if featureID:
                                    transID = featureID
                                else:
                                    transID = ID
                                if not Parent in genes:
                                    genes[Parent] = {'type': f.type, 'contig': Contig, 'location': (int(start),int(end)), 'strand': str, 'transcripts': [transID], 'proteins': [], 'product': product, 'note': note}                                                                  
                                else:
                                    genes[Parent]['transcripts'].append(transID)
                                    genes[Parent]['type'] = f.type
                                    genes[Parent]['product'] = product
                                    genes[Parent]['note'] = note
                                if not transID in transcript_parts:
                                    transcript_parts[transID] = {'mRNA': sortedExons, 'strand': str, 'location': (int(start),int(end))}
                                else:
                                    print('Duplicate transcript: %s' % ID)
                            elif f.type == 'CDS':
                                protSeq = f.qualifiers['translation'][0]
                                cdsTuples = []
                                phase = f.qualifiers['codon_start'][0]
                                if num_parts < 2: #only single CDS
                                    cdsTuples.append((int(start),int(end)))
                                else:
                                    for i in range(0, num_parts):
                                        ex_start = f.location.parts[i].nofuzzy_start + 1 
                                        ex_end = f.location.parts[i].nofuzzy_end
                                        cdsTuples.append((int(ex_start),int(ex_end)))
                                if strand == 1:
                                    sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
                                else:
                                    sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0], reverse=True)
                                #check for annotations
                                try:
                                    product = f.qualifiers['product'][0]
                                except KeyError:
                                    product = 'hypothetical protein'
                                note = []
                                dbxref = []
                                go_terms = []                                
                                #update dictionary
                                if featureID:
                                    protID = featureID
                                else:
                                    protID = ID
                                if not Parent in genes:
                                    genes[Parent] = {'contig': Contig, 'location': (int(start),int(end)), 'strand': str, 'transcripts': [], 'proteins': [protID], 'product': product}
                                else:
                                    genes[Parent]['proteins'].append(protID)
                                    genes[Parent]['product'] = product
                                if not protID in protein_parts:
                                    protein_parts[protID] = {'CDS': sortedCDS, 'strand': str, 'location': (int(start),int(end)), 'seq': protSeq, 'codon_start': int(phase), 'dbxref': dbxref, 'note': note, 'go_ontology': go_terms}
                                else:
                                    print('Duplicate protein: %s' % ID)
                                protout.write('>%s %s\n%s\n' % (Parent, protID, protSeq))  
    #now loop through each scaffold generating tbl and gffoutput
    sScaff = sorted(scaffolds.iteritems(), key=lambda x: x[1]['length'], reverse=True)
    orderedScaff = OrderedDict(sScaff)
    with open(tbl, 'w') as tblout:
        for k,v in orderedScaff.items():
            tblout.write('>Feature %s\n' % k)
            tblout.write('1\t%i\tREFERENCE\n' % scaffolds[k]['length'])
            tblout.write('\t\t\t%s\t%s\n' % ('CFMR', '12345'))
            for g in v['loci']: #now loop through each gene on the scaffold
                geneInfo = genes.get(g)
                if geneInfo['type'] == 'mRNA':
                    if geneInfo['strand'] == '+':
                        tblout.write('%s%i\t%s%i\tgene\n' % (geneInfo['5partial'], geneInfo['location'][0], geneInfo['3partial'], geneInfo['location'][1]))
                        tblout.write('\t\t\tlocus_tag\t%s\n' % g)
                    else:
                        tblout.write('%s%i\t%s%i\tgene\n' % (geneInfo['5partial'], geneInfo['location'][1], geneInfo['3partial'], geneInfo['location'][0]))
                        tblout.write('\t\t\tlocus_tag\t%s\n' % g)
                    #now add mRNA and CDS
                    for i in range(0,len(geneInfo['transcripts'])):              
                        proteinID = geneInfo['proteins'][i]
                        transcriptID = geneInfo['transcripts'][i]
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(transcript_parts[transcriptID]['mRNA']):
                                if num == 0 and num == len(transcript_parts[transcriptID]['mRNA']) - 1: #single exon, so slightly differnt method
                                    tblout.write('%s%s\t%s%s\tmRNA\n' % (geneInfo['5partial'], exon[0], geneInfo['3partial'], exon[1]))
                                elif num == 0:
                                    tblout.write('%s%s\t%s\tmRNA\n' % (geneInfo['5partial'], exon[0], exon[1]))
                                elif num == len(transcript_parts[transcriptID]['mRNA']) - 1: #this is last one
                                    tblout.write('%s\t%s%s\n' % (exon[0], geneInfo['3partial'], exon[1]))
                                else:
                                    tblout.write('%s\t%s\n' % (exon[0], exon[1]))
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tblout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (g, i+1))
                            for num, cds in enumerate(protein_parts[proteinID]['CDS']):
                                if num == 0 and num == len(protein_parts[proteinID]['CDS']) - 1: #single exon, so slightly differnt method
                                    tblout.write('%s%s\t%s%s\tCDS\n' % (geneInfo['5partial'], cds[0], geneInfo['3partial'], cds[1]))
                                elif num == 0:
                                    tblout.write('%s%s\t%s\tCDS\n' % (geneInfo['5partial'], cds[0], cds[1]))
                                elif num == len(protein_parts[proteinID]['CDS']) - 1: #this is last one
                                    tblout.write('%s\t%s%s\n' % (cds[0], geneInfo['3partial'], cds[1]))
                                else:
                                    tblout.write('%s\t%s\n' % (cds[0], cds[1]))
                            tblout.write('\t\t\tcodon_start\t%i\n' % protein_parts[proteinID]['codon_start'])
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tblout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (g, i+1))                                      
                        else:
                            for num, exon in enumerate(transcript_parts[transcriptID]['mRNA']):
                                if num == 0 and num == len(transcript_parts[transcriptID]['mRNA']) - 1: #single exon, so slightly differnt method
                                    tblout.write('%s%s\t%s%s\tmRNA\n' % (geneInfo['5partial'], exon[1], geneInfo['3partial'], exon[0]))
                                elif num == 0:
                                    tblout.write('%s%s\t%s\tmRNA\n' % (geneInfo['5partial'], exon[1], exon[0]))
                                elif num == len(transcript_parts[transcriptID]['mRNA']) - 1: #this is last one
                                    tblout.write('%s\t%s%s\n' % (exon[1], geneInfo['3partial'], exon[0]))
                                else:
                                    tblout.write('%s\t%s\n' % (exon[1], exon[0]))                 
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tblout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (g, i+1))
                            for num, cds in enumerate(protein_parts[proteinID]['CDS']):
                                if num == 0 and num == len(protein_parts[proteinID]['CDS']) - 1: #single exon, so slightly differnt method
                                    tblout.write('%s%s\t%s%s\tCDS\n' % (geneInfo['5partial'], cds[1], geneInfo['3partial'], cds[0]))
                                elif num == 0:
                                    tblout.write('%s%s\t%s\tCDS\n' % (geneInfo['5partial'], cds[1], cds[0]))
                                elif num == len(protein_parts[proteinID]['CDS']) - 1: #this is last one
                                    tblout.write('%s\t%s%s\n' % (cds[1], geneInfo['3partial'], cds[0]))
                                else:
                                    tblout.write('%s\t%s\n' % (cds[1], cds[0]))
                            tblout.write('\t\t\tcodon_start\t%i\n' % protein_parts[proteinID]['codon_start'])
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tblout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (g, i+1))
                elif geneInfo['type'] == 'tRNA':
                    for i in range(0,len(geneInfo['transcripts'])):              
                        transcriptID = geneInfo['transcripts'][i]
                        if geneInfo['strand'] == '+':
                            tblout.write('%s%i\t%s%i\tgene\n' % (geneInfo['5partial'], geneInfo['location'][0], geneInfo['3partial'], geneInfo['location'][1]))
                            tblout.write('\t\t\tlocus_tag\t%s\n' % g)

                            tblout.write('<%i\t>%i\ttRNA\n' % (geneInfo['location'][0], geneInfo['location'][1]))
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tblout.write('\t\t\tpseudo\n')
                            if geneInfo['note']:
                                tblout.write('\t\t\tnote\t%s\n' % geneInfo['note'])                             
                        else:
                            tblout.write('%s%i\t%s%i\tgene\n' % (geneInfo['5partial'], geneInfo['location'][1], geneInfo['3partial'], geneInfo['location'][0]))
                            tblout.write('\t\t\tlocus_tag\t%s\n' % g)
                            tblout.write('<%i\t>%i\ttRNA\n' % (geneInfo['location'][1], geneInfo['location'][0]))
                            tblout.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tblout.write('\t\t\tpseudo\n')
                            if geneInfo['note']:
                                tblout.write('\t\t\tnote\t%s\n' % geneInfo['note'])   
    return GeneCount

def gb2allout(input, GFF, Proteins, Transcripts, DNA):
    #this will not output any UTRs for gene models, don't think this is a problem right now....
    errors = []
    with open(GFF, 'w') as gff:
        gff.write("##gff-version 3\n")
        with open(Proteins, 'w') as proteins:
            with open(Transcripts, 'w') as transcripts:
                with open(DNA, 'w') as scaffolds:
                    with open(input, 'rU') as gbk:
                        for record in SeqIO.parse(gbk, 'genbank'):
                            scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                            for f in record.features:
                                #get info from features
                                if f.type == 'gene' or f.type == 'mRNA' or f.type == 'CDS' or f.type == 'tRNA':
                                    locusTag, ID, Parent = getID(f, f.type)
                                    if not locusTag:
                                        continue
                                else:
                                    continue
                                strand = f.location.strand
                                if strand == 1:
                                    strand = '+'
                                elif strand == -1:
                                    strand = '-'
                                start = f.location.nofuzzy_start + 1
                                end = f.location.nofuzzy_end
                                chr = record.id
                                #now do some type specific routines
                                if f.type == "gene":
                                    if ID:
                                        gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s;\n" % (chr, start, end, strand, locusTag, ID))
                                    else:
                                        gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s;\n" % (chr, start, end, strand, locusTag))
                                if f.type == "mRNA":
                                    feature_seq = f.extract(record.seq)
                                    if ID:
                                        transcripts.write(">%s %s\n%s\n" % (locusTag, ID, feature_seq))
                                    else:
                                        transcripts.write(">%s %s\n%s\n" % (locusTag, Parent, feature_seq))
                                if f.type == 'CDS':
                                    if ID:
                                        proteins.write(">%s %s\n%s\n" % (locusTag, ID, f.qualifiers['translation'][0].rstrip('*')))
                                    else:
                                        proteins.write(">%s %s\n%s\n" % (locusTag, Parent, f.qualifiers['translation'][0].rstrip('*')))
                                    try:
                                        product = f.qualifiers['product'][0]
                                    except KeyError:
                                        product = "hypothetical protein"
                                    num_exons = len(f.location.parts)
                                    current_phase = int(f.qualifiers['codon_start'][0]) - 1 #need to adjust NCBI to GFF3 notation
                                    gff.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s;\n" % (chr, start, end, strand, locusTag, Parent, product))
                                    if num_exons < 2: #only a single exon
                                        ex_start = str(f.location.nofuzzy_start + 1)
                                        ex_end = str(f.location.nofuzzy_end)
                                        gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s.exon1;Parent=%s-T1;\n" % (chr, ex_start, ex_end, strand, locusTag, locusTag))
                                        gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t0\tID=%s.cds;Parent=%s-T1;\n" % (chr, ex_start, ex_end, strand, locusTag, locusTag))
                                    else: #more than 1 exon, so parts are actually in correct orientation, so loop through
                                        for i in range(0,num_exons):
                                            ex_start = str(f.location.parts[i].nofuzzy_start + 1)
                                            ex_end = str(f.location.parts[i].nofuzzy_end)
                                            ex_num = i + 1
                                            gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s.exon%i;Parent=%s-T1;\n" % (chr, ex_start, ex_end, strand, locusTag, ex_num, locusTag))
                                            gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s.cds;Parent=%s-T1;\n" % (chr, ex_start, ex_end, strand, current_phase, locusTag, locusTag))
                                            current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
                                            if current_phase == 3:
                                                current_phase = 0

                                if f.type == 'tRNA':
                                    try:
                                        product = f.qualifiers['product'][0]
                                    except KeyError:
                                        product = "tRNA-XXX"
                                    chr = record.id
                                    gff.write("%s\tGenBank\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s;\n" % (chr, start, end, strand, locusTag, Parent, product))
                                    gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s.exon1;Parent=%s-T1;\n" % (chr, start, end, strand, locusTag, locusTag))
    if len(errors) > 0:
        log.debug("No Translation in GBK file: %s" % ','.join(errors))
        
def runGMAP(transcripts, genome, cpus, intron, tmpdir, output):
    #first build genome database
    build_log = os.path.join(tmpdir, 'gmap-build.log')
    with open(build_log, 'w') as logfile:
        subprocess.call(['gmap_build', '-D', tmpdir, '-d', 'genome', '-k', '13', genome], stdout = logfile, stderr = logfile)
    #now map transcripts
    map_log = os.path.join(tmpdir, 'gmap-map.log')
    with open(map_log, 'w') as logfile:
        with open(output, 'w') as out:
            subprocess.call(['gmap', '--cross-species', '-f', '3', '-K', str(intron), '-n', '1', '-t', str(cpus), '-B', '5', '-D', tmpdir, '-d', 'genome', transcripts], stdout = out, stderr = logfile)
    
def runBUSCO(input, Database, cpus, tmpdir, output):
    #run busco in protein mapping mode
    BUSCO = os.path.join(UTIL, 'funannotate-BUSCO2.py')
    cmd = [BUSCO, '-i', input, '-m', 'proteins', '-l', Database, '-o', 'busco', '-c', str(cpus), '-f']
    runSubprocess(cmd, tmpdir, log)
    #now parse output and write to annotation file
    with open(output, 'w') as out:
        with open(os.path.join(tmpdir, 'run_busco', 'full_table_busco.tsv'), 'rU') as busco:
            for line in busco:
                if line.startswith('#'):
                    continue
                col = line.split('\t')
                if col[1] == 'Complete' or col[1] == 'Duplicated': #if diploid these should show up, but problematic for drawing trees....
                    if col[2].endswith('-T1'):
                        ID = col[2]
                    else:
                        ID = col[2]+'-T1'
                    out.write("%s\tnote\tBUSCO:%s\n" % (ID, col[0]))

def dupBUSCO2gff(ID, base_folder, locationID):
    hmmerfolder = os.path.join(base_folder, 'hmmer_output')
    geneID = ''
    AugFile = ''
    GFFfile = os.path.join(base_folder, 'augustus_output', 'gffs', ID+'.gff')
    if geneID == '':
        for file in os.listdir(hmmerfolder):
            if file.startswith(ID):
                with open(os.path.join(hmmerfolder, file), 'rU') as hmmer:
                    for line in hmmer:
                        if not line.startswith('#'):
                            longID = line.split()[0]
                            longID = longID.replace(']', '')
                            partsID = longID.split('[')
                            if locationID == partsID[1]:
                                geneID = partsID[0]
                                AugFile = os.path.join(base_folder, 'augustus_output', 'predicted_genes', file)
                                break
    #so now should have gene name, get the GFF from augustus
    with open(GFFfile, 'w') as gffout:
        with open(AugFile, 'rU') as augustus:
            for pred in readBlocks(augustus, '# start gene'):
                if pred[0].startswith('# This output'):
                    continue
                if pred[0].startswith('##gff-version 3'):
                    continue
                if pred[0].startswith('# Please cite'):
                    continue
                if geneID in pred[0]:
                    for x in pred:
                        if not x.startswith('#'):
                            gffout.write(x)
     
                  
def parseBUSCO2genome(input, ploidy, ContigSizes, output):
    #input is BUSCO output, ploidy is integer, ContigSizes is dictionary, output is a bedfile, function returns dictionary
    busco_complete = {}
    hits = {}
    with open(output, 'w') as bedfile:
        with open(input, 'rU') as buscoinput:
            for line in buscoinput:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                if cols[1] == 'Complete' or cols[1] == 'Duplicated':
                    contig = cols[2]              
                    start = cols[3]
                    end = cols[4]
                    score = cols[5]
                    length = cols[6]
                    ID = contig+':'+start+'-'+end           
                    if cols[1] == 'Complete':
                        if not cols[0] in hits:
                            hits[cols[0]] = (ID,score,contig,start,end,length)
                    if ploidy > 1:
                        if cols[1] == 'Duplicated':
                            if not cols[0] in hits:
                                hits[cols[0]] = (ID,score,contig,start,end,length)
                                dupBUSCO2gff(cols[0], os.path.dirname(input), ID)
                            else:
                                oldscore = float(hits.get(cols[0])[1])
                                if float(score) > oldscore:
                                    hits[cols[0]] = (ID,score,contig,start,end,length)
                                    dupBUSCO2gff(cols[0], os.path.dirname(input), ID)
            for k,v in natsorted(hits.items()):
                #validate locations for bedfile, move 100 bp in each direction for bedfile
                start = int(v[3]) - 100
                if start < 1: #negative no good
                    start = 1
                end = int(v[4]) + 100
                if end > ContigSizes.get(contig): #check it doesn't go past contig length
                    end = ContigSizes.get(contig)
                bedfile.write('%s\t%i\t%i\t%s\n' % (contig,start,end,k))         
                busco_complete[k] = v[0]
    return busco_complete

def RepeatBlast(input, cpus, evalue, DataBase, tmpdir, output, diamond=True):
    #run blastp against repeats
    blast_tmp = os.path.join(tmpdir, 'repeats.xml')
    if diamond:
        blastdb = os.path.join(DataBase,'repeats.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input, '--threads', str(cpus), '--out', blast_tmp, '--db', blastdb, '--evalue', str(evalue), '--max-target-seqs', '1', '--outfmt', '5']
    else:
        blastdb = os.path.join(DataBase,'REPEATS')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    runSubprocess(cmd, '.', log)
    #parse results   
    with open(output, 'w') as out:
        with open(blast_tmp, 'rU') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                qlen = qresult.seq_len
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    length = 0
                    for i in range(0,len(hits[0].hsps)):
                        length += hits[0].hsps[i].aln_span
                    pident = hits[0].hsps[0].ident_num / float(length)
                    out.write("%s\t%s\t%f\t%s\n" % (ID, hits[0].id, pident, hits[0].hsps[0].evalue))

def eggnog2dict(annotations):
    #load in annotation dictionary
    EggNog = {}
    with open(annotations, 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    return EggNog
    
def number_present(s):
    return any(i.isdigit() for i in s)
    
def capfirst(x):
    return x[0].upper() + x[1:]
    
def item2index(inputList, item):
    #return the index of an item in the input list
    item_index = None
    for x in inputList:
        if item in x:
            item_index = inputList.index(x)
    return item_index

def getEggNogHeaders(input):
    IDi, DBi, OGi, Genei, COGi, Desci = (None,)*6
    with open(input, 'rU') as infile:
        for line in infile:
            line = line.replace('\n', '')
            if line.startswith('#query_name'): #this is HEADER
                headerCols = line.split('\t')
                IDi = item2index(headerCols, 'query_name')
                Genei = item2index(headerCols, 'predicted_gene_name')
                DBi = item2index(headerCols, 'Annotation_tax_scope')
                OGi = item2index(headerCols, 'OGs')
                COGi = item2index(headerCols, 'COG cat')
                Desci = item2index(headerCols, 'eggNOG annot')
                break
    return IDi, DBi, OGi, Genei, COGi, Desci
    
def parseEggNoggMapper(input, output):
    Definitions = {}
    #indexes from header file
    IDi, DBi, OGi, Genei, COGi, Desci = getEggNogHeaders(input)
    #take annotations file from eggnog-mapper and create annotations
    with open(output, 'w') as out:
        with open(input, 'rU') as infile:
            for line in infile:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                ID = cols[IDi]
                DB = cols[DBi].split('[')[0]
                OGs = cols[OGi].split(',')
                NOG = ''
                for x in OGs:
                    if DB in x:
                        NOG = 'ENOG41'+ x.split('@')[0]
                Gene = ''
                if cols[Genei] != '':
                    if not '_' in cols[Genei] and not '.' in cols[Genei] and number_present(cols[Genei]):
                        Gene = cols[Genei]
                Description = cols[Desci]
                if not ID.endswith('-T1'):
                    ID = ID+'-T1'
                if NOG == '':
                    continue
                if not NOG in Definitions:
                    Definitions[NOG] = Description
                out.write("%s\tnote\tEggNog:%s\n" % (ID, NOG))
                if cols[COGi] != '':
                    out.write("%s\tnote\tCOG:%s\n" % (ID, cols[COGi].replace(' ','')))
                if Gene != '':
                    product = Gene.lower()+'p'
                    product = capfirst(product)                  
                    out.write("%s\tname\t%s\n" % (ID.split('-T1')[0], Gene))
                    out.write("%s\tproduct\t%s\n" % (ID, product))
                    if Description != '':
                        out.write("%s\tnote\t%s\n" % (ID, Description))
    return Definitions

def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def fasta2chunks(input, chunks, tmpdir, output):
    #split the input fasta file into 20 chunks to process
    with open(input, 'rU') as seqs:
        SeqCount = countfasta(input)
        SeqRecords = SeqIO.parse(seqs, 'fasta')
        chunks = SeqCount / int(chunks)
        #divide into chunks, store in tmp file
        folder = os.path.join(tmpdir, output)
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            shutil.rmtree(folder)
            os.makedirs(folder)
        for i, batch in enumerate(batch_iterator(SeqRecords, chunks)) :
            filename = "chunk_%i.fa" % (i+1)
            tmpout = os.path.join(folder, filename)
            handle = open(tmpout, "w")
            count = SeqIO.write(batch, handle, "fasta")
            handle.close()

def signalP(input, tmpdir, output):
    #split input file into chunks, 20 should mean < 200 proteins per chunk
    fasta2chunks(input, 40, tmpdir, 'signalp_tmp')
    for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
        if file.startswith('chunk'):
            file = os.path.join(tmpdir, 'signalp_tmp', file)
            tmp_out = file.replace('.fa', '.signalp.out')
            cmd = ['signalp', '-t', 'euk', '-f', 'short', file]
            runSubprocess2(cmd, '.', log, tmp_out)
    #now concatenate all outputs
    if os.path.isfile(output):
        os.remove(output)            
    with open(output, 'a') as finalout:
        for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
            if file.endswith('.signalp.out'):
                file = os.path.join(tmpdir, 'signalp_tmp', file)
                with open(file) as infile:
                    finalout.write(infile.read())
    #cleanup tmp directory
    shutil.rmtree(os.path.join(tmpdir, 'signalp_tmp'))

def parseSignalP(sigP, secretome_annot):
    sigpDict = {}
    with open(sigP, 'rU') as results:
        for line in results:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            col = line.split(' ') #not tab delimited
            col = filter(None, col) #clean up empty spaces
            if col[9] == 'Y': #then there is signal peptide
                ID = col[0]
                if not ID.endswith('-T1'):
                    ID = ID + '-T1'
                end = int(col[2]) - 1
                sigpDict[ID] = end
    with open(secretome_annot, 'w') as secout:
         for k,v in natsorted(sigpDict.items()):   
            secout.write("%s\tnote\tSECRETED:SignalP(1-%s)\n" % (k, v))

def parsePhobiusSignalP(phobius, sigP, membrane_annot, secretome_annot):
    #give directory of annotate_misc, first get phobius results
    '''
    This is what phobius results look like
    ID  TM  SP  Prediction
    VE00_00001  0   0   o
    VE00_00002  2   0   i198-219o283-301i
    VE00_00003  0   0   o
    VE00_00004  0   Y   n8-18c23/24o
    VE00_00005  12  0   i49-69o89-107i119-138o144-167i179-200o212-234i280-299o319-341i348-366o378-398i410-430o442-465i
    '''
    pSecDict = {}
    pTMDict = {}
    sigpDict = {}
    #parsing short format phobius
    with open(phobius, 'rU') as input1:
        for line in input1:
            line = line.replace('\n', '')
            if line.startswith('ID\t'):
                continue
            cols = line.split('\t')
            geneID = cols[0]
            if not geneID.endswith('-T1'):
                geneID = geneID + '-T1'
            if int(cols[1]) > 0: #then found TM domain
                annot = cols[3]
                if '/' in annot:
                    annotation = annot.split('/')[-1]
                if not geneID in pTMDict:
                    pTMDict[geneID] = 'TransMembrane:'+cols[1]+' ('+annot+')'
            if cols[2] == 'Y': #then sig pep discovered
                location = cols[3].split('/')[0]
                clevage = location.split('c')[-1]
                if not geneID in pSecDict:
                    pSecDict[geneID] = clevage
    if sigP: #will be passed FALSE if signalP data missing
        #parse signalp output and turn into annotation file
        with open(sigP, 'rU') as results:
            for line in results:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                col = line.split(' ') #not tab delimited
                col = filter(None, col) #clean up empty spaces
                if col[9] == 'Y': #then there is signal peptide
                    ID = col[0]
                    if not ID.endswith('-T1'):
                        ID = ID + '-T1'
                    end = int(col[2]) - 1
                    #save as secreted only if also found in phobius
                    if ID in pSecDict:
                        sigpDict[ID] = end
    else:
        sigpDict = pSecDict
    #write annotation files
    with open(membrane_annot, 'w') as memout:
        for k,v in natsorted(pTMDict.items()):
            memout.write("%s\tnote\t%s\n" % (k, v))
    with open(secretome_annot, 'w') as secout:
         for k,v in natsorted(sigpDict.items()):   
            secout.write("%s\tnote\tSECRETED:SignalP(1-%s)\n" % (k, v))
                
def RepeatModelMask(input, cpus, tmpdir, output, debug):
    log.info("Loading sequences and soft-masking genome")
    outdir = os.path.join(tmpdir, 'RepeatModeler')
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    #lets run RepeatModeler here to get repeat library
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    log.info("Soft-masking: building RepeatModeler database")
    with open(debug, 'a') as debug_log:
        subprocess.call(['BuildDatabase', '-name', 'Repeats', input], cwd=outdir, stdout = debug_log, stderr=debug_log)
    log.info("Soft-masking: generating repeat library using RepeatModeler")
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatModeler', '-e', 'ncbi', '-database', 'Repeats', '-pa', str(cpus)], cwd=outdir, stdout = debug_log, stderr=debug_log)
    #find name of folder
    for i in os.listdir(outdir):
        if i.startswith('RM_'):
            RP_folder = i
    library = os.path.join(tmpdir, 'repeatmodeler.lib.fa')
    library = os.path.abspath(library)
    try:
        os.rename(os.path.join(outdir, RP_folder, 'consensi.fa.classified'), library)
    except OSError:
        pass
    #now soft-mask the genome for gene predictors
    outdir2 = os.path.join(tmpdir, 'RepeatMasker')
    if os.path.isdir(outdir2):
        shutil.rmtree(outdir2)
    os.makedirs(outdir2)
    if not os.path.isfile(library):
        log.info("Soft-masking: running RepeatMasker with default library (RepeatModeler found 0 models)")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-gff','-species', 'fungi','-pa', str(cpus), '-xsmall', '-dir','.', input], cwd=outdir2, stdout=debug_log, stderr = debug_log)
    else:
        log.info("Soft-masking: running RepeatMasker with custom library")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-gff','-lib', library, '-pa', str(cpus), '-xsmall', '-dir', '.', input], cwd=outdir2, stdout=debug_log, stderr = debug_log)
    for file in os.listdir(outdir2):
        if file.endswith('.masked'):
            shutil.copyfile(os.path.join(outdir2, file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            cmd = ['rmOutToGFF3.pl', file]
            runSubprocess2(cmd, outdir2, log, rm_gff3)

def RepeatMask(input, library, cpus, tmpdir, output, debug):
    FNULL = open(os.devnull, 'w')
    outdir = os.path.join(tmpdir, 'RepeatMasker')
    #now soft-mask the genome for gene predictors
    log.info("Soft-masking: running RepeatMasker with custom library")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatMasker', '-e', 'ncbi', '-lib', os.path.abspath(library), '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stderr = debug_log, stdout=debug_log, cwd = tmpdir)
    for file in os.listdir(outdir):
        if file.endswith('.masked'):
            os.rename(os.path.join(outdir, file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            cmd = ['rmOutToGFF3.pl', file]
            runSubprocess2(cmd, outdir, log, rm_gff3)

def RepeatMaskSpecies(input, species, cpus, tmpdir, output, debug):
    FNULL = open(os.devnull, 'w')
    outdir = os.path.join(tmpdir, 'RepeatMasker')
    #now soft-mask the genome for gene predictors
    log.info("Soft-masking: running RepeatMasker using %s species" % species)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatMasker', '-e', 'ncbi', '-species', species, '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stderr = debug_log, stdout=debug_log, cwd = tmpdir)
    for file in os.listdir(outdir):
        if file.endswith('.masked'):
            os.rename(os.path.join(outdir, file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            cmd = ['rmOutToGFF3.pl', file]
            runSubprocess2(cmd, outdir, log, rm_gff3)

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())
   
def CheckAugustusSpecies(input):
    #get the possible species from augustus
    augustus_list = []
    for i in os.listdir(os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], 'species')):
        if not i.startswith('.'):
            augustus_list.append(i)
    augustus_list = set(augustus_list)
    if input in augustus_list:
        return True
    else:
        return False

def SortRenameHeaders(input, output):
    #sort records and write temp file
    with open(output, 'w') as out:
        with open(input, 'rU') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = 'scaffold_' + str(counter)
                counter +=1
            SeqIO.write(records, out, 'fasta')

def RunGeneMarkES(input, cpus, tmpdir, output, fungus):
    #make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    contigs = os.path.abspath(input)
    log.info("Running GeneMark-ES on assembly")
    if fungus:
        cmd = ['gmes_petap.pl', '--ES', '--fungus', '--soft_mask', '5000', '--cores', str(cpus), '--sequence', contigs]
    else:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--cores', str(cpus), '--sequence', contigs]
    runSubprocess3(cmd, outdir, log)
    #rename results and grab mod file
    try:
        os.rename(os.path.join(outdir,'output','gmhmm.mod'), os.path.join(tmpdir, 'gmhmm.mod'))
    except OSError:
        log.error("GeneMark-ES failed, please check logfiles.")
        sys.exit(1)
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as out:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' gene models')
        
def RunGeneMark(input, mod, cpus, tmpdir, output, fungus):
    #make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    contigs = os.path.abspath(input)
    mod = os.path.abspath(mod)
    log.info("Running GeneMark-ES on assembly")
    if fungus:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--ini_mod', mod, '--fungus', '--cores', str(cpus), '--sequence', contigs]
    else:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--ini_mod', mod, '--cores', str(cpus), '--sequence', contigs]
    runSubprocess3(cmd, outdir, log)
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as out:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' gene models')

def MemoryCheck():
    import psutil
    mem = psutil.virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)

def systemOS():
    if sys.platform == 'darwin':
        system_os = 'MacOSX '+ platform.mac_ver()[0]
    elif sys.platform == 'linux':
        linux_version = platform.linux_distribution()
        system_os = linux_version[0]+ ' '+linux_version[1]
    else:
        system_os = sys.platform
    return system_os

def SystemInfo():
    system_os = systemOS()
    python_vers = str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])   
    log.info("OS: %s, %i cores, ~ %i GB RAM. Python: %s" % (system_os, multiprocessing.cpu_count(), MemoryCheck(), python_vers))    

def runtRNAscan(input, tmpdir, output):
    tRNAout = os.path.join(tmpdir, 'tRNAscan.out')
    tRNAlenOut = os.path.join(tmpdir, 'tRNAscan.len-filtered.out')
    if os.path.isfile(tRNAout): #tRNAscan can't overwrite file, so check first
        os.remove(tRNAout)
    cmd = ['tRNAscan-SE', '-o', tRNAout, input]
    runSubprocess(cmd, '.', log)
    #enforce NCBI length rules
    with open(tRNAlenOut, 'w') as lenOut:
        with open(tRNAout, 'rU') as infile:
            for line in infile:
                if line.startswith('Sequence') or line.startswith('Name') or line.startswith('--------'):
                    lenOut.write('%s' % line)
                else:
                    cols = line.split('\t')
                    start = cols[2]
                    end = cols[3]
                    if int(start) < int(end):
                        length = abs(int(end) - int(start))
                    else:
                        length = abs(int(start) - int(end))
                    if length < 50 or length > 150:
                        continue
                    else:
                        lenOut.write('%s' % line)

    # now convert to GFF3
    trna2gff = os.path.join(UTIL, 'trnascan2gff3.pl')
    with open(output, 'w') as out:
        subprocess.call(['perl', trna2gff, '--input', tRNAlenOut], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' tRNA gene models')


def runtbl2asn(folder, template, discrepency, organism, isolate, strain, parameters, version):
    '''
    function to run NCBI tbl2asn
    '''
    # get funannotate version
    fun_version = get_version()
    # input should be a folder
    if not os.path.isdir(folder):
        log.error("tbl2asn error: %s is not a directory, exiting" % folder)
        sys.exit(1)
    # based on organism, isolate, strain, construct meta info for -j flag
    if not organism:
        log.error("tbl2asn error: organism not specified")
        sys.exit(1)       
    meta = "[organism=" + organism + "]"
    if isolate:
        isolate_meta = "[isolate=" + isolate + "]"
        meta = meta + " " + isolate_meta
    if strain:
        strain_meta = "[strain=" + strain + "]"
        meta = meta + " " + strain_meta
    cmd = ['tbl2asn', '-y', '"Annotated using '+fun_version+'"', '-N', str(version), '-p', folder, '-t', template, '-M', 'n', '-Z', discrepency, '-j', '"'+meta+'"', '-V', 'b', '-c', 'fx', '-T', '-a', 'r10u']
    # check for custom parameters
    if parameters:
        params = parameters.split(' ')
        cmd = cmd + params
    runSubprocess(cmd, '.', log)
    return ' '.join(cmd)


def gb2smurf(input, prot_out, smurf_out):
    with open(smurf_out, 'w') as smurf:
        with open(prot_out, 'w') as proteins:
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        name = re.sub('[^0-9]','', record.name)
                        if f.type == "CDS":
                            proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0].rstrip('*')))
                            locus_tag = f.qualifiers.get("locus_tag", ["No ID"])[0]
                            product_name = f.qualifiers.get("product", ["No Description"])[0]
                            mystart = f.location.start
                            myend = f.location.end
                            strand = f.location.strand
                            if strand == 1:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(mystart), int(myend), product_name))
                            else:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(myend), int(mystart), product_name))

def GAGprotClean(input, output):
    '''
    gag.py v1 had headers like:
    >>evm.model.Contig100.1 protein
    gag.py v2 has headers like:
    >protein|evm.model.scaffold_1.169 ID=evm.model.scaffold_1.169|Parent=evm.TU.scaffold_1.169|Name=EVM%20prediction%20scaffold_1.169
    '''
    with open(output, 'w') as outfile:
        with open(input, 'ru') as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                if rec.id.startswith('protein|'):
                    ID = rec.id.replace('protein|', '').split(' ')[0]
                else:
                    ID = rec.id.split(' ')[0]
                rec.id = ID
                rec.name = ''
                rec.description = ''
                SeqIO.write(rec, outfile, 'fasta')
                          
def OldRemoveBadModels(proteins, gff, length, repeats, BlastResults, tmpdir, output):
    #first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
    cmd = ['bedtools', 'intersect', '-f', '0.9', '-a', gff, '-b', repeats]
    runSubprocess2(cmd, '.', log, repeat_temp)
    #now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    remove = []
    reason = {}
    #parse the results from bedtools and add to remove list
    with open(repeat_temp, 'rU') as input:
        for line in input:
            if "\tgene\t" in line:
                ninth = line.split('ID=')[-1]
                ID = ninth.split(";")[0]
                remove.append(ID)
                if not ID in reason:
                    reason[ID] = 'remove_reason=repeat_overlap;'
    #parse the results from BlastP search of transposons
    with open(BlastResults, 'rU') as input:
        for line in input:
            col = line.split('\t')
            remove.append(col[0])
            if not col[0] in reason:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_match;'
            else:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_overalap|repeat_match;'
 
    #I'm only seeing these models with GAG protein translations, so maybe that is a problem? skip enforcing start with M
    with open(proteins, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            ID = rec.id.replace('evm.model.', 'evm.TU.')
            if len(Seq) < int(length):
                remove.append(ID)
                if not ID in reason:     
                    reason[ID] = 'remove_reason=seq_too_short;'
            if 'XX' in Seq:
                remove.append(ID)
                if not rec.id in reason:
                    reason[ID] = 'remove_reason=model_span_gap;'
    remove = [w.replace('evm.TU.','') for w in remove]
    remove = [w.replace('evm.model.','') for w in remove]
    remove = set(remove)
    if len(remove) > 0:
        remove_match = re.compile(r'\b\evm.(.*?:%s)[\.;]\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as out2:
                with open(gff, 'rU') as GFF:
                    for line in GFF:
                        if '\tstart_codon\t' in line:
                            continue
                        if '\tstop_codon\t' in line:
                            continue
                        matchLine = remove_match.search(line)
                        if not matchLine:
                            line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                            out.write(line)           
                        else:
                            #print matchLine.group()
                            #print line
                            if "\tgene\t" in line:
                                bad_ninth = line.split('ID=')[-1]
                                bad_ID = bad_ninth.split(";")[0]                
                                bad_reason = reason.get(bad_ID)
                                if bad_reason:
                                    line = line.replace('\n', ';'+bad_reason+'\n')
                                    #print bad_reason
                                else:
                                    log.debug("%s was removed in removeBadModels function for unknown reason, please check manually" % bad_ID)
                                    line = line.replace('\n', ';remove_reason=unknown;\n')
                                    #print 'uknown'
                            out2.write(line)
    else: #if nothing to remove, just print out GFF
        with open(output, 'w') as out:
            with open(gff, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                    out.write(line)

def RemoveBadModels(proteins, gff, length, repeats, BlastResults, tmpdir, output):
    #first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
    cmd = ['bedtools', 'intersect', '-f', '0.9', '-a', gff, '-b', repeats]
    runSubprocess2(cmd, '.', log, repeat_temp)
    #now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    reason = {}
    tooShort = 0
    repeat = 0
    gapspan = 0
    #parse the results from bedtools and add to remove list
    with open(repeat_temp, 'rU') as input:
        for line in input:
            if "\tgene\t" in line:
                ninth = line.split('ID=')[-1]
                ID = ninth.split(";")[0]
                if not ID in reason:
                    reason[ID] = 'remove_reason=repeat_overlap;'
                    repeat += 1
    #parse the results from BlastP search of transposons
    with open(BlastResults, 'rU') as input:
        for line in input:
            col = line.split('\t')
            if not col[0] in reason:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_match;'
                repeat += 1
            else:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_overlap|repeat_match;'
    #Look for models that are too short
    with open(proteins, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            ID = rec.id.replace('evm.model.', 'evm.TU.')
            if len(Seq) < int(length):
                if not ID in reason:     
                    reason[ID] = 'remove_reason=seq_too_short;'
                    tooShort += 1
            if 'XX' in Seq:
                if not rec.id in reason:
                    reason[ID] = 'remove_reason=model_span_gap;'
                    gapspan += 1
    #now read the EVM gene models in Blocks so you can parse gene ID
    numTotal = len(reason)
    if numTotal > 0:
        log.info("Found {:,} gene models to remove: {:,} too short; {:,} span gaps; {:,} transposable elements".format(numTotal,tooShort,gapspan,repeat))
        with open(output, 'w') as out:
            with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as out2:
                with open(gff, 'rU') as GFF:
                    for gene_model in readBlocks(GFF, '\n'):
                        if len(gene_model) > 1:
                            if gene_model[0].startswith('\n'):
                                ID = gene_model[1].split('ID=')[-1].split(';')[0]
                            else:
                                ID = gene_model[0].split('ID=')[-1].split(';')[0]
                            if ID in reason:
                                out2.write('#%s removed; %s\n' % (ID, reason.get(ID)))
                                for line in gene_model:
                                    if not line.startswith('\n'):
                                        out2.write('%s' % (line))
                            else:
                                for line in gene_model:
                                    line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                                    out.write('%s' % (line))

                               
def CleantRNAtbl(GFF, TBL, output):
    #clean up genbank tbl file from gag output
    #try to read through GFF file, make dictionary of tRNA genes and products
    TRNA = {}
    matches = []
    with open(GFF, 'rU') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            line = line.replace('\n', '')
            scaffold, source, feature, start, end, score, orientation, phase, info = line.split('\t')
            if feature == 'tRNA':
                ID = info.split(';')[0].replace('ID=', '')
                ID = ID.replace('-T1', '')
                product = info.split('product=')[-1]
                TRNA[ID] = product
                matches.append(product)
    matches = set(matches)
    tRNAmatch = re.compile(r'\t\t\tproduct\t%s\n' % '|'.join(matches))
    with open(output, 'w') as out:
        with open(TBL, 'rU') as input:
            for line in input:
                if line.startswith('\t\t\tlocus_tag\t'):
                    out.write(line)
                    geneID = line.split('locus_tag\t')[-1].replace('\n', '')
                    if geneID in TRNA:
                        CurrentProduct = TRNA.get(geneID)
                        if 'tRNA-Xxx' == CurrentProduct:
                            out.write("\t\t\tpseudo\n")       
                elif line.startswith("\t\t\tproduct\ttRNA-Xxx"):
                    out.write(line)
                    out.write("\t\t\tpseudo\n")
                    input.next()
                    input.next()
                elif tRNAmatch.search(line):
                    out.write(line)
                    input.next()
                    input.next()
                else: #otherwise just write line
                    out.write(line)

def getFailedProductNames(input, GeneDict):
    #input is NCBI tbl2asn discrepency report, parse to get suspect product names
    failed = {}
    with open(input, 'rU') as discrep:
        for block in readBlocks(discrep, 'DiscRep_'):
            if 'DiscRep_SUB:SUSPECT_PRODUCT_NAMES::' in block[0]:
                reason = []
                for item in block:
                    if item.startswith('DiscRep_SUB:'):
                        bad = item.split('::')[-1].rstrip()
                        if 'features' in bad.lower():
                            bad = bad.split('features ')[-1]
                        reason.append(bad)
                    elif item.startswith('genome:'):
                        gene = item.split('\t')[-1].strip()
                        if gene.startswith('DiscRep'):
                            continue
                        if gene in GeneDict:
                            hit = GeneDict.get(gene)
                            if not hit[0] in failed:
                                failed[hit[0]] = (hit[1], gene, reason)
    return failed


def ParseErrorReport(input, Errsummary, val, Discrep, output, keep_stops):
    errors = []
    gapErrors = []
    remove = []
    with open(Errsummary) as summary:
        for line in summary:
            if 'ERROR' in line:
                if 'SEQ_DESCR.OrganismIsUndefinedSpecies' in line or 'SEQ_DESCR.BadOrgMod' in line or 'SEQ_FEAT.MissingTrnaAA' in line or 'SEQ_INST.TerminalNs' in line: #there are probably other errors you are unaware of....
                    pass
                elif 'SEQ_FEAT.NoStop' in line:
                    if keep_stops:
                        pass
                    else:
                        err = line.split(" ")[-1].rstrip()
                        errors.append(err)
                elif 'SEQ_FEAT.FeatureBeginsOrEndsInGap' in line:
                    err = line.split(" ")[-1].rstrip()
                    gapErrors.append(err)
                else:
                    err = line.split(" ")[-1].rstrip()
                    errors.append(err)
    #parse the discrepency report and look for overlapping genes, so far, all have been tRNA's in introns, so just get those for now.
    with open(Discrep, 'rU') as discrep:
        #process discrepency report into blocks, then look for block headers where overlapping genes are, remove only tRNA models right now
        for block in readBlocks(discrep, 'DiscRep_'):
            if 'DiscRep_ALL:OVERLAPPING_GENES::' in block[0] or 'DiscRep_SUB:RNA_CDS_OVERLAP::' in block[0]:
                for item in block:
                    if item.startswith('genome:tRNA'):
                        gene = item.split('\t')[-1].replace('\n', '')
                        if gene.startswith('DiscRep'):
                            continue
                        tRNA = gene + '_tRNA'
                        exon = gene + '_exon'
                        remove.append(gene)
                        remove.append(tRNA)
                        remove.append(exon)
            if 'DiscRep_ALL:FIND_OVERLAPPED_GENES::' in block[0]:
                for item in block:
                    gene = item.split('\t')[-1].replace('\n', '')
                    if gene.startswith('DiscRep'):
                        continue
                    tRNA = gene + '_tRNA'
                    exon = gene + '_exon'
                    remove.append(gene)
                    remove.append(tRNA)
                    remove.append(exon)

    if len(errors) < 1 and len(remove) < 1: #there are no errors, then just remove stop/start codons and move on
        with open(output, 'w') as out:
            with open(input, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue                 
                    out.write(line)
    else:
        with open(val) as validate:
            for line in validate:
                if any(x in line for x in errors):
                    mRNA = line.split("ncbi|")[-1].replace(']', '').rstrip()
                    gene = mRNA.replace('evm.model', 'evm.TU')
                    exon = mRNA + '.exon'
                    mRNA = mRNA + ';'
                    remove.append(mRNA)
                    remove.append(gene)
                    remove.append(exon)
                #this is only picking up tRNAs right now, which "probably" is all that it needs to.....but u never know
                if any(x in line for x in gapErrors):
                    cols = line.split(' ')
                    if 'Gene:' in cols:
                        gene = line.split('Gene: ')[-1]
                        gene = gene.split(' ')[0]
                        tRNA = gene + '_tRNA'
                        exon = gene + '_exon'
                        remove.append(gene)
                        remove.append(tRNA)
                        remove.append(exon)
        #make sure no empty strings
        remove = list(filter(None, remove))
        remove = set(remove)
        remove_match = re.compile(r'\b(?:%s)+\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(input, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue   
                    if not remove_match.search(line):
                        if '\tgene\t' in line:
                            line = line.replace('Name=;', '')
                        out.write(line)
                        
def ParseAntiSmash(input, tmpdir, output, annotations):
    log.info("Now parsing antiSMASH results, finding SM clusters")
    global bbDomains, bbSubType, BackBone
    BackBone = {}; SMCOGs = {}; bbSubType = {}; bbDomains = {}; smProducts = {}
    backboneCount = 0; clusterCount = 0; cogCount = 0
    #parse antismash genbank to get clusters in bed format and slice the record for each cluster prediction
    with open(output, 'w') as antibed:
        with open(input, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == "source":
                        record_start = f.location.start
                        record_end = f.location.end
                    if f.type == "cluster":
                        clusterCount += 1
                        chr = record.id
                        start = f.location.start
                        end = f.location.end
                        clusternum = f.qualifiers.get("note")[0].replace("Cluster number: ", "")
                        antibed.write("%s\t%s\t%s\tCluster_%s\t0\t+\n" % (chr, start, end, clusternum))
                    Domains = []
                    if f.type == "CDS":
                        ID = f.qualifiers.get('locus_tag')[0]                    
                        if f.qualifiers.get('sec_met'):            
                            for k, v in f.qualifiers.items():
                                if k == 'sec_met':
                                    for i in v:
                                        if i.startswith('Type:'):
                                            type = i.replace('Type: ', '')
                                            backboneCount += 1
                                            BackBone[ID] = type
                                        if i.startswith('NRPS/PKS subtype:'):
                                            subtype = i.replace('NRPS/PKS subtype: ', '')
                                            bbSubType[ID] = subtype
                                        if i.startswith('NRPS/PKS Domain:'):
                                            doms = i.replace('NRPS/PKS Domain: ', '')
                                            doms = doms.split('. ')[0]
                                            Domains.append(doms)
                                bbDomains[ID] = Domains
                        for k,v in f.qualifiers.items():
                            if k == 'note':
                                for i in v:
                                    if i.startswith('smCOG:'):
                                        COG = i.replace('smCOG: ', '')
                                        COG = COG.split(' (')[0]
                                        SMCOGs[ID] = COG
                                        cogCount += 1
                                    elif not i.startswith('smCOG tree'):
                                        notes = i
                                        smProducts[ID] = notes
                            
    log.info("Found %i clusters, %i biosynthetic enyzmes, and %i smCOGs predicted by antiSMASH" % (clusterCount, backboneCount, cogCount))
    #now generate the annotations to add to genome
    with open(annotations, 'w') as out:
        #add product annotations - use bbSubType --> BackBone
        for k, v in natsorted(BackBone.items()):
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            if k in bbSubType:
                hit = bbSubType.get(k)
                if hit == 'NRPS':
                    hit = 'Nonribosomal Peptide Synthase (NRPS)'
                if hit == 'Type I Iterative PKS':
                    hit = 'Type I Iterative Polyketide synthase (PKS)'
            else:
                hit = v
            if hit == 'terpene':
                hit = 'terpene cyclase'
            elif hit == 'other':
                hit = 'putative secondary metabolism biosynthetic enzyme'
            elif hit == 'indole':
                hit = 'aromatic prenyltransferase (DMATS family)'
            elif hit == 'alkaloid' or hit == 'lignan' or hit == 'saccharide' or hit == 'polyketide':
                hit = 'putative ' + hit + ' biosynthetic cluster'
            elif hit == 'putative':
                hit = 'putative uncategorized biosynthetic cluster'
            elif '-' in hit:
                hit = 'putative '+ hit + ' biosynthetic cluster'
            if hit != 'none':
                out.write("%s\tproduct\t%s\n" % (ID, hit))          
        #add annots from smProducts
        for k, v in smProducts.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            if v != 'none' and not 'BLAST' in v:
                sys.stdout.write("%s\tproduct\t%s\n" % (ID, v))               
        #add smCOGs into note section
        for k, v in SMCOGs.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            if v != 'none':
                out.write("%s\tnote\t%s\n" % (ID, v))
              
def GetClusterGenes(input, GFF, output, annotations):
    global dictClusters
    #pull out genes in clusters from GFF3, load into dictionary
    cmd = ['bedtools', 'intersect','-wo', '-a', input, '-b', GFF]
    runSubprocess2(cmd, '.', log, output)
    dictClusters = {}
    with open(output, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            if cols[8] != 'gene':
                continue
            gene = cols[14].replace('ID=', '')
            if ';' in gene:
                gene = gene.replace(';', '')
            ID = cols[3]
            if ID not in dictClusters:
                dictClusters[ID] = [gene]
            else:
                dictClusters[ID].append(gene)
    with open(annotations, 'w') as annotout: 
        for k, v in dictClusters.items():
            for i in v:
                if not i.endswith('-T1'):
                    ID = i + ('-T1')
                else:
                    ID = i
                annotout.write("%s\tnote\tantiSMASH:%s\n" % (ID, k))

def splitFASTA(input, outputdir):
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    with open(input, 'rU') as InputFasta:
        SeqRecords = SeqIO.parse(InputFasta, 'fasta')
        for record in SeqRecords:
            name = str(record.id)
            outputfile = os.path.join(outputdir, name+'.fa')
            with open(outputfile, 'w') as output:
                SeqIO.write(record, output, 'fasta')

def genomeStats(input):
    from Bio.SeqUtils import GC
    lengths = []
    GeeCee = []
    Genes = 0
    tRNA = 0
    Prots = 0
    locus_tag = ''
    organism = None
    isolate = None
    strain = None
    uniqueIso = None
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            lengths.append(len(record.seq))
            GeeCee.append(str(record.seq))
            organism = record.annotations['organism'].replace(' Unclassified.', '')
            for f in record.features:
                if f.type == "source":
                    isolate = f.qualifiers.get("isolate", [None])[0]
                    strain = f.qualifiers.get("strain", [None])[0]
                if f.type == "CDS":
                    Prots += 1
                if f.type == "gene":
                    Genes += 1
                    if Genes == 1:
                        locus_tag = f.qualifiers.get("locus_tag")[0].split('_')[0]
                if f.type == "tRNA":
                    tRNA += 1
    if strain:
        log.info("working on %s %s" % (organism, strain))
        uniqueIso = strain
    elif isolate:
        log.info("working on %s %s" % (organism, isolate))
        uniqueIso = isolate
    else:
        log.info("working on %s" % organism)
    GenomeSize = sum(lengths)
    LargestContig = max(lengths)
    ContigNum = len(lengths)
    AvgContig = int(round(GenomeSize / ContigNum))
    pctGC = round(GC("".join(GeeCee)), 2)
    
    #now get N50
    lengths.sort()
    nlist = []
    for x in lengths:
        nlist += [x]*x
    if len(nlist) % 2 == 0:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos-1]) / 2)
    else:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
    #return values in a list
    return [organism, uniqueIso, locus_tag, "{0:,}".format(GenomeSize)+' bp', "{0:,}".format(LargestContig)+' bp', "{0:,}".format(AvgContig)+' bp', "{0:,}".format(ContigNum), "{0:,}".format(N50)+' bp', "{:.2f}".format(pctGC)+'%', "{0:,}".format(Genes), "{0:,}".format(Prots), "{0:,}".format(tRNA)]

def MEROPS2dict(input):
    dict = {}
    with open(input, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                cols = line.split(' ')
                ID = cols[0].replace('>', '')
                family = cols[1].replace('\n', '')
                dict[ID] = family
    return dict

def getEggNogfromNote(input):
    dict = {}
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k,v in f.qualifiers.items():
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    if not ID in dict:
                                        dict[ID] = hit
    return dict
                                      
def getStatsfromNote(input, word, Database):
    dict = {}
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k,v in f.qualifiers.items():
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if hit.startswith('MER'): #change to family name
                                        hit = meropsDict.get(hit)
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict
    
def getSMBackbones(input):
    dict = {'NRPS': 0, 'PKS': 0, 'Hybrid': 0}
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    product = f.qualifiers['product'][0]
                    if not product == 'hypothetical protein':
                        ID = f.qualifiers['locus_tag'][0]
                        if product == "Hybrid PKS-NRPS":
                            dict['Hybrid'] += 1
                        if product == "Nonribosomal Peptide Synthase (NRPS)":
                            dict['NRPS'] += 1
                        if 'Polyketide synthase (PKS)' in product:
                            dict['PKS'] += 1
    return dict              

def parseGOterms(input, folder, genome):
    with open(os.path.join(folder, 'associations.txt'), 'a') as assoc:
        with open(os.path.join(folder, genome+'.txt'), 'w') as terms:
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == 'CDS':
                            try:
                                ID = f.qualifiers['locus_tag'][0]
                            except KeyError:
                                log.debug("%s has no locus_tag, skipping")
                                continue
                            GOS = []
                            for k,v in f.qualifiers.items():
                                if k == 'note':
                                    notes = v[0].split('; ')
                                    for i in notes:
                                        if i.startswith('GO'):
                                            go_term = i.split(' ')[1]
                                            GOS.append(go_term)
                            if GOS:
                                assoc.write("%s\t%s\n" % (ID, ";".join(GOS)))
                                terms.write("%s\n" % ID)       

def getStatsfromDbxref(input, word):
    dict = {}
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k,v in f.qualifiers.items():
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict

def getGBKannotation(input, Database):
    '''
    Function will loop through GBK file pulling out funannotate functional annotation
    and returning a list of dictionaries for each annotation class
    '''
    #convert merops on the fly, need database
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    SMs = {'NRPS': 0, 'PKS': 0, 'Hybrid': 0}
    pfams = {}
    iprs = {}
    nogs = {}
    cogs = {}
    merops = {}
    cazys = {}
    secreted = {}
    membrane = {}
    buscos = {}
    secmet = {}
    phibase = {}  # added by NW
    with open(input, 'rU') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    product = f.qualifiers['product'][0]
                    if product == "Hybrid PKS-NRPS":
                        SMs['Hybrid'] += 1
                    if product == "Nonribosomal Peptide Synthase (NRPS)":
                        SMs['NRPS'] += 1
                    if 'Polyketide synthase (PKS)' in product:
                        SMs['PKS'] += 1
                    for k,v in f.qualifiers.items():
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith('PFAM:'):
                                    hit = i.replace('PFAM:', '')
                                    if not hit in pfams:
                                        pfams[hit] = [ID]
                                    else:
                                        pfams[hit].append(ID)
                                elif i.startswith('InterPro:'):
                                    hit = i.replace('InterPro:', '')
                                    if not hit in iprs:
                                        iprs[hit] = [ID]
                                    else:
                                        iprs[hit].append(ID)
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:              
                                if i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    if not ID in nogs:
                                        nogs[ID] = hit
                                elif i.startswith('BUSCO:'):
                                    hit = i.replace('BUSCO:', '')
                                    if not hit in buscos:
                                        buscos[hit] = [ID]
                                    else:
                                        buscos[hit].append(ID)
                                elif i.startswith('MEROPS:'): #change to family name
                                    hit = i.replace('MEROPS:', '')
                                    hit = meropsDict.get(hit)
                                    if not hit in merops:
                                        merops[hit] = [ID]
                                    else:
                                        merops[hit].append(ID)
                                elif i.startswith('CAZy:'):
                                    hit = i.replace('CAZy:', '')
                                    if not hit in cazys:
                                        cazys[hit] = [ID]
                                    else:
                                        cazys[hit].append(ID)
                                elif i.startswith('COG:'):
                                    hit = i.replace('COG:', '')
                                    hits = hit.split(',')
                                    for x in hits:
                                        if not x in cogs:
                                            cogs[x] = [ID]
                                        else:
                                            cogs[x].append(ID)                                
                                elif i.startswith('SECRETED:'):
                                    hit = i.replace('SECRETED:', '')
                                    if not hit in secreted:
                                        secreted[hit] = [ID]
                                    else:
                                        secreted[hit].append(ID)                                 
                                elif i.startswith('TransMembrane:'):
                                    hit = i.replace('TransMembrane:', '')
                                    if not hit in membrane:
                                        membrane[hit] = [ID]
                                    else:
                                        membrane[hit].append(ID)
                                elif i.startswith('antiSMASH:'):
                                    hit = i.replace('antiSMASH:', '')
                                    if not hit in secmet:
                                        secmet[hit] = [ID]
                                    else:
                                        secmet[hit].append(ID)
                                elif i.startswith('PHI-Base:'):  # added by NW
                                    hit = i.replace('PHI-Base:', '')
                                    hit = hit.replace('#', ' ')
                                    if not hit in phibase:
                                        phibase[hit] = [ID]
                                    else:
                                        phibase[hit].append(ID)
    return [pfams, iprs, nogs, buscos, merops, cazys, cogs, secreted, membrane, secmet, SMs, phibase]
                                    
def annotationtable(input, Database, output):
    '''
    Function will create a tsv annotation table from GenBank file
    trying to capture all annotation in a parsable tsv file or 
    something that could be imported into excel
    '''
    # convert merops on the fly, need database
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    # input should be fully annotation GBK file from funannotate
    with open(output, 'w') as outfile:
        header = ['GeneID', 'Feature', 'Contig', 'Start', 'Stop', 'Strand', 'Name', 'Product', 'BUSCO', 'PFAM',
                  'InterPro', 'EggNog', 'COG', 'GO Terms', 'Secreted', 'Membrane', 'Protease', 'CAZyme', 'PHI',
                  'Notes', 'Translation']  # added 'PHI'
        outfile.write('%s\n' % '\t'.join(header))
        for record in SeqIO.parse(input, 'genbank'):
            Contig = record.id
            for f in record.features:
                if f.type == 'tRNA':
                    ID = f.qualifiers['locus_tag'][0]
                    Start = f.location.nofuzzy_start
                    End = f.location.nofuzzy_end
                    strand = f.location.strand
                    if strand == 1:
                        Strand = '+'
                    elif strand == -1:
                        Strand = '-'
                    Product = f.qualifiers['product'][0]
                    result = [ID,'tRNA',Contig,str(Start),str(End),Strand,'',Product,'','','','','','','','','','','','','']  # added ,'' for PHI column
                    outfile.write('%s\n' % '\t'.join(result))
                if f.type == 'CDS':
                    ID = f.qualifiers['locus_tag'][0]
                    Start = f.location.nofuzzy_start
                    End = f.location.nofuzzy_end
                    strand = f.location.strand
                    if strand == 1:
                        Strand = '+'
                    elif strand == -1:
                        Strand = '-'
                    Product = f.qualifiers['product'][0]
                    try:
                        Name = f.qualifiers['gene'][0]
                    except KeyError:
                        Name = ''
                    try:
                        Translation = f.qualifiers['translation'][0]
                    except KeyError:
                        Translation = ''
                    pfams = []
                    iprs = []
                    GOS = []
                    nogs = []
                    cogs = []
                    merops = []
                    cazys = []
                    secreted = []
                    membrane = []
                    phi = []  # added by NW
                    therest = []
                    buscos = []
                    for k, v in f.qualifiers.items():
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith('PFAM:'):
                                    hit = i.replace('PFAM:', '')
                                    pfams.append(hit)
                                elif i.startswith('InterPro:'):
                                    hit = i.replace('InterPro:', '')
                                    iprs.append(hit)
                        elif k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('GO'):
                                    go_term = i.split(' ')[1]
                                    GOS.append(go_term)               
                                elif i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    nogs.append(hit)
                                elif i.startswith('BUSCO:'):
                                    hit = i.replace('BUSCO:', '')
                                    buscos.append(hit)
                                elif i.startswith('MEROPS:'): #change to family name
                                    hit = i.replace('MEROPS:', '')
                                    if hit in meropsDict:
                                        hit = meropsDict.get(hit)
                                        merops.append(hit)
                                    else:
                                        log.error("MEROPS database inconsistency: %s not found" % hit)
                                elif i.startswith('CAZy:'):
                                    hit = i.replace('CAZy:', '')
                                    cazys.append(hit)
                                elif i.startswith('COG:'):
                                    hit = i.replace('COG:', '')
                                    hits = hit.split(',')
                                    for x in hits:
                                        desc = x + ':'+ COGS.get(x)
                                        cogs.append(desc)                                
                                elif i.startswith('SECRETED:'):
                                    hit = i.replace('SECRETED:', '')
                                    secreted.append(hit)                                   
                                elif i.startswith('TransMembrane:'):
                                    hit = i.replace('TransMembrane:', '')
                                    membrane.append(hit)
                                elif i.startswith('PHI-Base:'):
                                    hit = i.replace('PHI-Base:', '')
                                    hit = hit.replace('#', ' ')
                                    phi.append(hit)  # added by NW
                                else: #capture everything else
                                    hit = i
                                    therest.append(hit)
                    result = [ID, 'CDS', Contig, str(Start), str(End), Strand, Name, Product, ';'.join(buscos),
                              ';'.join(pfams), ';'.join(iprs), ';'.join(nogs), ';'.join(cogs), ';'.join(GOS),
                              ';'.join(secreted), ';'.join(membrane), ';'.join(merops), ';'.join(cazys), ';'.join(phi),
                              ';'.join(therest), Translation]  # added PHI column results
                    outfile.write('%s\n' % '\t'.join(result))


def ncbiCheckErrors(error, validation, genename, fixOut):
    ncbi_error = 0
    actual_error = 0
    with open(error, 'rU') as errors:
        for line in errors:
            line = line.strip()
            if 'ERROR' in line:
                num = line.split(' ')[0]
                ncbi_error = ncbi_error + int(num)
    #if errors in summary, then parse validation report, only get errors with gene names
    if ncbi_error > 0:
        #see if we can get the gene models that need to be fixed
        needFixing = {}
        with open(validation, 'rU') as validationFile:
            for line in validationFile:
                line = line.strip()
                if line.startswith('ERROR') and genename+'_' in line:
                    actual_error += 1
                    parts = line.split(' ')
                    for x in parts:
                        if genename+'_' in x:
                            ID = x.split('|')[-1]
                    if '-' in ID:
                        ID = ID.split('-')[0]
                    reason = line.split(' FEATURE:')[0]
                    reason = reason.split('] ')[-1]
                    if not ID in needFixing:
                        needFixing[ID] = reason
        if actual_error > 0:
            log.info("There are %i gene models that need to be fixed." % actual_error)
            print('-------------------------------------------------------')
            with open(fixOut, 'w') as fix:
                fix.write('#GeneID\tError Message\n')
                for k,v in natsorted(needFixing.items()):
                    fix.write('%s\t%s\n' % (k,v))
                    print('%s\t%s' % (k,v))
    return actual_error

                  
def convert2counts(input):
    import pandas as pd
    Counts = []
    for i in range(0, len(input)):
        dict = {}
        for k, v in input[i].items():
            dict[k] = len(v)
        Counts.append(dict)
    df = pd.DataFrame(Counts)
    df.fillna(0, inplace=True)  # fill in zeros for missing data
    return df

def gb2proteinortho(input, folder, name):
    history = []
    gffOut = os.path.join(folder, name+'.gff')
    FastaOut = os.path.join(folder, name+'.faa')
    Transcripts = os.path.join(folder, name+'.transcripts.fa')
    with open(gffOut, 'w') as gff:
        with open(FastaOut, 'w') as fasta:
            with open(Transcripts, 'w') as transcripts:
                with open(input, 'rU') as input:
                    SeqRecords = SeqIO.parse(input, 'genbank')
                    for record in SeqRecords:
                        for f in record.features:
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
                            if f.type == 'CDS':
                                try:
                                    locusID = f.qualifiers['locus_tag'][0]
                                except KeyError:
                                    log.debug("%s has no locus_tag, skipping")
                                    continue
                                try:  #saw in a genome downloaded from Genbank that a few models don't have protID?  
                                    protID = f.qualifiers['protein_id'][0]
                                except KeyError:
                                    protID = 'ncbi:'+locusID+'-T1'                       
                                start = f.location.nofuzzy_start
                                end = f.location.nofuzzy_end
                                strand = f.location.strand
                                if strand == 1:
                                    strand = '+'
                                elif strand == -1:
                                    strand = '-'
                                translation = f.qualifiers['translation'][0].rstrip('*')
                                product = f.qualifiers['product'][0]
                                chr = record.id
                                if '.' in chr:
                                    chr = chr.split('.')[0]
                                if not protID in history:
                                    history.append(protID)
                                    gff.write("%s\tNCBI\tCDS\t%s\t%s\t.\t%s\t.\tID=%s;Alias=%s;Product=%s;\n" % (chr, start, end, strand, locusID, protID, product))
                                    fasta.write(">%s\n%s\n" % (locusID, translation))

def drawStackedBar(panda, type, labels, ymax, output, colors=False):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from stackedBarGraph import StackedBarGrapher as StackedBarGrapher
    #stackedbargraph from summary data
    SBG = StackedBarGrapher()
    #labels
    d_labels = panda.index.values
    #y-ticks
    ticks = np.linspace(0,ymax,6)
    ticks = list(ticks)
    nums = [ int(x) for x in ticks ]
    vals = [ str(x) for x in nums ]
    yticks = [nums,vals]
    #colors
    if not colors:
        color_palette = sns.hls_palette(len(panda.columns), l=.4, s=.8).as_hex()
        color_palette = [ str(x).upper() for x in color_palette ]
    else:
        color_palette = colors
    #set up plot
    sns.set_style('darkgrid')
    sns.set_context('paper')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    YLabel = "Number of "+type
    SBG.stackedBarPlot(ax,panda,color_palette,xLabels=panda.index.values,endGaps=True,gap=0.25,xlabel="Genomes",ylabel=YLabel,yTicks=yticks) 
    plt.title(type+" summary")
    #get the legend
    legends = [] 
    i = 0 
    for column in panda.columns: 
        legends.append(mpatches.Patch(color=color_palette[i], label=panda.columns.values[i]+ ": " + labels.get(panda.columns.values[i]))) 
        i+=1 
    lgd = ax.legend(handles=legends, fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.ylim([0,ymax]) 
    #set the font size - i wish I knew how to do this proportionately.....but setting to something reasonable.
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    #setup the plot
    fig.subplots_adjust(bottom=0.4)
    fig.savefig(output, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig) 

def drawHeatmap(df, color, output, labelsize, annotate):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #get size of table
    width = len(df.columns) / 2
    height = len(df.index) / 4
    fig, ax = plt.subplots(figsize=(width,height))
    cbar_ax = fig.add_axes(shrink=0.4)
    if annotate:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, fmt="d", annot_kws={"size": 4}, annot=True)
    else:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, annot=False)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    for item in ax.get_yticklabels():
        item.set_fontsize(int(labelsize))
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)

def donutplot(df, LongName, output, colors=False):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import seaborn as sns
    # create data
    longnames=[]
    for x in df.columns.tolist():
        if x in LongName:
            longnames.append(LongName.get(x))
        else:
            longnames.append(x)
    names = df.columns.tolist() 
    data = df.values.tolist()
    species = df.index.values
    #get size of table
    categories = len(df.columns)
    total = len(df.index)
    Rows = total // 2
    Rows += total % 2
    Position = range(1,total+1)
    #get colors figured out
    if not colors:
        color_palette = pref_colors
    else:
        color_palette = colors
    #draw figure
    if len(species) < 3:
        fig = plt.figure(1,figsize=(8,4))
    else:
        fig = plt.figure(1,figsize=(8,8))
    for k in range(total):
        ax = fig.add_subplot(Rows,2,Position[k])
        # Create a circle for the center of the plot
        my_circle=plt.Circle( (0,0), 0.7, color='white')
        plt.pie(data[0], labels=names, colors=color_palette)
        p=plt.gcf()
        p.gca().add_artist(my_circle)
        plt.title(species[k])
    patches = [ mpatches.Patch(color=color_palette[i], label="{:s}".format(longnames[i]) ) for i in range(len(longnames)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1,0.5), bbox_transform=fig.transFigure, loc="center left", ncol=1)
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)

def drawbarplot(df, output):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #num = len(df.columns) + 1
    sns.set(style="darkgrid")
    fig = plt.figure()
    #colors
    colorplot = sns.husl_palette(len(df), l=.5).as_hex()
    #colorplot = sns.hls_palette(len(df), l=.4, s=.8).as_hex()
    colorplot = [ str(x).upper() for x in colorplot ]
    ax = sns.barplot(data=df, palette=colorplot)
    plt.xlabel('Genomes')
    plt.ylabel('Secreted Proteins')
    plt.xticks(rotation=90)
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig) 
 
def distance2mds(df, distance, type, output):
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from sklearn.metrics.pairwise import pairwise_distances
        from sklearn.manifold import MDS
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #run distance metric on matrix and then plot using NMDS
    num = len(df.index)
    data = np.array(df).astype(int)
    bc_dm = pairwise_distances(data, metric=distance)
    mds = MDS(n_components=2, metric=False, max_iter=999, dissimilarity='precomputed', n_init=10, verbose=0)
    result = mds.fit(bc_dm)
    coords = result.embedding_
    stress = 'stress=' + '{0:.4f}'.format(result.stress_)
    #get axis information and make square plus some padding
    xcoords = abs(maxabs(coords[:,0])) + 0.1
    ycoords = abs(maxabs(coords[:,1])) + 0.1
    #setup plot
    fig = plt.figure()
    #colors
    colorplot = sns.husl_palette(len(df), l=.5).as_hex()
    colorplot = [ str(x).upper() for x in colorplot ]
    for i in range(0,num):
        plt.plot(coords[i,0], coords[i,1], 'o', markersize=14, color=colorplot[i], label=df.index.values[i])
    plt.xlabel('NMDS axis 1')
    plt.ylabel('NMDS axis 2')
    plt.ylim(-ycoords,ycoords)
    plt.xlim(-xcoords,xcoords)
    '''
    if num < 13: #if number too large, don't plot
    '''
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title('NMDS analysis of '+type+' domains')
    plt.annotate(stress, xy=(1,0), xycoords='axes fraction', fontsize=12, ha='right', va='bottom')
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)

def singletons(poff, name):
    with open(poff, 'rU') as input:
        count = 0
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('#'):
                header = line
                species = header.split('\t')[3:]
                i = species.index(name.replace(' ', '_')) + 3
                continue
            col = line.split('\t')
            if col[0] == '1' and col[i] != '*':
                count += 1
        return count

def orthologs(poff, name):
    with open(poff, 'rU') as input:
        count = 0
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('#'):
                header = line
                species = header.split('\t')[3:]
                i = species.index(name.replace(' ', '_')) + 3
                continue
            col = line.split('\t')
            if col[0] != '1' and col[i] != '*':
                count += 1
        return count

def iprxml2dict(xmlfile, terms):
    import xml.etree.cElementTree as cElementTree
    iprDict = {}
    for event, elem in cElementTree.iterparse(xmlfile):
        if elem.tag == 'interpro':
            ID = elem.attrib['id']
            if ID in terms:
                for x in elem.getchildren():
                    if x.tag == 'name':
                        description = x.text
                iprDict[ID] = description
                elem.clear()
            else:
                elem.clear()
    return iprDict


def pfam2dict(file):
    pfamDict = {}
    with open(file, 'rU') as input:
        for line in input:
            if line.startswith('PF'): #just check to be sure
                line = line.replace('\n', '')
                cols = line.split('\t')
                ID = cols[0]
                desc = cols[4]
                pfamDict[ID] = desc
    return pfamDict

def dictFlip(input):
    # flip the list of dictionaries
    outDict = {}
    for x in input:  
        for k, v in natsorted(x.iteritems()):
            for i in v:
                if i in outDict:
                    outDict[i].append(k)
                else:
                    outDict[i] = [k]
    return outDict

def busco_dictFlip(input):
    # flip the list of dictionaries
    output = []
    for x in input:
        outDict = {}
        for k,v in natsorted(x.iteritems()):
            for i in v:
                if i in outDict:
                    outDict[i].append(k)
                else:
                    outDict[i] = [k]
        output.append(outDict)
    return output


def dictFlipLookup(input, lookup):
    outDict = {}
    for x in input:
        for k,v in natsorted(x.iteritems()):
            #lookup description in another dictionary
            if not lookup.get(k) is None:
                result = k+': '+lookup.get(k)
            else:
                result = k+': No description'
            for i in v:
                if i in outDict:
                    outDict[i].append(str(result))
                else:
                    outDict[i] = [str(result)]
    return outDict

def copyDirectory(src, dest):
    import shutil
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

buscoTree='eukaryota (303)\n\tmetazoa (978)\n\t\tnematoda (982)\n\t\tarthropoda (1066)\n\t\t\tinsecta (1658)\n\t\t\tendopterygota (2442)\n\t\t\thymenoptera (4415)\n\t\t\tdiptera (2799)\n\t\tvertebrata (2586)\n\t\t\tactinopterygii (4584)\n\t\t\ttetrapoda (3950)\n\t\t\taves (4915)\n\t\t\tmammalia (4104)\n\t\teuarchontoglires (6192)\n\t\t\tlaurasiatheria (6253)\n\tfungi (290)\n\t\tdikarya (1312)\n\t\t\tascomycota (1315)\n\t\t\t\tpezizomycotina (3156)\n\t\t\t\t\teurotiomycetes (4046)\n\t\t\t\t\tsordariomycetes (3725)\n\t\t\t\t\tsaccharomycetes (1759)\n\t\t\t\t\t\tsaccharomycetales (1711)\n\t\t\tbasidiomycota (1335)\n\t\tmicrosporidia (518)\n\tembryophyta (1440)\n\tprotists (215)\n\t\talveolata_stramenophiles (234)\n'
        
busco_links = {
'fungiv1': 'http://busco.ezlab.org/v1/files/fungi_buscos.tar.gz',
'fungi': 'http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz',
'microsporidia': 'http://busco.ezlab.org/v2/datasets/microsporidia_odb9.tar.gz',
'dikarya': 'http://busco.ezlab.org/v2/datasets/dikarya_odb9.tar.gz',
'ascomycota': 'http://busco.ezlab.org/v2/datasets/ascomycota_odb9.tar.gz',
'pezizomycotina' :'http://busco.ezlab.org/v2/datasets/pezizomycotina_odb9.tar.gz',
'eurotiomycetes' : 'http://busco.ezlab.org/v2/datasets/eurotiomycetes_odb9.tar.gz',
'sordariomycetes' : 'http://busco.ezlab.org/v2/datasets/sordariomyceta_odb9.tar.gz',
'saccharomycetes' : 'http://busco.ezlab.org/v2/datasets/saccharomyceta_odb9.tar.gz',
'saccharomycetales' : 'http://busco.ezlab.org/v2/datasets/saccharomycetales_odb9.tar.gz',
'basidiomycota' : 'http://busco.ezlab.org/v2/datasets/basidiomycota_odb9.tar.gz',
'eukaryota' : 'http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz',
'protists' : 'http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz',
'alveolata_stramenophiles' : 'http://busco.ezlab.org/v2/datasets/alveolata_stramenophiles_ensembl.tar.gz',
'metazoa' : 'http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz',
'nematoda' : 'http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz',
'arthropoda' : 'http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz',
'insecta' : 'http://busco.ezlab.org/v2/datasets/insecta_odb9.tar.gz',
'endopterygota' : 'http://busco.ezlab.org/v2/datasets/endopterygota_odb9.tar.gz',
'hymenoptera' : 'http://busco.ezlab.org/v2/datasets/hymenoptera_odb9.tar.gz',
'diptera' : 'http://busco.ezlab.org/v2/datasets/diptera_odb9.tar.gz',
'vertebrata' : 'http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz',
'actinopterygii' : 'http://busco.ezlab.org/v2/datasets/actinopterygii_odb9.tar.gz',
'tetrapoda' : 'http://busco.ezlab.org/v2/datasets/tetrapoda_odb9.tar.gz',
'aves' : 'http://busco.ezlab.org/v2/datasets/aves_odb9.tar.gz',
'mammalia' : 'http://busco.ezlab.org/v2/datasets/mammalia_odb9.tar.gz',
'euarchontoglires' : 'http://busco.ezlab.org/v2/datasets/euarchontoglires_odb9.tar.gz',
'lauraiatheria' : 'http://busco.ezlab.org/v2/datasets/laurasiatheria_odb9.tar.gz',
'embryophyta' : 'http://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz'}
   
def download_buscos(name, Database):
    if name in busco_links:
        log.info("Downloading %s busco models" % name)
        address = busco_links.get(name)
        filename = address.split('/')[-1]
        if name == 'fungiv1':
            foldername = 'fungi'
        else:
            foldername = filename.split('.')[0]
        cmd = ['wget', '-c', '--tries=0', '--read-timeout=20', address]
        runSubprocess(cmd, '.', log)
        cmd = ['tar', '-zxf', filename]
        runSubprocess(cmd, '.', log)
        copyDirectory(os.path.abspath(foldername), os.path.join(Database, name))
        shutil.rmtree(foldername)
        os.remove(filename)
    else:
        log.error("%s not a valid BUSCO database" % name)
        validBusco = list(busco_links.keys())
        log.error("Valid BUSCO DBs: %s" % (', '.join(validBusco)))
        sys.exit(1)
    
def fasta2dict(Fasta):
    answer = dict()
    with open(Fasta, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'fasta')
        for record in SeqRecords:
            if record.id in answer:
                print "WARNING - duplicate key!"
            else:
                answer[record.id] = str(record.seq)
    return answer 

def ortho2phylogeny(folder, df, num, dict, cpus, bootstrap, tmpdir, outgroup, sp_file, name, sc_buscos):
    import random, pylab
    from Bio import Phylo
    from Bio.Phylo.Consensus import get_support
    if outgroup:
        #load species fasta ids into dictionary
        OutGroup = {}
        with open(sp_file, 'rU') as sp:
            for rec in SeqIO.parse(sp, 'fasta'):
                OutGroup[rec.id] = rec.seq
    #single copy orthologs are in a dataframe, count and then randomly select
    num_species = len(df.columns)
    species = df.columns.values
    if len(df) == 0:
        log.error("0 single copy BUSCO orthologs found, skipping phylogeny")
        return
    if len(df) < int(num):
        number = len(df)
        log.info("Found %i single copy BUSCO orthologs, will use all to infer phylogeny" % (len(df)))
        subsampled = df
    else:
        number = int(num) 
        log.info("Found %i single copy BUSCO orthologs, will randomly select %i to infer phylogeny" % (len(df), number))
        subsampled = df.sample(n=number)

    if outgroup: #passed a list to extract from parent script
        busco_list = sc_buscos

    #since you checked for BUSCO id across all previously, loop through first set and print BUSCOs to file
    with open(os.path.join(tmpdir, 'phylogeny.buscos.used.txt'), 'w') as busco_out:                
        with open(os.path.join(tmpdir, 'phylogeny.concat.fa'), 'w') as proteinout:
            if outgroup:
                proteinout.write(">%s\n" % name)
                for y in busco_list:
                    proteinout.write("%s" % (OutGroup.get(y)))
                proteinout.write('\n')
            for i in range(0,num_species):
                proteinout.write(">%s\n" % species[i])
                proteins = fasta2dict(os.path.join(folder, species[i]+'.faa'))
                for row in subsampled[species[i]].iteritems():
                    proteinout.write("%s" % proteins.get(row[1]))
                    busco_out.write("%s\t%s\n" % (dict[i].get(row[1]), row[1]))
                proteinout.write('\n')
    cmd = ['mafft', '--quiet', os.path.join(tmpdir,'phylogeny.concat.fa')]
    runSubprocess2(cmd, '.', log, os.path.join(tmpdir,'phylogeny.mafft.fa'))
    cmd = ['trimal', '-in', os.path.join(tmpdir,'phylogeny.mafft.fa'), '-out', os.path.join(tmpdir, 'phylogeny.trimal.phylip'), '-automated1', '-phylip']
    runSubprocess(cmd, '.', log)
    if int(cpus) == 1:
        if not outgroup:
            cmd = ['raxmlHPC-PTHREADS', '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
        else:
            cmd = ['raxmlHPC-PTHREADS', '-f', 'a', '-m', 'PROTGAMMAAUTO', '-o', name, '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
    else:
        if not outgroup:
            cmd = ['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
        else:
            cmd = ['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-o', name, '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
    #run RAxML
    runSubprocess(cmd, tmpdir, log)
    #parse with biopython and draw
    trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.nwk'), 'newick'))
    best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.nwk'), 'newick')
    support_tree = get_support(best, trees)
    Phylo.draw(support_tree, do_show=False)
    pylab.axis('off')
    pylab.savefig(os.path.join(tmpdir, 'RAxML.phylogeny.pdf'), format='pdf', bbox_inches='tight', dpi=1000) 

def getTrainResults(input):  
    with open(input, 'rU') as train:
        for line in train:
            if line.startswith('nucleotide level'):
                line = line.replace(' ', '')
                values1 = line.split('|') #get [1] and [2]
            if line.startswith('exon level'):
                line = line.replace(' ', '') #get [6] and [7]
                values2 = line.split('|')
            if line.startswith('gene level'):
                line = line.replace(' ', '')
                values3 = line.split('|') #get [6] and [7]
        return (values1[1], values1[2], values2[6], values2[7], values3[6], values3[7])

def trainAugustus(AUGUSTUS_BASE, train_species, trainingset, genome, outdir, cpus, optimize):
    RANDOMSPLIT = os.path.join(AUGUSTUS_BASE, 'scripts', 'randomSplit.pl')
    OPTIMIZE = os.path.join(AUGUSTUS_BASE, 'scripts', 'optimize_augustus.pl')
    NEW_SPECIES = os.path.join(AUGUSTUS_BASE, 'scripts', 'new_species.pl')
    aug_cpus = '--cpus='+str(cpus)
    species = '--species='+train_species
    aug_log = os.path.join(outdir, 'logfiles', 'augustus_training.log')
    trainingdir = 'tmp_opt_'+train_species
    with open(aug_log, 'w') as logfile:
        if not CheckAugustusSpecies(train_species):
            subprocess.call(['perl', NEW_SPECIES, species], stdout = logfile, stderr = logfile)
        #run etraining again to only use best models from EVM for training
        subprocess.call(['etraining', species, trainingset], stderr = logfile, stdout = logfile)
        subprocess.call(['perl', RANDOMSPLIT, trainingset, '200']) #split off 200 models for testing purposes
        if os.path.isfile(os.path.join(outdir, 'predict_misc', 'busco.training.gb.train')):
            with open(os.path.join(outdir, 'predict_misc', 'augustus.initial.training.txt'), 'w') as initialtraining:
                subprocess.call(['augustus', species, trainingset+'.test'], stdout=initialtraining)
            train_results = getTrainResults(os.path.join(outdir, 'predict_misc', 'augustus.initial.training.txt'))
            log.info('Initial training: '+'{0:.2%}'.format(float(train_results[4]))+' genes predicted exactly and '+'{0:.2%}'.format(float(train_results[2]))+' of exons predicted exactly')
            if optimize:
                #now run optimization
                subprocess.call(['perl', OPTIMIZE, species, aug_cpus, '--onlytrain='+trainingset+'.train', trainingset+'.test'], stderr = logfile, stdout = logfile)
                #run etraining again
                subprocess.call(['etraining', species, trainingset], stderr = logfile, stdout = logfile)
                with open(os.path.join(outdir, 'predict_misc', 'augustus.final.training.txt'), 'w') as finaltraining:
                    subprocess.call(['augustus', species, trainingset+'.test'], stdout=finaltraining)
                train_results = getTrainResults(os.path.join(outdir, 'predict_misc', 'augustus.final.training.txt'))
                log.info('Optimized training: '+'{0:.2%}'.format(float(train_results[4]))+' genes predicted exactly and '+'{0:.2%}'.format(float(train_results[2]))+' of exons predicted exactly')
                #clean up tmp folder
                shutil.rmtree(trainingdir)
            else:
                if float(train_results[4]) < 0.50:
                    log.info("Accuracy seems low, you can try to improve by passing the --optimize_augustus option.")
        else:
            log.error("AUGUSTUS training failed, check logfiles")
            sys.exit(1)

def checkgoatools(input):
    with open(input, 'rU') as goatools:
        count = -1
        result = False
        headercount = 0
        for line in goatools:
            count += 1
            if line.startswith('GO\tNS'):
                header = line.replace('\n', '')
                headercount = count
            if line.startswith('GO:'):
                result = True
    return (result, headercount)

def translatemRNA(input, output):
    with open(output, 'w') as outfile:
        with open(input, 'rU') as fasta:
            for rec in SeqIO.parse(fasta, 'fasta'):
                rec.seq = rec.seq.translate()
                SeqIO.write(rec, outfile, 'fasta')

def alignMAFFT(input, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['mafft', '--quiet', input], stderr = FNULL, stdout = outfile)

def align2Codon(alignment, transcripts, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['perl', os.path.join(UTIL,'pal2nal.pl'), alignment, transcripts, '-output', 'fasta'], stderr=FNULL, stdout = outfile)
    if getSize(output) < 1:
        os.remove(output)
        log.debug('dNdS Error: pal2nal failed for %s' % alignment)

def counttaxa(input):
    ct = 0
    with open(input, 'rU') as tree:
        line = tree.readline()
        ct = line.count(',')+1
    return ct

def getMatchFileName(pattern, directory):
    result = None
    for f in os.listdir(directory):
        if pattern in f:
            result = os.path.join(directory, f)
    return result

def drawPhyMLtree(fasta, tree):
    FNULL = open(os.devnull, 'w')
    fc = countfasta(fasta)
    #need to convert to phylip format
    base = os.path.basename(fasta).split('.')[0]
    dir = os.path.dirname(fasta)
    tmp1 = os.path.join(dir, base+'.draw2tree.phylip')
    subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
    #draw tree
    subprocess.call(['phyml', '-i', tmp1], stdout = FNULL, stderr = FNULL)
    tmp2 = getMatchFileName(base+'.draw2tree.phylip_phyml_tree', dir)
    #check that num taxa in tree = input
    tc = counttaxa(tmp2)
    if tc != fc: #something failed...
        log.debug('dNdS Error: phyml tree failed for %s' % fasta)
        #retry
        subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
        subprocess.call(['phyml', '-i', tmp1], stdout = FNULL, stderr = FNULL)
    #rename and clean
    os.rename(tmp2, tree)
    SafeRemove(tmp1)
    stats = getMatchFileName(base+'.draw2tree.phylip_phyml_stats', dir)
    SafeRemove(stats)

def simplestTreeEver(fasta, tree):
    with open(tree, 'w') as outfile:
        with open(fasta, 'rU') as input:
            ids = []
            for rec in SeqIO.parse(input, 'fasta'):
                ids.append(rec.id)
            outfile.write('(%s,%s);' % (ids[0], ids[1]))

def rundNdSexhaustive(folder):
    FNULL = open(os.devnull, 'w')
    #setup intermediate files
    tmpdir = os.path.dirname(folder)
    name = os.path.basename(folder)
    transcripts = os.path.join(tmpdir, name+'.transcripts.fa')
    prots = os.path.join(tmpdir, name+'.proteins.fa')
    aln = os.path.join(tmpdir, name+'.aln')
    codon = os.path.join(tmpdir, name+'.codon.aln')
    tree = os.path.join(tmpdir, name+'.tree')
    log = os.path.join(tmpdir, name+'.log')
    finallog = os.path.join(tmpdir, name, name+'.log')
    if not checkannotations(finallog):
        num_seqs = countfasta(transcripts)
        #Translate to protein space
        translatemRNA(transcripts, prots)
        #align protein sequences
        alignMAFFT(prots, aln)
        #convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                #now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            #now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(tree), '--models', 'M0', 'M1', 'M2', 'M7', 'M8', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd = tmpdir, stdout = logfile, stderr = logfile)
    #clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file), os.path.join(tmpdir, name, file))            
              

def rundNdSestimate(folder):
    FNULL = open(os.devnull, 'w')
    #setup intermediate files
    tmpdir = os.path.dirname(folder)
    name = os.path.basename(folder)
    transcripts = os.path.join(tmpdir, name+'.transcripts.fa')
    prots = os.path.join(tmpdir, name+'.proteins.fa')
    aln = os.path.join(tmpdir, name+'.aln')
    codon = os.path.join(tmpdir, name+'.codon.aln')
    tree = os.path.join(tmpdir, name+'.tree')
    log = os.path.join(tmpdir, name+'.log')
    finallog = os.path.join(tmpdir, name, name+'.log')
    if not checkannotations(finallog):
        num_seqs = countfasta(transcripts)
        #Translate to protein space
        translatemRNA(transcripts, prots)
        #align protein sequences
        alignMAFFT(prots, aln)
        #convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                #now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            #now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(tree), '--models', 'M0', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd = tmpdir, stdout = logfile, stderr = logfile)
    #clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file), os.path.join(tmpdir, name, file))            

def get_subdirs(a_dir):
    return [os.path.join(a_dir, name) for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]                       

def get_subdirs2(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]                       

def parsedNdS(folder):
    results = {}
    hits = get_subdirs2(folder)
    for x in hits:
        finallog = os.path.join(folder, x, x+'.log')
        #parse logfile to get omega
        dnds = 'NA'
        m1m2p = 'NA'
        m7m8p = 'NA'
        if os.path.isfile(finallog):
            with open(finallog, 'rU') as input:
                for line in input:
                    line = line.strip()
                    if 'M7' in line and 'M8' in line and '|' in line:
                        m7m8p = line.split('|')[-1].strip()
                        m7m8p = m7m8p.replace('*','')
                        m7m8p = '{0:.5f}'.format(float(m7m8p))  
                    elif 'M1' in line and 'M2' in line and '|' in line:
                        m1m2p = line.split('|')[-1].lstrip()
                        m1m2p = m1m2p.replace('*','')
                        m1m2p = '{0:.5f}'.format(float(m1m2p))
                    elif line.startswith('- Model M0'):
                        nextline = next(input)             
                        dnds = nextline.split('tree: ')[1].rstrip()        
        results[x] =  (dnds, m1m2p, m7m8p)
    return results     


def chunkIt(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0
  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg
  return out

def getBlastDBinfo(input):
    '''
    function to return a tuple of info using blastdbcmd
    tuple: (name, date, #sequences)
    '''
    cmd = ['blastdbcmd', '-info', '-db', input]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr.split('\n')[0]
    results = stdout.split('\n\n')
    results = [x for x in results if x]
    #parse results which are now in list, look for starts with Database and then Date
    Name, Date, NumSeqs = (None,)*3
    for x in results:
        if x.startswith('Database:'):
            hit = x.split('\n\t')
            Name = hit[0].replace('Database: ', '')
            NumSeqs = hit[1].split(' sequences;')[0].replace(',', '')
        if x.startswith('Date:'):
            Date = x.split('\t')[0].replace('Date: ', '')
    return (Name, Date, NumSeqs)


HEADER = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
    <meta name="funannotate comparative genomics output" content="">
    <meta name="Jonathan Palmer" content="">
    <title>Funannotate</title>
    <!-- Bootstrap core CSS -->
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <!-- Custom styles for this template -->
    <link href="css/starter-template.css" rel="stylesheet">
    <script src="js/ie-emulation-modes-warning.js"></script>
  </head>
  <body>
    <nav class="navbar navbar-inverse navbar-fixed-top">
      <div class="container">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
          </button>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li><a href="stats.html">Stats</a></li>
            <li><a href="phylogeny.html">Phylogeny</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterPro</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
            <li><a href="cogs.html">COGs</a></li>
            <li><a href="signalp.html">SignalP</a></li>
            <li><a href="tf.html">TFs</a></li>
            <li><a href="secmet.html">SecMet</a></li>
            <li><a href="go.html">GO</a></li>
            <li><a href="citation.html">Cite</a></li>
          </ul>
        </div><!--/.nav-collapse -->
      </div>
    </nav>
'''
ORTHOLOGS = '''
    <div class="container">
      <div class="table">
        <h2 class="sub-header">Orthologous protein groups</h2>
          <div class="table-responsive">
'''
INDEX = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Funannotate Results</h2>
         <br>
         <p><a href='stats.html'>Genome Summary Stats</a></p>
         <p><a href='phylogeny.html'>Maximum likelihood Phylogeny (RAxML)</a></p>
         <p><a href='merops.html'>MEROPS Protease Stats</a></p>
         <p><a href='cazy.html'>CAZyme carbohydrate activating enzyme Stats</a></p>
         <p><a href='cogs.html'>COGs Stats</a></p>
         <p><a href='signalp.html'>Secreted proteins (SignalP)</a></p>
         <p><a href='interpro.html'>InterProScan Domain Stats</a></p>
         <p><a href='tf.html'>Transcription Factor Summary</a></p>
         <p><a href='secmet.html'>Secondary Metabolism Cluster Summary</a></p>
         <p><a href='pfam.html'>PFAM Domain Stats</a></p>
         <p><a href='go.html'>Gene Ontology Enrichment Analysis</a></p>
         <p><a href='orthologs.html'>Orthologous proteins</a></p>
         <br>
'''
SUMMARY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Genome Summary Stats</h2>
          <div class="table-responsive">
'''
PHYLOGENY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">RAxML Maximum Likelihood Phylogeny</h2>
        <a href='phylogeny/RAxML.phylogeny.pdf'><img src="phylogeny/RAxML.phylogeny.pdf" height="500" /></a></div>
'''
NOPHYLOGENY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Number of species too low to generate phylogeny</h2>
'''
MEROPS = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">MEROPS Protease Families per Genome Results</h2>
        <div class='row'>
        <div class="col-sm-7"><a href='merops/MEROPS.graph.pdf'><img src="merops/MEROPS.graph.pdf" height="350" /></a></div>
        <div class="col-sm-5"><a href='merops/MEROPS.heatmap.pdf'><img src="merops/MEROPS.heatmap.pdf" height="500" /></a></div>
        </div>
        <div class="table-responsive">
'''
INTERPRO = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">InterProScan Domains per Genome Results</h2>
        <div class='row'>
        <a href='interpro/InterProScan.nmds.pdf'><img src="interpro/InterProScan.nmds.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
PFAM = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">PFAM Domains per Genome Results</h2>
        <div class='row'>
        <a href='pfam/PFAM.nmds.pdf'><img src="pfam/PFAM.nmds.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
SIGNALP = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Secreted Proteins per Genome Results</h2>
        <div class='row'>
        <a href='signalp/signalp.pdf'><img src="signalp/signalp.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
TF = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Fungal Transcription Factors per Genome Results</h2>
        <div class='row'>
        <a href='tfs/TF.heatmap.pdf'><img src="tfs/TF.heatmap.pdf" height="800" /></a></div>
        <div class="table-responsive">
'''
SECMET = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Secondary Metabolism Clusters per Genome Results</h2>
        <div class='row'>
        <a href='secmet/SM.graph.pdf'><img src="secmet/SM.graph.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''

CAZY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">CAZyme Families per Genome Results</h2>
        <div class='row'>
        <div class="col-sm-7"><a href='cazy/CAZy.graph.pdf'><img src="cazy/CAZy.graph.pdf" height="350" /></a></div>
        <div class="col-sm-5"><a href='cazy/CAZy.heatmap.pdf'><img src="cazy/CAZy.heatmap.pdf" height="600" /></a></div>
        </div>
        <div class="table-responsive">
'''

COG = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Clusters of Orthologous Groups (COGs) per Genome Results</h2>
        <div class='row'>
        <a href='cogs/COGS.graph.pdf'><img src="cogs/COGS.graph.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''

GO = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">GO ontology enrichment Results</h2>
        <div class='row'>
'''
MISSING = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">These data are missing from annotation.</h2>
'''
CITATION = '''
    <div class="container">
      <div class="starter-template">
        <h3 class="sub-header">If you found Funannotate useful please cite:</h3>
        <p>Palmer JM. 2016. Funannotate: a fungal genome annotation and comparative genomics pipeline. <a href="https://github.com/nextgenusfs/funannotate">https://github.com/nextgenusfs/funannotate</a>.</p>
'''
FOOTER = '''
          </div>  
      </div>

    </div><!-- /.container -->


    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script>window.jQuery || document.write('<script src="js/jquery.min.js"><\/script>')</script>
    <script src="js/bootstrap.min.js"></script>
    <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
    <script src="js/ie10-viewport-bug-workaround.js"></script>
  </body>
</html>

'''
HEADER2 = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="funannotate comparative genomics output" content="">
    <meta name="Jonathan Palmer" content="">
    <title>Funannotate</title>
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <link href="css/starter-template.css" rel="stylesheet">
    <script src="js/ie-emulation-modes-warning.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/t/bs/dt-1.10.11/datatables.min.css"/>
    <script type="text/javascript" src="https://cdn.datatables.net/t/bs/dt-1.10.11/datatables.min.js"></script>
  </head>
  <body>
    <nav class="navbar navbar-inverse navbar-fixed-top">
      <div class="container-fluid">
        <div class="navbar-header">
            <span class="sr-only">Toggle navigation</span>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div class="navbar-header">
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li class="active"><a href="stats.html">Stats</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterProScan</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
            <li><a href="signalp.html">SignalP</a></li>
            <li><a href="go.html">GO ontology</a></li>
            <li><a href="citation.html">Citation</a></li>
            <li class="dropdown">
          <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Genomes <span class="caret"></span></a>
          <ul class="dropdown-menu">
'''
