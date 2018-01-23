# python3

import os, logging, subprocess, re

from subprocess import Popen, PIPE, STDOUT
from multiprocessing.dummy import Pool as ThreadPool

CDBYANKEXT = '.cidx'

logging.basicConfig()

# this could get updated/overwritten with a config file perhaps?
# otherwise assuming apps are in the path
# and if wanted to generalize this where there might be a different file name...
Apps = { 'cdbfasta': 'cdbfasta',
         'cdbyank':  'cdbyank'}


def init_allseqdb(dbfolder,dbpath, dbext, force):
    dbpathidx = dbpath+CDBYANKEXT
    makedb = (force or not os.path.exists(dbpath))
    makeidx = (force or not os.path.exists(dbpathidx))

    # check and see if we need to make the db file if any of the seqfiles is newer
    if not makedb:
        dbtime = os.path.getmtime(dbpath)    
        for file in os.listdir(dbfolder):
            if file.endswith(dbext):
                ftime = os.path.getmtime(os.path.join(dbfolder,file))
                if ftime > dbtime:
                    logging.info("db %s older than filetime %s of %s" % (dbtime,ftime,file))
                    print("Remaking db %s as file %s is more recent" % (dbpath, file))
                    makedb = True
                    break

    if makedb:
        # concat
        with open(dbpath,"w") as seqdb:
            for file in os.listdir(dbfolder):
                if file.endswith(dbext):
                    with open(os.path.join(dbfolder,file)) as input:
                        for line in input:
                            seqdb.write(line)
        makeidx = True
    elif not makeidx:
        dbtime = os.path.getmtime(dbpath)
        dbxtime = os.path.getmtime(dbpathidx)
        if dbxtime < dbtime:
            logging.info("Remaking Index for %s, index is older than dbfile" % (dbpath))
            print("Remaking Index for %s" % (dbpath))
            makeidx = True
                            
    if makeidx:
        subprocess.call([Apps["cdbfasta"],dbpath])

    return dbpathidx


def run_cdbyank (fileset) :
    index = fileset[0]
    outfile = fileset[1]
    names = fileset[2]
    p = subprocess.Popen([Apps["cdbyank"],index,
                          "-o",outfile],stdin=PIPE)
    # only call this once with the complete list of IDs otherwise 
    # the process gets closed
    p.communicate(input=names.encode())


def make_unaln_files (search_dir, best_extension, cutoff, 
                      dbidx, outdir, outext, force, threads=2):
    orthologs = {}
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    if not os.path.exists(dbidx):
        print("No dbidx %s exists for reading" % (dbidx))
        return -1
    
    for file in os.listdir(search_dir):
        if file.endswith("."+best_extension):
            with open(os.path.join(search_dir,file),"r") as fh:
                for line in fh:
                    row = line.strip().split("\t")
                    if float(row[2]) <= float(cutoff):
                        if row[0] in orthologs:
                            orthologs[row[0]].append(row[1])
                        else:
                            orthologs[row[0]] = [ row[1] ]        

    pool = ThreadPool(threads)
    fileset = []
    for orth in orthologs:
        outfile = "%s.%s" % (os.path.join(outdir,orth),outext)
        if force or (not os.path.exists(outfile)):        
            fileset.append( [dbidx, outfile,
                             "\n".join(orthologs[orth]) + "\n"])

    results = pool.map(run_cdbyank, fileset)
    
    # close the pool and wait for the work to finish 
    pool.close() 
    pool.join() 

                
# sequences = read_fasta
def read_fasta (file):
    seq = ""
    fasta_pat = re.compile("^>(\S+)")
    seqs = []
    with open(file, 'r') as f:
        seenheader = 0
        id = ""
        for line in f:
            line = line.strip()
            m = fasta_pat.match(line)

            if m:
                if seenheader:
                    seq = re.sub("\s+","",seq)
                    seqs.append([id,seq])
                    id = ""
                    seq = ""

                id = m.group(1)
                seenheader = 1
            else:
                seq += line
        # fencepost
        seqs.append([id,seq])
    return seqs
