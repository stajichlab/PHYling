# python3

import os, logging, subprocess, re

from subprocess import Popen, PIPE, STDOUT
from multiprocessing.dummy import Pool as ThreadPool

CDBYANKEXT = '.cidx'
SFETCHEXT  = '.ssi'

logging.basicConfig()

# this could get updated/overwritten with a config file perhaps?
# otherwise assuming apps are in the path
# and if wanted to generalize this where there might be a different file name...
Apps = { 'cdbfasta': 'cdbfasta',
         'cdbyank':  'cdbyank',
         'sfetch' : 'esl-sfetch',
}


def init_allseqdb(dbfolder,dbpath,dbext,
                  index_type='sfetch',force=False):
    dbpathidx = ""
    index_type = index_type.lower()
    if index_type == 'cdbfasta':
        dbpathidx = dbpath+CDBYANKEXT
    elif index_type == 'sfetch':
        dbpathidx = dbpath+SFETCHEXT
    else:
        print("need an expected index type not %s" % (index_type))
        exit(-1)

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
        if index_type == "cdbfasta":
            subprocess.call([Apps["cdbfasta"],dbpath])
        elif index_type == "sfetch":
            subprocess.call([Apps["sfetch"],
                             "--index",
                             dbpath])


def run_sfetch(fileset):
    dbfile = fileset[0]
    outfile = fileset[1]
    names = fileset[2]
    p = subprocess.Popen([Apps["sfetch"],"-o", outfile,
                          "-f",dbfile,"-"],stdin=PIPE)
    # only call this once with the complete list of IDs otherwise
    # the process gets closed
    p.communicate(input=names.encode())

def run_cdbyank (fileset):
    index = fileset[0]+CDBYANKEXT
    outfile = fileset[1]
    names = fileset[2]
    p = subprocess.Popen([Apps["cdbyank"],index,
                          "-o",outfile],stdin=PIPE)
    # only call this once with the complete list of IDs otherwise
    # the process gets closed
    p.communicate(input=names.encode())


def make_unaln_files (search_dir, best_extension, cutoff,
                      dbpath, outdir, outext, force=False,
                      index_type="sfetch",
                      threads=2,multi=False):

    orthologs = {}
    dbidx = ""
    index_type = index_type.lower()
    if index_type == "cdbfasta":
        dbidx = dbpath + CDBYANKEXT
    elif index_type == "sfetch":
        dbidx = dbpath + SFETCHEXT

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
                    HMMname = row.pop(0) # take first col as HMM name

                    if best_extension == "best":
                        if float(row[1]) <= float(cutoff):
                            if HMMname in orthologs:
                                orthologs[HMMname].append(row[0])
                            else:
                                orthologs[HMMname] = [ row[0] ]
                    elif best_extension == "best_multi":
                        if not multi:
                            row = [ row[0] ]

                        for hit in row:
                            hit_dat = hit.split(",")

                            if HMMname in orthologs:
                                orthologs[HMMname].append(hit_dat[0])
                            else:
                                orthologs[HMMname] = [ hit_dat[0] ]

    pool = ThreadPool(threads)
    fileset = []
    for orth in orthologs:
        outfile = "%s.%s" % (os.path.join(outdir,orth),outext)
        if force or (not os.path.exists(outfile)):
            fileset.append( [dbpath, outfile,
                             "\n".join(orthologs[orth]) + "\n"])

    if index_type == "cdbfasta":
        results = pool.map(run_cdbyank, fileset)
    elif index_type == "sfetch":
        results = pool.map(run_sfetch, fileset)

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
