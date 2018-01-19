# python3

import os, logging, subprocess, re

from subprocess import Popen, PIPE, STDOUT

logging.basicConfig()

# this could get updated/overwritten with a config file perhaps?
# otherwise assuming apps are in the path
# and if wanted to generalize this where there might be a different file name...
Apps = { 'cdbfasta': 'cdbfasta',
         'cdbyank':  'cdbyank'}


def init_allseqdb(dbfolder,dbname, dbext, force):
    dbpath = os.path.join(dbfolder,dbname)
    dbpathidx = dbpath+".cidx"
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


def make_unaln_files (search_dir, best_extension, cutoff, 
                      dbidx, outdir, outext, force):
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
                    if row[2] <= cutoff:
                        if row[0] in orthologs:
                            orthologs[row[0]].append(row[1])
                        else:
                            orthologs[row[0]] = [ row[1] ]

    for orth in orthologs:
        print(orth)
        outfile = "%s.%s" % (os.path.join(outdir,orth), outext)
        if force or (not os.path.exists(outfile)):
            p = subprocess.Popen([Apps["cdbyank"],dbidx,
                                  "-o",outfile],stdin=PIPE)
            instr = ""
            for gene in orthologs[orth]:
                instr += "%s\n"%(gene)
            # only call this once with the complete list of IDs otherwise the process gets closed
            p.communicate(input=instr.encode())
                
# sequences = read_fasta
def read_fasta (file):
    seq = ""
    fasta_pat = re.compile(">(\S+)")
    seqs = []
    with open(file, 'r') as f:
        seenheader = 0
        id = ""
        for line in f:
            line = line.strip()
            m = fasta_pat.match(line);
            if m:
                id = m.group(1)
                if seenheader:
                    seq = re.sub("\s+","",seq)
                    seqs.append([id,seq])
                    id = ""
                    seq = ""
                seenheader = 1
            else:
                seq += line
        seqs.append([id,seq])
        
    return seqs
