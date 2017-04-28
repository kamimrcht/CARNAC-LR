import sys
import gzip


def get_read(sequencefile,offset):
    sequencefile.seek(offset)
    read = ""
    line = sequencefile.readline()
    if not line: 
        print "cannot read read at offset", offset
        exit(1)
    read += line  # include header
    read += sequencefile.readline()  # include sequence
    if read[0] == '>': return read
    if read[0] != '@': 
        print "read offset", offset, "does not start with @ or >"
        exit(1)
    read += sequencefile.readline() # include header2
    read += sequencefile.readline() # include quality


def index_bank_offsets(bank_file_name):
    print "indexing reads..."
    read_offsets = []
    if "gz" in bank_file_name:
        sequencefile = gzip.open(bank_file_name,"r")
    else: 
        sequencefile = open(bank_file_name,"r")
    line = sequencefile.readline()
    real_bank_names = []  # for fof
    fof = false
    if not line: 
        print "Can't open file", bank_file_name
        exit(1)
    if line[0] != '@' and line[0] != '>':
        if "fa" or "fasta" or "fq" or "fastq" in line[0]:  # file of file
            real_bank_names.append(line[0])
        else:
            print "File", bank_file_name, "not correctly formatted"
            exit(1)
    if fof:
        for bank_name in real_bank_names:
            if "gz" in bank_file_name:
                sequencefile = gzip.open(bank_name,"r")
            else: 
                sequencefile = open(bank_name,"r")
            return_read_offsets(sequencefile, read_offsets)
    else:
        return_read_offsets(sequencefile, read_offsets)



def return_read_offsets(sequencefile, read_offsets):
    linesperread=2  # fasta by default
    if line[0] =='@':
        linesperread=4 # fastq
    
    sequencefile.seek(0)
    while True:
        offset = sequencefile.tell()
        line = sequencefile.readline()
        if not line: break
        read_offsets.append(offset)
        for i in range(linesperread-1):
            line = sequencefile.readline()
    sequencefile.close()
    return read_offsets



 
def convert_clustering_output(read_offsets, clustering_output_file_name, read_bank_queries_file_name):
    print "Converting clustering output to reads"
    #OPEN clustering OUTPUT
    if "gz" in clustering_output_file_name:
        clusterfile = gzip.open(clustering_output_file_name,"r")
    else: 
        clusterfile = open(clustering_output_file_name,"r")
    
    #OPEN BANK OUTPUT
    if "gz" in read_bank_queries_file_name:
        bankfile = gzip.open(read_bank_queries_file_name,"r")
    else: 
        bankfile = open(read_bank_queries_file_name,"r")
        
    #one cluster per line, each read index separated by spaces
    for line in clusterfile.readlines():
        line=line.rstrip()
        query_read_ids = line.split(' ')
        for read in query_read_ids:
            print get_read(bankfile,read_offsets[int(read)])
    clusterfile.close()
    bankfile.close()
 

if len(sys.argv) < 2 :
     print "USAGE"
     print " Used to transform the clustering output into a set of reads from the query. Each cluster will be converted into a separated fasta file."
     print "COMMAND"
     print  sys.argv[0], "<query set file (fasta/q format, each sequence on ONE line, gzipped or not)> <clustering output file>"
     sys.exit(1)


read_bank_queries_file_name = sys.argv[1]
clustering_output_file_name = sys.argv[2]
read_offsets = index_bank_offsets(read_bank_queries_file_name)
convert_clustering_output(read_offsets, clustering_output_file_name, read_bank_queries_file_name)

