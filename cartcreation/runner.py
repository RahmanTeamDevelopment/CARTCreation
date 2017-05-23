from core import Transcript
import gzip

# Read transcript database from file
def read_transcript_database(fn):
    ret = dict()
    for line in gzip.open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = line
    return ret

# Read RefSeqScan output file
def read_refseqscan_output(fn):
    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[3]
    return ret

# Check if NM transcript does not have multiple or no mapping
def check_mapping(refseq_db, build, nm_id):
    cols = refseq_db[nm_id].split()
    mapping_str = cols[6] if build == 'GRCh37' else cols[7]
    return mapping_str not in ['multiple_mapping', 'no_mapping']

########################################################################################################################

def run(options):

    # Read in RefSeq transcript database
    refseq_db = read_transcript_database(options.refsdb)

    # Read in refseq_scan output
    refseqscan = read_refseqscan_output(options.refss)

    # Initialize output files
    out_genepred = open(options.out + '.genepred', 'w')
    out_gff = open(options.out + '.gff', 'w')
    out_gff.write('##gff-version 3\n\n')

    # Iterate through input file records
    for line in open(options.input):

        line = line.strip()
        if line == '' or line.startswith('#'): continue
        cols = line.split()
        if len(cols) != 3: continue

        hgnc_id = cols[0]
        cart_id = cols[1]
        nm_id = cols[2]

        if nm_id not in refseq_db:
            # TODO Handle this case
            continue

        # Check mapping
        if not check_mapping(refseq_db, options.build, nm_id):
            # TODO Handle this case
            continue

        # Creat Transcript object
        t = Transcript(refseq_db[nm_id], options.build)
        t.cartid = cart_id
        t.hgncid = hgnc_id

        # Write to output files
        t.output_genepred(out_genepred)
        t.output_gff3(out_gff)

    # Close output files
    out_genepred.close()
    out_gff.close()