from utils.transcripts import TranscriptDB


def read_refseqscan_results(fn):
    """Read RefSeqScan output file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[3]
    return ret


def read_excluded_transcripts(fn):
    """Read excluded transcripts from txt file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line.startswith('#'): continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def output_genepred(transcript, outfile):
    """Output transcript in GenePred format"""

    # TODO write this output function
    record = [transcript.cartid, transcript.chrom, transcript.strand, str(transcript.start), str(transcript.end), str(transcript.coding_start),
              str(transcript.coding_end), str(len(transcript.exons))]
    outfile.write('\t'.join(record) + '\n')


def output_gff3(transcript, outfile):
    """Output transcript in GFF3 format"""

    attr = ';'.join(['ID=' + transcript.cartid, 'HGNCID=' + transcript.hgncid])
    outfile.write('\t'.join([transcript.chrom, '.', 'transcript', str(transcript.start + 1), str(transcript.end + 1), '.', transcript.strand, '.', attr]) + '\n')

    # Exons
    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        exon_id = 'EXON' + transcript.cartid[-5:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + exon_id, 'Parent=' + transcript.cartid])
        outfile.write('\t'.join([transcript.chrom, '.', 'exon', str(exon.start + 1), str(exon.end + 1), '.', transcript.strand, '.', attr]) + '\n')

    # CDS
    cds_regs = transcript.cds_regions()
    cdspos = 0
    for i in range(len(cds_regs)):
        cds_reg = cds_regs[i]

        cds_id = 'CDS' + transcript.cartid[-5:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + cds_id, 'Parent=' + transcript.cartid])

        outfile.write('\t'.join(
            [transcript.chrom, '.', 'CDS', str(cds_reg[0] + 1), str(cds_reg[1] + 1), '.', transcript.strand, str(cdspos % 3), attr]) + '\n')
        cdspos += cds_reg[1] - cds_reg[0] + 1


def main(options):
    """Main function"""

    # Read in RefSeq transcript database
    refseq_db = TranscriptDB(options.refsdb)
    refseq_db.read()

    # Read excluded transcripts from file
    excluded = read_excluded_transcripts(options.input[:-3]+'_excluded.txt')

    # Read in refseq_scan output
    refseqscan = read_refseqscan_results(options.refss)

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

        if not refseq_db.contains(nm_id):
            # TODO Handle this case
            continue

        if nm_id in excluded:
            # TODO Handle this case
            continue

        # ....
        t = refseq_db.by_id(nm_id)

        # ...
        t.cartid = cart_id
        t.hgncid = hgnc_id

        # ...
        t.output_genepred(out_genepred)
        t.output_gff3(out_gff)


    # Close output files
    out_genepred.close()
    out_gff.close()