

def read_refseqscan_results(fn):
    """Read RefSeqScan output file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[2]
    return ret


def read_excluded_transcripts(fn):
    """Read excluded transcripts from txt file"""

    ret = dict()
    for line in open(fn):
        line = line.strip()
        if line.startswith('#'):
            continue
        cols = line.split()
        ret[cols[0]] = cols[1]
    return ret


def output_genepred(transcript, outfile):
    """Output transcript in GenePred format"""

    record = [
                transcript.cartid,
                transcript.chrom,
                transcript.strand,
                str(transcript.start),
                str(transcript.end),
                str(min(transcript.coding_start, transcript.coding_end)),
                str(max(transcript.coding_start, transcript.coding_end)+1),
                str(len(transcript.exons)),
                ''.join([str(e.start)+',' for e in transcript.exons]),
                ''.join([str(e.end)+',' for e in transcript.exons]),
                transcript.hgnc_id
            ]

    outfile.write('\t'.join(record) + '\n')


def output_gff3(transcript, outfile):
    """Output transcript in GFF3 format"""

    attr = ';'.join(['ID=' + transcript.cartid, 'HGNCID=' + transcript.hgnc_id])
    outfile.write('\t'.join([transcript.chrom, '.', 'transcript', str(transcript.start + 1), str(transcript.end), '.', transcript.strand, '.', attr]) + '\n')

    # Exons
    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        exon_id = 'EXON' + transcript.cartid[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + exon_id, 'Parent=' + transcript.cartid])
        outfile.write('\t'.join([transcript.chrom, '.', 'exon', str(exon.start + 1), str(exon.end), '.', transcript.strand, '.', attr]) + '\n')

    # CDS
    cds_regs = transcript.cds_regions()
    cdspos = 0
    for i in range(len(cds_regs)):
        cds_reg = cds_regs[i]

        cds_id = 'CDS' + transcript.cartid[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + cds_id, 'Parent=' + transcript.cartid])

        if cdspos % 3 == 0:
            phase = 0
        elif cdspos % 3 == 1:
            phase = 2
        else:
            phase = 1

        outfile.write('\t'.join([transcript.chrom, '.', 'CDS', str(cds_reg[0] + 1), str(cds_reg[1]), '.', transcript.strand, str(phase), attr]) + '\n')
        cdspos += cds_reg[1] - cds_reg[0]
