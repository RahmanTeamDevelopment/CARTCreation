from utils.transcripts import TranscriptDB
import sys
import helper


def main(options):
    """Main function"""

    # Read in RefSeq transcript database
    refseq_db = TranscriptDB(options.refsdb)
    refseq_db.read()
    print 'RefSeq transcript database has been read successfully.'

    # Read excluded transcripts from file
    excluded = helper.read_excluded_transcripts(options.refsdb[:-3]+'_excluded.txt')

    # Read in refseq_scan output
    refseqscan = helper.read_refseqscan_results(options.refss)
    print 'RefSeqScan output file has been read successfully.'

    # Initialize output files
    our_list = open(options.out + '.txt', 'w')
    our_list.write('\t'.join(['#HGNC_ID', 'CART_ID', 'NM', 'GenomeDifference', 'MissingReason']) + '\n')
    out_genepred = open(options.out + '.genepred', 'w')
    out_genepred.write('\t'.join(['#NAME', 'CHROM', 'STRAND', 'START', 'END', 'CDS_START', 'CDS_END', 'EXON_COUNT', 'EXON_STARTS', 'EXON_ENDS', 'HGNC_ID']) + '\n')
    out_gff = open(options.out + '.gff', 'w')
    out_gff.write('##gff-version 3\n\n')

    # Iterate through input file records
    sys.stdout.write('\nProcessing input file ...')
    sys.stdout.flush()
    for line in open(options.input):

        line = line.strip()
        if line == '' or line.startswith('#'): continue
        cols = line.split()
        if len(cols) != 3: continue

        hgnc_id = cols[0]
        cart_id = cols[1]
        nm_id = cols[2]

        # Transcripts excluded from the RefSeq database by RefSeqDB
        if nm_id in excluded:
            if excluded[nm_id] in ['missing_cds', 'missing_sequence', 'missing_version', 'missing_hgncid']:
                reason = 'incomplete_refseq_record'
            else:
                reason = excluded[nm_id]
            our_list.write('\t'.join([hgnc_id, cart_id, nm_id, 'missing', reason]) + '\n')
            continue

        # Transcripts missing from the RefSeq db
        if not refseq_db.contains(nm_id):
            our_list.write('\t'.join([hgnc_id, cart_id, nm_id, 'missing', 'missing_from_refseq_database']) + '\n')
            continue

        # Retrieve transcript from database
        t = refseq_db.by_id(nm_id)

        # Checking HGNC ID mismatch
        if hgnc_id != t.hgnc_id:
            our_list.write('\t'.join([hgnc_id, cart_id, nm_id, 'missing', 'hgnc_id_mismatch']) + '\n')
            continue

        # Assign CART ID
        t.cartid = cart_id

        # Output in GenePred and GFF3 formats
        helper.output_genepred(t, out_genepred)
        helper.output_gff3(t, out_gff)

        # Output to list file
        our_list.write('\t'.join([hgnc_id, cart_id, nm_id, refseqscan[nm_id], '.']) + '\n')

    print ' - Done.'
    print '\nOutput files written:'
    print ' - ' + options.out + '.txt'
    print ' - ' + options.out + '.genepred'
    print ' - ' + options.out + '.gff'

    # Close output files
    our_list.close()
    out_genepred.close()
    out_gff.close()