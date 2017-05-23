
# Class representing a transcript
class Transcript(object):

    def __init__(self, line, build):
        cols = line.split()
        self.id = cols[0]
        self.hgncid = None
        self.cartid = None

        mapping_str = cols[6] if build == 'GRCh37' else cols[7]
        mapping = mapping_str.split(';')
        self.chrom = mapping[0]
        self.strand = mapping[1]

        if self.strand == '+':
            self.coding_start = int(mapping[3])
            self.coding_end = int(mapping[4]) - 1
        else:
            self.coding_start = int(mapping[4]) - 1
            self.coding_end = int(mapping[3])

        self.exons = []
        for exonstr in mapping[2].split(','):
            self.exons.append(Exon(exonstr))

        if self.strand == '+':
            self.start = self.exons[0].start
            self.end = self.exons[-1].end
        else:
            self.start = self.exons[-1].start
            self.end = self.exons[0].end

    def cds_regions(self):
        ret = []
        coding_5prime = self.coding_start if self.strand == '+' else self.coding_end
        coding_3prime = self.coding_end if self.strand == '+' else self.coding_start
        for i in range(len(self.exons)):
            exon = self.exons[i]
            if exon.end < coding_5prime or exon.start > coding_3prime:
                continue
            cds_start = max(exon.start, coding_5prime)
            cds_end = min(exon.end, coding_3prime)
            ret.append((cds_start, cds_end))
        return ret

    def output_genepred(self, outfile):
        # TODO write this output function
        record = [self.cartid, self.chrom, self.strand, str(self.start), str(self.end), str(self.coding_start), str(self.coding_end), str(len(self.exons))]
        outfile.write('\t'.join(record)+'\n')

    def output_gff3(self, outfile):
        attr = ';'.join(['ID=' + self.cartid, 'HGNCID=' + self.hgncid])
        outfile.write('\t'.join([self.chrom, '.', 'transcript', str(self.start+1), str(self.end+1), '.', self.strand, '.', attr])+'\n')

        # Exons
        for i in range(len(self.exons)):
            exon = self.exons[i]
            exon_id = 'EXON' + self.cartid[-5:] + '.' + str(i + 1)
            attr = ';'.join(['ID=' + exon_id, 'Parent=' + self.cartid])
            outfile.write('\t'.join([self.chrom, '.', 'exon', str(exon.start+1), str(exon.end+1), '.', self.strand, '.', attr])+'\n')

        # CDS
        cds_regs = self.cds_regions()
        cdspos = 0
        for i in range(len(cds_regs)):
            cds_reg = cds_regs[i]

            cds_id = 'CDS' + self.cartid[-5:] + '.' + str(i + 1)
            attr = ';'.join(['ID=' + cds_id, 'Parent=' + self.cartid])

            outfile.write('\t'.join([self.chrom, '.', 'CDS', str(cds_reg[0]+1), str(cds_reg[1]+1), '.', self.strand, str(cdspos % 3), attr]) + '\n')
            cdspos += cds_reg[1] - cds_reg[0] + 1


# Class representing an exon
class Exon(object):

    def __init__(self, s):
        [s, e] = s.split('-')
        self.start = int(s)
        self.end = int(e) - 1
