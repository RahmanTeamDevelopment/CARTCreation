
# Class representing a transcript
class Transcript(object):

    # Constructor
    def __init__(self):
        self.id = ''
        self.version = '.'
        self.gene = ''
        self.chrom = ''
        self.strand = ''
        self.exons = []
        self.cds_exons = []
        self.utr5_exons = []
        self.utr3_exons = []
        self.cds = []
        self.seq = ''

    # Read from RefSeq db record
    def read(self, line, build):
        cols = line.split()
        self.id = cols[0]
        self.version = cols[1]
        self.gene = cols[2]
        self.seq = cols[5]

        if build == 'GRCh37': mapping_str = cols[6]
        else: mapping_str = cols[7]
        mapping = mapping_str.split(';')

        self.chrom = mapping[0]
        self.strand = 1 if mapping[1] == '+' else -1

        if self.strand == 1: self.cds = [int(mapping[3]) + 1, int(mapping[4])]
        else: self.cds = [int(mapping[4]), int(mapping[3]) + 1]

        self.exons = []
        for exonstr in mapping[2].split(','):
            [s,e] = exonstr.split('-')
            self.exons.append([int(s)+1, int(e)])
        self.split_exons()

    # Split exons to UTR5, CDS and UTR3
    def split_exons(self):
        if self.strand == 1:
            for e in self.exons:
                if e[1] < self.cds[0]:
                    self.utr5_exons.append(e)

                elif e[0] > self.cds[1]:
                    self.utr3_exons.append(e)

                elif e[0] < self.cds[0] and self.cds[0] <= e[1] <= self.cds[1]:
                    self.utr5_exons.append([e[0], self.cds[0] - 1])
                    self.cds_exons.append([self.cds[0], e[1]])

                elif self.cds[0] <= e[0] <= self.cds[1] and self.cds[1] < e[1]:
                    self.cds_exons.append([e[0], self.cds[1]])
                    self.utr3_exons.append([self.cds[1] + 1, e[1]])

                elif self.cds[0] <= e[0] <= self.cds[1] and self.cds[0] <= e[1] <= self.cds[1]:
                    self.cds_exons.append(e)

                elif e[0] < self.cds[0] and self.cds[1] < e[1]:
                    self.utr5_exons.append([e[0], self.cds[0] - 1])
                    self.cds_exons.append([self.cds[0], self.cds[1]])
                    self.utr3_exons.append([self.cds[1] + 1, e[1]])
        else:
            for e in self.exons:
                if self.cds[0] < e[0]:
                    self.utr5_exons.append(e)

                elif e[1] < self.cds[1]:
                    self.utr3_exons.append(e)

                elif e[0] < self.cds[1] and self.cds[1] <= e[1] <= self.cds[0]:
                    self.utr3_exons.append([e[0], self.cds[1] - 1])
                    self.cds_exons.append([self.cds[1], e[1]])

                elif self.cds[1] <= e[0] <= self.cds[0] and self.cds[0] < e[1]:
                    self.cds_exons.append([e[0], self.cds[0]])
                    self.utr5_exons.append([self.cds[0] + 1, e[1]])

                elif self.cds[1] <= e[0] <= self.cds[0] and self.cds[1] <= e[1] <= self.cds[0]:
                    self.cds_exons.append(e)

                elif e[0] < self.cds[1] and self.cds[0] < e[1]:
                    self.utr3_exons.append([e[0], self.cds[1] - 1])
                    self.cds_exons.append([self.cds[1], self.cds[0]])
                    self.utr5_exons.append([self.cds[0] + 1, e[1]])

    # Return 5' start
    def utr5_start(self):
        if self.strand == 1:
            if len(self.utr5_exons) > 0: return self.utr5_exons[0][0]
            else: return self.cds_exons[0][0]
        else:
            if len(self.utr5_exons) > 0: return self.utr5_exons[0][1]
            else: return self.cds_exons[0][1]

    # Return 3' end
    def utr3_end(self):
        if self.strand == 1:
            if len(self.utr3_exons) > 0: return self.utr3_exons[-1][1]
            else: return self.cds_exons[-1][1]
        else:
            if len(self.utr3_exons) > 0: return self.utr3_exons[-1][0]
            else: return self.cds_exons[-1][0]

    # Return UTR_ends
    def UTR_ends(self):
        return [self.utr5_start(), self.utr3_end()]

    # Return length of UTR5 exonic content
    def utr5_exonic_content_length(self):
        ret = 0
        for e in self.utr5_exons: ret += e[1]-e[0]+1
        return ret

    # Return length of UTR3 exonic content
    def utr3_exonic_content_length(self):
        ret = 0
        for e in self.utr3_exons: ret += e[1]-e[0]+1
        return ret
