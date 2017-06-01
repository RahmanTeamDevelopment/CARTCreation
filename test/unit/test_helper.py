"""Unit tests for helper.py"""


from unittest import TestCase
from utils.transcripts import Transcript, Exon
from cartcreation import helper
import uuid
import os


class TestHelper(TestCase):


    def setUp(self):
        self.outfn = str(uuid.uuid4())
        self.out = open(self.outfn, 'w')

        self.t1 = Transcript(chrom='12', strand='+', start=1000, end=2000)
        self.t1.coding_start = 1100
        self.t1.coding_end = 1900
        self.t1.hgnc_id = 'HGNC:1234'
        self.t1.cartid = 'CART12345'
        self.t1.exons = [Exon('1100-1300'),Exon('1500-1600'),Exon('1800-1900')]


    def tearDown(self):
        self.out.close()
        os.remove(self.outfn)


    def test_output_genepred(self):

        helper.output_genepred(self.t1, self.out)

        self.out.close()

        for line in open(self.outfn):
            line = line.strip()
            break

        assert line.split() == [
            'CART12345',
            '12',
            '+',
            '1000',
            '2000',
            '1100',
            '1900',
            '3',
            '1100,1500,1800,',
            '1300,1600,1900,',
            'HGNC:1234'
        ]


    def test_output_gff3(self):
        helper.output_gff3(self.t1, self.out)

        self.out.close()

        correct = [
            ['12', '.', 'transcript', '1001', '2000', '.', '+', '.', 'ID=CART12345;HGNCID=HGNC:1234'],
            ['12', '.',	'exon',	'1101',	'1300',	'.', '+', '.', 'ID=EXON12345.1;Parent=CART12345'],
            ['12', '.', 'exon', '1501', '1600', '.', '+', '.', 'ID=EXON12345.2;Parent=CART12345'],
            ['12', '.', 'exon', '1801', '1900', '.', '+', '.', 'ID=EXON12345.3;Parent=CART12345'],
            ['12', '.',	'CDS', '1101', '1300',	'.', '+', '0', 'ID=CDS12345.1;Parent=CART12345'],
            ['12', '.', 'CDS', '1501', '1600', '.', '+', '1', 'ID=CDS12345.2;Parent=CART12345'],
            ['12', '.', 'CDS', '1801', '1900', '.', '+', '0', 'ID=CDS12345.3;Parent=CART12345']
        ]

        i = 0
        for line in open(self.outfn):
            line = line.strip()
            assert line.split() == correct[i]
            i += 1
