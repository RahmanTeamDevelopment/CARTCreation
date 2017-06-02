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

        self.t_forward = Transcript(chrom='12', strand='+', start=1000, end=2000)
        self.t_forward.coding_start = 1100
        self.t_forward.coding_end = 1900
        self.t_forward.hgnc_id = 'HGNC:1234'
        self.t_forward.cartid = 'CART12345'
        self.t_forward.exons = [Exon('1100-1300'), Exon('1500-1600'), Exon('1800-2000')]

        self.t_reverse = Transcript(chrom='9', strand='-', start=5000, end=10000)
        self.t_reverse.coding_start = 9000
        self.t_reverse.coding_end = 5500
        self.t_reverse.hgnc_id = 'HGNC:2468'
        self.t_reverse.cartid = 'CART54321'
        self.t_reverse.exons = [Exon('8000-10000'), Exon('7000-7500'), Exon('6000-6900'), Exon('5000-5800')]


    def tearDown(self):
        self.out.close()
        os.remove(self.outfn)


    def test_output_genepred_forward(self):

        helper.output_genepred(self.t_forward, self.out)

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
            '1901',
            '3',
            '1100,1500,1800,',
            '1300,1600,2000,',
            'HGNC:1234'
        ]


    def test_output_gff3_forward(self):

        helper.output_gff3(self.t_forward, self.out)

        self.out.close()

        correct = [
            ['12', '.', 'transcript', '1001', '2000', '.', '+', '.', 'ID=CART12345;HGNCID=HGNC:1234'],
            ['12', '.',	'exon',	'1101',	'1300',	'.', '+', '.', 'ID=EXON12345.1;Parent=CART12345'],
            ['12', '.', 'exon', '1501', '1600', '.', '+', '.', 'ID=EXON12345.2;Parent=CART12345'],
            ['12', '.', 'exon', '1801', '2000', '.', '+', '.', 'ID=EXON12345.3;Parent=CART12345'],
            ['12', '.',	'CDS', '1101', '1300',	'.', '+', '0', 'ID=CDS12345.1;Parent=CART12345'],
            ['12', '.', 'CDS', '1501', '1600', '.', '+', '1', 'ID=CDS12345.2;Parent=CART12345'],
            ['12', '.', 'CDS', '1801', '1901', '.', '+', '0', 'ID=CDS12345.3;Parent=CART12345']
        ]

        i = 0
        for line in open(self.outfn):
            line = line.strip()
            assert line.split() == correct[i]
            i += 1



    def test_output_genepred_reverse(self):

        helper.output_genepred(self.t_reverse, self.out)

        self.out.close()

        for line in open(self.outfn):
            line = line.strip()
            break

        assert line.split() == [
            'CART54321',
            '9',
            '-',
            '5000',
            '10000',
            '5500',
            '9001',
            '4',
            '5000,6000,7000,8000,',
            '5800,6900,7500,10000,',
            'HGNC:2468'
        ]


    def test_output_gff3_reverse(self):

        helper.output_gff3(self.t_reverse, self.out)

        self.out.close()

        correct = [
            ['9', '.', 'transcript', '5001', '10000', '.', '-', '.', 'ID=CART54321;HGNCID=HGNC:2468'],
            ['9', '.', 'exon', '8001', '10000', '.', '-', '.', 'ID=EXON54321.1;Parent=CART54321'],
            ['9', '.', 'exon', '7001', '7500', '.', '-', '.', 'ID=EXON54321.2;Parent=CART54321'],
            ['9', '.', 'exon', '6001', '6900', '.', '-', '.', 'ID=EXON54321.3;Parent=CART54321'],
            ['9', '.', 'exon', '5001', '5800', '.', '-', '.', 'ID=EXON54321.4;Parent=CART54321'],
            ['9', '.', 'CDS', '8001', '9001', '.', '-', '0', 'ID=CDS54321.1;Parent=CART54321'],
            ['9', '.', 'CDS', '7001', '7500', '.', '-', '1', 'ID=CDS54321.2;Parent=CART54321'],
            ['9', '.', 'CDS', '6001', '6900', '.', '-', '2', 'ID=CDS54321.3;Parent=CART54321'],
            ['9', '.', 'CDS', '5501', '5800', '.', '-', '2', 'ID=CDS54321.4;Parent=CART54321']
        ]

        i = 0
        for line in open(self.outfn):
            line = line.strip()
            assert line.split() == correct[i]
            i += 1
