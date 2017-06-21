#!env/bin/python

from optparse import OptionParser
from cartcreation.main import main
import datetime


# Version
ver = '0.1.1'

# Command line argument parsing
descr = 'CARTCreation v'+ver
parser = OptionParser(usage='CARTCreation/env/bin/cart_creation <options>', version=ver, description=descr)
parser.add_option("-i", dest='input', action='store', help="Input file")
parser.add_option("-d", dest='refsdb', action='store', help="RefSeq transcript database file")
parser.add_option("-s", dest='refss', action='store', help="refseq_scan output file")
parser.add_option("-o", dest='out', action='store', help="Output file prefix")
(options, args) = parser.parse_args()

print '\n'+'='*100
print 'CARTCreation v'+ver+' started: '+str(datetime.datetime.now())+'\n'

main(options)

print '\nFinished: '+str(datetime.datetime.now())
print '='*100+'\n'
