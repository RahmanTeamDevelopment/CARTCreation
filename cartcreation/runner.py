
def run(options):

    for line in open(options.input):

        line = line.strip()
        if line == '' or line.startswith('#'): continue
        cols = line.split()
        if len(cols) != 3: continue

        hgnc_id = cols[0]
        cart_id = cols[1]
        nm_id = cols[2]



