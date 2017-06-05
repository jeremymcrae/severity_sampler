
def get_severity(cadd, chrom, pos, alt):
    ''' get CADD scores for a specific alt at a specific site
    
    See downloadable CADD files here: http://cadd.gs.washington.edu/download
    
    Args:
        cadd: pysam.TabixFile for quick fetching of CADD scores
        chrom: chromosome
        pos: nucleotide position
        alt: alternate allele
    
    Returns:
        CADD score at the given site for the given alt allele
    '''
    
    for line in cadd.fetch(chrom, pos-1, pos):
        _, temp_pos, _, temp_alt, _, scaled = line.split('\t')
        if int(pos) == pos and temp_alt == alt:
            return float(scaled)
