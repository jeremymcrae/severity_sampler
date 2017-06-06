"""
Copyright (c) 2017 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from itertools import groupby, count

def as_range(g):
    l = list(g)
    return l[0], l[-1]

def get_severity(cadd, chrom, rates):
    ''' get CADD scores for a specific alt at a specific site
    
    See downloadable CADD files here: http://cadd.gs.washington.edu/download
    
    Args:
        cadd: pysam.TabixFile for quick fetching of CADD scores
        chrom: chromosome
        rates: Weighted Choice object for all sites in a gene
    
    Returns:
        list of 9CADD score at the given site for the given alt allele
    '''
    
    try:
        positions = sorted(set([ x['position'] for x in rates ]))
    except KeyError:
        positions = sorted(set([ x['pos'] for x in rates ]))
    
    scores = {}
    for _, g in groupby(positions, key=lambda n, c=count(): n-next(c)):
        # define the start and end of a continuous set of positions, so we can
        # load cadd scores by iterating through the cadd file efficiently
        start, end = as_range(g)
        for line in cadd.fetch(chrom, start-1, end):
            _, pos, _, alt, _, score = line.split('\t')
            scores[(int(pos), alt)] = float(score)
    
    # match the cadd scores to the order of sites in the rates object
    try:
        return [ scores[(x['position'], x['alt'])] for x in rates ]
    except KeyError:
        return [ scores[(x['pos'], x['alt'])] for x in rates ]
