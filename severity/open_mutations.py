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

MISSENSE_CQ = set(["missense_variant", "stop_lost", "inframe_deletion",
    "inframe_insertion", "coding_sequence_variant", "protein_altering_variant"])

LOF_CQ = set(["stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "frameshift_variant", "initiator_codon_variant",
    "start_lost", "conserved_exon_terminus_variant"])

def parse_header(header):
    ''' parse the column positions of required fields from the header
    
    Args:
        header: string for header line
    
    Returns:
        dictionary of positions, indexed by column type
    '''
    header = header.strip().split('\t')
    
    chrom = header.index('chrom')
    pos = header.index('pos')
    ref = header.index('ref')
    alt = header.index('alt')
    symbol = header.index('symbol')
    cq = header.index('consequence')
    
    return {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'symbol': symbol,
        'cq': cq}

def parse(line, indices):
    ''' parse a variant line
    
    Args:
        line: string for variant line
        indices: dictionary of column positions for various required fields
    
    Returns:
        tuple of chrom, pos, ref, alt, symbol and consequence
    '''
    line = line.strip().split('\t')
    
    chrom = line[indices['chrom']]
    pos = int(line[indices['pos']])
    ref = line[indices['ref']]
    alt = line[indices['alt']]
    symbol = line[indices['symbol']]
    cq = line[indices['cq']]
    
    return chrom, pos, ref, alt, symbol, cq

def open_mutations(path, indels=False):
    ''' load mutations from a file.
    
    Args:
        path: path to table of de novo mutations.
        indels: whether to include indels or not.
    
    Returns:
        list of de novo mutations observed in current gene. Each entry is a
        dict with 'chrom', 'position', 'ref', 'alt', and 'consequence' keys.
    '''
    
    genes = {}
    with open(path) as handle:
        indices = parse_header(handle.readline())
        for line in handle:
            chrom, pos, ref, alt, symbol, cq = parse(line, indices)
            
            # ignore noncoding
            if cq not in MISSENSE_CQ | LOF_CQ:
                continue
            
            if not indels and (len(ref) > 1 or len(alt) > 1):
                continue
            
            if symbol not in genes:
                genes[symbol] = []
            
            genes[symbol].append({'chrom': chrom, 'pos': pos, 'ref': ref,
                'alt': alt, 'consequence': cq})
    
    return genes
