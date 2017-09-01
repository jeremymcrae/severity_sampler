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

import argparse

import pysam

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates

from severity.open_mutations import open_mutations
from severity.regional_constraint import load_regional_constraint
from severity.check_gene import analyse_gene
from severity.weights import weights as WEIGHTS

def get_options():
    parser = argparse.ArgumentParser('')
    parser.add_argument('--de-novos',
        default='/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt',
        help='path to table of de novo mutations. Table must contain columns '
            'named chrom, pos, ref, alt, symbol, and consequence.'),
    parser.add_argument('--cadd',
        default='/lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz',
        help='Path to tabix-indexed CADD scores for all SNVs.')
    parser.add_argument('--constraint',
        help='Path to table of regional constraint.')
    parser.add_argument('--cache', default='cache',
        help='Path to cache transcript coordinates and sequence from Ensembl.')
    parser.add_argument('--genome-build', default='grch37',
        help='Genome build for coordinates from Ensembl.')
    
    parser.add_argument('-o', '--output', default='results.txt',
        help='Path to write output results to.')
    
    return parser.parse_args()

def main():
    args = get_options()
    
    ensembl = EnsemblRequest(args.cache, args.genome_build)
    cadd = pysam.TabixFile(args.cadd)
    
    constraint = load_regional_constraint(args.constraint)
    
    # open de novo mutations
    all_de_novos = open_mutations(args.de_novos)
    
    mut_dict = load_mutation_rates()
    
    output = open(args.output, 'w')
    output.write('symbol\tseverity_p_value\n')
    for symbol in all_de_novos:
        if symbol in ['', '.']:
            continue
        
        print(symbol)
        de_novos = all_de_novos[symbol]
        p_value = analyse_gene(ensembl, mut_dict, cadd, symbol, de_novos, constraint, WEIGHTS)
        line = '{}\t{}\n'.format(symbol, p_value)
        output.write(line)

if __name__ == '__main__':
    main()
