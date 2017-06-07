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
from denovonear.load_gene import load_gene
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.weights import WeightedChoice

from severity.open_mutations import open_mutations
from severity.open_severity import get_severity
from severity.simulation import analyse

def get_options():
    parser = argparse.ArgumentParser('')
    parser.add_argument('--de-novos',
        default='/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt',
        help='path to table of de novo mutations. Table must contain columns ' \
            'named chrom, pos, ref, alt, symbol, and consequence.'),
    parser.add_argument('--cadd',
        default='/lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz',
        help='Path to tabix-indexed CADD scores for all SNVs.')
    parser.add_argument('--cache', default='cache',
        help='Path to cache transcript coordinates and sequence from Ensembl.')
    parser.add_argument('--genome-build', default='grch37',
        help='Genome build for coordinates from Ensembl.')
    
    parser.add_argument('-o', '--output', default='results.txt',
        help='Path to write output results to.')
    
    return parser.parse_args()

def get_site_sampler(transcripts, mut_dict):
    ''' get per position and alt allele mutation probability sampler.
    
    We need to be able to sample each site within a gene, where the probability
    of sampling a given site is equal to the sequence-context derived mutation
    probability. We use the denovonear.weights.WeightedChoice for this, which
    wraps around a cpp class for quick sampling. We use the SiteRates class to
    derive the per site/allele probabilities for different consequence
    categories. We combine the categories of interest intoa single object, so we
    can sample across the full transcript at once. This also allows for multiple
    transcripts for a single gene, by taking the union of transcripts.
    
    Args:
        transcripts: list of Transcript objects for a gene.
        mut_dict: list of sequence-context mutation probabilities.
    
    Returns:
        denovonear.WeightedChoice object, containing the mutation probabilities
        per position and alt allele.
    '''
    
    all_rates = WeightedChoice()
    combined_tx = None
    for tx in transcripts:
        
        rates = SiteRates(tx, mut_dict, masked_sites=combined_tx, cds_coords=False)
        if combined_tx is None:
            combined_tx = tx
        else:
            combined_tx += tx
        
        for cq in ['synonymous', 'nonsense', 'missense', 'splice_lof']:
            all_rates.append(rates[cq])
    
    return all_rates

def analyse_gene(ensembl, mut_dict, cadd, symbol, de_novos):
    ''' analyse the severity of de novos found in a gene
    
    Args:
        ensembl: EnsemblRequest object, for transcript coordinates and sequence
        mut_dict: list of sequence-context mutation probabilities.
        cadd: pysam.TabixFile object for CADD scores (SNVs only)
        symbol: HGNC symbol for current gene
        de_novos: list of de novo mutations observed in current gene. Each entry
            is a dict with 'position', 'ref', 'alt', and 'consequence' keys.
    
    Returns:
        p-value for the observed total severity with respect to a null
        distribution of severities for the gene.
    '''
    
    sites = [ x['position'] for x in de_novos ]
    try:
        # create gene/transcript for de novo mutations
        transcripts = load_gene(ensembl, symbol, sites)
    except IndexError:
        return 'NA'
    
    # get per site/allele mutation rates
    rates = get_site_sampler(transcripts, mut_dict)
    
    chrom = transcripts[0].get_chrom()
    
    # get per site/allele severity scores
    severity = get_severity(cadd, chrom, rates)
    
    # TODO: weight severity scores
    
    # get summed score for observed de novos
    observed = sum(( get_severity(cadd, chrom, de_novos) ))
    
    # simulate distribution of summed scores within transcript
    return analyse(rates, severity, observed, len(de_novos), 1000000)

def main():
    args = get_options()
    
    ensembl = EnsemblRequest(args.cache, args.genome_build)
    cadd = pysam.TabixFile(args.cadd)
    
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
        p_value = analyse_gene(ensembl, mut_dict, cadd, symbol, de_novos)
        line = '{}\t{}\n'.format(symbol, p_value)
        output.write(line)

if __name__ == '__main__':
    main()
