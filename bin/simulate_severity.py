
import time

import pysam

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import load_gene
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.weights import WeightedChoice

from severity.open_mutations import open_mutations
from severity.open_severity import get_severity
from severity.simulation import analyse

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
        continue
    
    # get per site/allele mutation rates
    rates = get_site_sampler(transcripts, mut_dict)
    
    chrom = transcripts[0].get_chrom()
    
    # get per site/allele severity scores
    severity = get_severity(cadd, chrom, rates)
    
    # TODO: weight severity scores
    
    # get summed score for observed de novos
    observed = sum(( get_severity(cadd, chrom, de_novos) ))
    
    start = time.time()
    # simulate distribution of summed scores within transcript
    p_value = analyse(rates, severity, observed, len(de_novos), 1000000)
    elapsed = time.time() - start
    
    return p_value, elapsed

def main():
    ensembl = EnsemblRequest(cache_folder='cache', genome_build='grch37')
    
    cadd_path = '/lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz'
    cadd = pysam.TabixFile(cadd_path)
    
    # open de novo mutations
    path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt'
    all_de_novos = open_mutations(path)
    
    mut_dict = load_mutation_rates()
    
    output = open('severity_results.txt', 'w')
    output.write('symbol\tp_value\tn_de_novos\tcpu_time\n')
    for symbol in all_de_novos:
        
        if symbol in ['', '.']:
            continue
        
        print(symbol)
        de_novos = all_de_novos[symbol]
        p_value, elapsed = analyse_gene(ensembl, mut_dict, cadd, symbol, de_novos)
        line = '{}\t{}\t{}\t{}\n'.format(symbol, p_value, len(de_novos), elapsed)
        output.write(line)

if __name__ == '__main__':
    main()
