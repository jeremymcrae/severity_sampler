
import pysam

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import load_gene
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.weights import WeightedChoice

from severity.open_mutations import open_mutations
from severity.open_severity import get_severity
from severity.simulate import analyse

def get_site_sampler(transcripts, mut_dict):
    '''
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
    
    return all_rates, combined_tx

def main():
    ensembl = EnsemblRequest(cache_folder='cache', genome_build='grch37')
    
    cadd_path = '/lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz'
    cadd = pysam.TabixFile(cadd_path)
    
    # open de novo mutations
    path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2016-04-19.txt'
    all_de_novos = open_mutations(path)

    mut_dict = load_mutation_rates()

    for symbol in all_de_novos:
        
        if symbol in ['', '.']:
            continue
        
        print(symbol)
        de_novos = all_de_novos[symbol]
        sites = [ x['position'] for x in de_novos ]
        try:
            # create gene/transcript for de novo mutations
            transcripts = load_gene(ensembl, symbol, sites)
        except IndexError:
            continue
        
        # get per site/allele mutation rates
        rates, tx = get_site_sampler(transcripts, mut_dict)
        
        chrom = transcripts[0].get_chrom()
        
        # get per site/allele severity scores
        severity = [ get_severity(cadd, chrom, x['pos'], x['alt']) for x in rates ]
        
        # TODO: weight severity scores
        
        # get summed score for observed de novos
        observed = sum(( get_severity(cadd, chrom, x['position'], x['alt']) for x in de_novos ))
        
        # simulate distribution of summed scores within transcript
        analyse(deref(rates.thisptr), severity, observed, len(de_novos), 100000000)

if __name__ == '__main__':
    main()
