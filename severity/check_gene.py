
from denovonear.load_gene import load_gene
from denovonear.site_specific_rates import SiteRates
from denovonear.weights import WeightedChoice

from severity.open_severity import get_severity
from severity.simulation import analyse
from severity.regional_constraint import get_constrained_positions

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
    
    consequences = ['nonsense', 'missense', 'splice_lof']
    all_rates = {}
    for cq in consequences:
        all_rates[cq] = WeightedChoice()
    
    combined_tx = None
    for tx in transcripts:
        
        rates = SiteRates(tx, mut_dict, masked_sites=combined_tx, cds_coords=False)
        if combined_tx is None:
            combined_tx = tx
        else:
            combined_tx += tx
        
        for cq in consequences:
            all_rates[cq].append(rates[cq])
    
    return all_rates

def analyse_gene(ensembl, mut_dict, cadd, symbol, de_novos, constraint, weights):
    ''' analyse the severity of de novos found in a gene
    
    Args:
        ensembl: EnsemblRequest object, for transcript coordinates and sequence
        mut_dict: list of sequence-context mutation probabilities.
        cadd: pysam.TabixFile object for CADD scores (SNVs only)
        symbol: HGNC symbol for current gene
        de_novos: list of de novo mutations observed in current gene. Each entry
            is a dict with 'position', 'ref', 'alt', and 'consequence' keys.
        weights: dictionary of objects to weight CADD severity scores. We have
            different weights for protein-truncating and protein-altering
            variants, and within the protein-altering variants, different
            weights for variants in constrained and unconstrained regions.
    
    Returns:
        p-value for the observed total severity with respect to a null
        distribution of severities for the gene.
    '''
    
    sites = [ x['pos'] for x in de_novos ]
    try:
        # create gene/transcript for de novo mutations
        transcripts = load_gene(ensembl, symbol, sites)
    except IndexError:
        return 'NA'
    
    # get per site/allele mutation rates
    rates_by_cq = get_site_sampler(transcripts, mut_dict)
    
    chrom = transcripts[0].get_chrom()
    
    # get per site/allele severity scores, weighted by enrichment of missense
    # in known dominant at different severity thresholds
    constrained = get_constrained_positions(ensembl, constraint, symbol)
    severity = get_severity(cadd, chrom, rates_by_cq, weights, constrained)
    
    # convert the rates per site per consequence to rates per site
    rates = WeightedChoice()
    for cq in sorted(rates_by_cq):
        rates.append(rates_by_cq[cq])
    
    # get summed score for observed de novos
    observed = sum(( get_severity(cadd, chrom, de_novos, weights, constrained) ))
    
    # simulate distribution of summed scores within transcript
    return analyse(rates, severity, observed, len(de_novos), 1000000)
