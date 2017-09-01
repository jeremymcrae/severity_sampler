

import pandas
pandas.set_option('display.width', 170)
pandas.set_option('display.max_rows', 200)

import matplotlib
matplotlib.use('agg')

from matplotlib import pyplot
from numpy import log10

from mupit.combine_analyses import fishersMethod

matplotlib.rcParams.update({'font.size': 16})

enrichment_path = '/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.all.2015-10-12.txt'

enrichment = pandas.read_table(enrichment_path)
enrichment = enrichment[['hgnc', 'chrom', 'ddd.p_lof', 'ddd.p_func',
    'ddd.p_missense_clust', 'ddd.p_combined', 'ddd.p_min']]

severity_path = '/nfs/users/nfs_j/jm33/severity_sampler/results.cadd_scores.txt'
severity = pandas.read_table(severity_path)
recode = dict(zip(severity['symbol'], severity['severity_p_value']))

enrichment['p_severity'] = enrichment['hgnc'].map(recode)

# combine the functional enrichment, severity and clustering
p_values = enrichment[["ddd.p_func", "p_severity", 'ddd.p_missense_clust']]
enrichment["p_func_sev_clust"] = p_values.apply(fishersMethod, axis=1)

# check our top hitter
# NOTE: ARID1B has a much lower lower p-value now for two reasons:
# NOTE: 1) the severity p-values are currently capped at 1e-6. Additional
# NOTE:    simulations may get a more accurate p-value estimate. Aalthough, if
# NOTE:    the p-value is below 1e-6, it's likely to be genome-wide significant
# NOTE:    anyway, so there might not be much point in going deeper.
# NOTE: 2) The severity p-value is for SNVs only, so the functional impact
# NOTE:    appears to be smaller for genes with large numbers of frameshift
# NOTE:    mutations.
enrichment[enrichment['hgnc'] == 'ARID1B']

# check some anomalous variants. Most genes missing from previous significance
# are because they have numerous frameshift_variants, which are skipped in the
# severity simulations, reducing the apparent severity.
enrichment[enrichment['hgnc'] == 'PPM1D']
enrichment[enrichment['hgnc'] == 'AUTS2']

enrichment.sort_values('p_func_sev_clust').head(200)

thresh = 1e-6
combined = set(enrichment['hgnc'][enrichment['p_func_sev_clust'] < thresh])
func_clust = set(enrichment['hgnc'][enrichment['ddd.p_combined'] < thresh])
standard = set(enrichment['hgnc'][(enrichment['ddd.p_min'] < thresh)])

standard - combined
combined - standard

fig = pyplot.figure(figsize=(6, 6))
ax = fig.gca()

e = ax.plot(-log10(enrichment['ddd.p_min']), -log10(enrichment['p_func_sev_clust']),
    # color='blue',
    linestyle='None', marker='.')

e = ax.set_xlabel('-log10(P min) in DDD 4K')
e = ax.set_ylabel('-log10(P combined) in DDD 4K')

high = max(ax.get_xlim()[1], ax.get_ylim()[1])
low = min(ax.get_xlim()[0], ax.get_ylim()[0])

e = ax.plot([low, high], [low, high], color='grey', linestyle='--')

e = ax.spines['top'].set_visible(False)
e = ax.spines['right'].set_visible(False)

e = ax.spines['left'].set_smart_bounds(True)
e = ax.spines['bottom'].set_smart_bounds(True)

fig.savefig('comparison.pdf', format='pdf', transparent=True, bbox_inches='tight',
    pad_inches=0, frameon=False)
