from pprint import pprint
import msea
from msea import SetLibrary

## Download gene set databse
gmt_filepath = \
    'https://bitbucket.org/wangz10/msea/raw/aee6dd184e9bde152b4d7c2f3c7245efc1b80d23/msea/data/human_genes_associated_microbes/set_library.gmt'

d_gmt = msea.read_gmt(gmt_filepath)
print('Number of microbe-sets:', len(d_gmt))
# Look at a couple of reference sets in the library
#pprint(list(d_gmt.items())[:3])

microbe_set_input = set(['Colwellia',
                        'Deinococcus',
                        'Idiomarina',
                        'Neisseria',
                        'Pseudidiomarina',
                        'Pseudoalteromonas'])
# this can be done using the `msea.enrich` function
msea_result = msea.enrich(microbe_set_input, d_gmt=d_gmt, universe=1000)
print(' check the top enriched reference microbe-sets')
#print(msea_result.head())

# Step 2: perform MSEA with adjustment of expected ranks for reference sets. Sometimes certain reference microbe-sets in a library are more likely to be enriched by chance. We can adjust this by empirically estimating the null distributions of the ranks of the reference sets, then uses z-score to quantify if observed ranks are significantly different from the expected ranks.
set_library = SetLibrary.load(gmt_filepath)

# The SetLibrary.get_empirical_ranks() method helps compute the expected ranks and store the means and standard deviations of the ranks from the null distribution:
set_library.get_empirical_ranks(n=100)
print('Randomly sampled, without replacement, a universe of microbes under consideration. The random microbe-sets were then used for computing enrichment using Fisher’s exact test to estimate the expected ranks for the annotated microbe-sets. This procedure was repeated 10,000 times to compute the averages and standard deviations of the ranks for each annotated microbe-sets to compute a z-score for any future observed ranks from real microbe-set inputs')
print(set_library.rank_means.shape, set_library.rank_stds.shape)

# perform MSEA with this adjustment:
msea_result_adj = set_library.enrich(microbe_set_input, adjust=True, universe=1000)
print('perform MSEA with adjustment')
print(msea_result_adj.head())