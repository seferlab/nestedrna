import numpy as np
from scipy.stats import hypergeom


# hypergeometric test:
# N: The number of items in the population.
# k: The number of items in the population that are classified as successes.
# n: The number of items in the sample.
# x: The number of items in the sample that are classified as successes.


def hypergeometric_test(N, k, n, x):
    alpha = 0.95
    Prb = hypergeom.cdf(x, N, k, n)  # cumulative prob to draw this number
    rv = hypergeom(N, k, n)
    xx = x
    pmf_x = rv.pmf(xx)  # prob to draw this number
    interval_hg = hypergeom.interval(alpha, N, k, n, loc=0)  # confidence interval
    hypergeom_mean = hypergeom.mean(N, k, n, loc=0)

    print(f'\nWe drew # {int(x)} peaks of a certain HMM state,\n'
          f'cumulative probability is {np.round(Prb, 4)}')

    print(f'\nConf.interval low={int(interval_hg[0])},\n'
          f'Hypergeometric mean is {np.round(hypergeom_mean, 3)},\n'
          f'Conf.interval high={int(interval_hg[1])}')

    print(f'\nEnrichment folds for the given case is: {np.round(x/hypergeom_mean, 3)}')


N, k, n, x = 27263, 730, 766, 11
print(f'''Hypergeometric test for a given case: 
          N: {N} The number of items in the population
          k: {k} The number of items in the population that are classified as successes.
          n: {n} The number of items in the sample.
          x: {x} The number of items in the sample that are classified as successes.''')
hypergeometric_test(N, k, n, x)
