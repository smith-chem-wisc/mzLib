"""Reference values for TargetDecoyPep, from an independent implementation (scikit-learn's IsotonicRegression).

Writes, beside this script:
  target_decoy_pep_isotonic.tsv   value, weight, fitted          (the weighted isotonic regression alone)
  target_decoy_pep_search.tsv     score, decoy_fraction, pep     (a simulated search, end to end)

The simulated search has correct targets ~ N(3, 1), incorrect targets and decoys ~ N(0, 1), scores rounded to two
decimals so that ties occur, and some candidates half target, half decoy (decoy_fraction 0.5).

    python -I make_target_decoy_pep_fixtures.py
"""
import os

import numpy as np
from sklearn.isotonic import IsotonicRegression

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261103)


def isotonic(values, weights):
    return IsotonicRegression(increasing=True).fit_transform(np.arange(len(values)), values, sample_weight=weights)


# 1. The regression alone: values with violations, unequal positive weights.
values = np.round(rng.uniform(0, 1, 200), 6)
weights = np.round(rng.uniform(0.5, 5, 200), 6)
fitted = isotonic(values, weights)
with open(os.path.join(HERE, 'target_decoy_pep_isotonic.tsv'), 'w', newline='\n') as f:
    f.write('value\tweight\tfitted\n')
    for v, w, y in zip(values, weights, fitted):
        f.write(f'{float(v)!r}\t{float(w)!r}\t{float(y)!r}\n')

# 2. A simulated search.
correct = rng.normal(3, 1, 4000)
incorrect = rng.normal(0, 1, 1000)
decoys = rng.normal(0, 1, 1000)
half = rng.normal(0, 1, 50)
scores = np.round(np.concatenate([correct, incorrect, decoys, half]), 2)
fraction = np.concatenate([np.zeros(5000), np.ones(1000), np.full(50, 0.5)])
perm = rng.permutation(len(scores))
scores, fraction = scores[perm], fraction[perm]

distinct = np.unique(scores)[::-1]                      # best first
rate = np.array([fraction[scores == s].mean() for s in distinct])
size = np.array([np.sum(scores == s) for s in distinct], dtype=float)
fit = isotonic(rate, size)
block_pep = np.where(fit >= 0.5, 1.0, fit / (1.0 - fit))
pep_of = dict(zip(distinct, block_pep))
with open(os.path.join(HERE, 'target_decoy_pep_search.tsv'), 'w', newline='\n') as f:
    f.write('score\tdecoy_fraction\tpep\n')
    for s, d in zip(scores, fraction):
        f.write(f'{float(s)!r}\t{float(d)!r}\t{float(pep_of[s])!r}\n')
