# Copyright IonQ Inc., all rights reserved

from itertools import combinations
from collections import defaultdict
import numpy as np


def collapse_histograms(hists):
    """
    Count occurences of each bitstring across all variants
    and their frequencies across all measurements

    hists - list of variant histograms
    occurences[k] - number of variants that measured bitstring k
    frequencies[k] - frequency of measuring k across all variants
    """
    nn = len(hists)
    occurences = defaultdict(lambda: 0)
    frequencies = defaultdict(lambda: 0)
    for h in hists:
        for k, v in h.items():
            occurences[k] += 1
            frequencies[k] += v / nn

    return dict(occurences), dict(frequencies)


def make_table(hists, bitstrings):
    """
    For given histograms and selected bitstrings,
    build a table of frequencies skipping the variants
    that do not contain any of the selected bitstrings

    hists - list of variant histograms
    bitstrings - selected bitrstrings
    table - frequencies of the selected bitstrings by variant
    """
    table = []
    for h in hists:
        variant = [h.get(k, 0) for k in bitstrings]
        if sum(variant) > 0:
            table.append(variant)
    return np.array(table)


def approximate_plurality_vote(hists, t=7):
    """
    Calculate relative probabilities to find each bitstring in at least t variants
    which can be found for each bitstring by summing the products of
    the probabilities to find the bitstring in the sampled variants (ll) times
    the probabiltiies not to find the bitstring in the remaining variants (rr)
    over (a) all samplings below the threshold and subtracting that from one
    which is equivalent to (b) summing over all samplings above the threshold
    which is faster if the threshold is smaller than half of the number of variants

    a) prob_k = 1 - sum_{p in P(nn,m < t)} prod_{i in p} h_ik * prod_{j not in p} (1-h_jk)
    b) prob_k = sum_{p in P(nn,t <= m)} prod_{i in p} h_ik * prod_{j not in p} (1-h_jk)
    where nn is the number of variants and P(nn,m) are all samplings of m out of nn

    hists - list of variant histograms
    t - bitstring filtering threshold
    """
    occurrences, avg_frequencies = collapse_histograms(hists)

    # adjust the threshold
    max_occurrence = max(occurrences.values())
    t_max = min(t, max_occurrence)

    # return a bitstring-wise average if the threshold <= 2
    if t_max <= 2:
        return avg_frequencies

    # filter the bitstrings and return the result if there is only one
    bitstrings = [k for k, v in occurrences.items() if v >= t_max]
    if len(bitstrings) == 1:
        return {bitstrings[0]: 1}

    # continue with non-trivial cases and build a table of frequencies
    table = make_table(hists, bitstrings)
    nn, nk = table.shape

    # pick the best summation method
    if t_max < nn / 2:
        # use formula (a)
        mm = range(t_max)
        pp = np.ones(nk)
        s = -1
    else:
        # use formula (b)
        mm = range(t_max, nn + 1)
        pp = np.zeros(nk)
        s = 1

    # find the sum of all contributions for each bitstring
    for m in mm:
        ll = np.asarray(list(combinations(range(nn), m)), dtype=int)
        rr = np.asarray(list(combinations(range(nn), nn - m))[::-1], dtype=int)
        pp += s * np.sum(
            np.prod(table[ll, :], axis=1) * np.prod(1 - table[rr, :], axis=1), axis=0
        )

    # return normilized distribution
    norm = sum(pp)
    return {k: pp[i] / norm for i, k in enumerate(bitstrings) if pp[i] > 0}


def test_uni_vrs(st=4, du=0.02, nv=25):
    """
    Generate identical variant histograms with
    a skewed uniform distribution

    st - number of states
    du - skew per state
    nv - number of variants
    """
    u = 1 / st
    dx = (st - 1) / 2
    vv = [{str(x): u + (x - dx) * du for x in range(st)}] * nv
    return vv


def main():
    vv = test_uni_vrs()
    for t in range(2, 8):
        res = approximate_plurality_vote(vv, t=t)
        print({k: round(res[k], 5) for k in res})


if __name__ == "__main__":
    main()
