from reference import *
from collections import Counter

def ideal_spectrum_amino(peptide: str, cyclic: bool) -> list:
    """Finds the ideal spectrum of a peptide

    :param peptide: the peptide (with 1-char amino acid codes in order)
    :type peptide: str
    :param cyclic: whether to include cyclic sub-peptides
    :type cyclic: bool
    :returns: all observed weights in the ideal spectrum, least->greatest
    :rtype: list (of ints)
    """

    masses = [amino_to_weight[amino] for amino in peptide]
    return ideal_spectrum(masses, cyclic)

def ideal_spectrum(peptide: list, cyclic: bool) -> list:
    """Finds the ideal spectrum of a peptide

    :param peptide: the peptide (with amino masses in order)
    :type peptide: list (of ints)
    :param cyclic: whether to include cyclic sub-peptides
    :type cyclic: bool
    :returns: all observed weights in the ideal spectrum, least->greatest
    :rtype: list (of ints)
    """

    cum_mass = [0]
    mass_len = 1
    # fill cum_mass with the cumalitive mass of peptide prefixes
    for mass in peptide:
        cum_mass.append(cum_mass[mass_len - 1] + mass)
        mass_len += 1

    # the empty peptide is included in the spectrum
    spect = [0]
    # for every pair of cum_mass values, add difference to spect
    spect += [cum_mass[j] - cum_mass[i]
             for i in range(mass_len) for j in range(i + 1, mass_len)]

    if cyclic:
        # the final cumalative mass value is simply the total mass
        total_mass = cum_mass[mass_len - 1]
        # for every pair of cum_mass values that doesn't contain an end
        # add the cyclic peptide, or the total mass minus the above difference
        spect += [total_mass - (cum_mass[j] - cum_mass[i])
                  for i in range(1, mass_len)
                  for j in range(i + 1, mass_len - 1)]

    spect.sort()
    return spect

def count_peptides_with_mass(mass: int, pre_computed: dict={}) -> int:
    """Counts the number of peptides with a given mass

    :param mass: the mass of the peptides to count
    :type mass: int
    :param pre_computed: values already computed, with mass:return_val
                         defaults to empty
    :type pre_computed: dict (int: int)
    :returns: the number of peptides with this mass
    :rtype: int
    """

    if mass < 0:
        return 0
    elif mass == 0:
        return 1
    elif mass in pre_computed:
        return pre_computed[mass]
    else:
        count = 0
        for m in amino_masses:
            count += count_peptides_with_mass(mass - m, pre_computed)
        pre_computed[mass] = count
        return count

def expand(peptides: list, masses: list=amino_masses) -> list:
    """Expands each peptide in in a list by all possible next-masses

    :param peptides: the peptides to expand
    :type peptides: list (of lists (of ints))
    :param masses: the list of masses to expand by (default amino_masses)
    :type masses: list (of ints)
    :returns: the expanded list
    :rtype: list (of lists (of ints))
    """
    
    expanded = []
    for peptide in peptides:
        for mass in masses:
            expanded.append(peptide + [mass])
    return expanded

def contains(big: list, small: list) -> bool:
    """Determines if one list is a superset of another

    By using Counter and not converting to a set, this is sensative
    to the number of elements (if small has 2 'a', big must have at
    least 2 'a' and not just one)

    :param big: the list that may contain another
    :type big: list
    :param small: the list that may be contained
    :type small: list
    :returns: if big is a supersert of small
    :rtype: bool
    """
    
    return not Counter(small) - Counter(big)

def sequence(spect: list) -> list:
    """Determines possible amino acid compositions given an ideal spectrum

    :param spect: the ideal cyclic spectrum of the peptide
    :type spect: list (of ints)
    :returns: all possible peptides as a lists of masses
    :rtype: list (of lists (of ints))
    """

    all_peptides = [[]]
    good_peptides = []
    total_mass = spect[len(spect) - 1]
    while all_peptides:
        # branch step
        all_peptides = expand(all_peptides)
        # iterate over list backwards for easier removal
        for i in range(len(all_peptides) - 1, -1, -1):
            # if masses match, this is a possible end-answer
            if sum(all_peptides[i]) == total_mass:
                # only check cyclic spectrum if peptide is possible answer
                if ideal_spectrum(all_peptides[i], True) == spect:
                    good_peptides.append(all_peptides[i])
                # no matter what, stop considering this peptide
                del all_peptides[i]
            # remove intermediate peptides if not consistant with spect
            elif not contains(spect, ideal_spectrum(all_peptides[i], False)):
                del all_peptides[i]
    return good_peptides

def score_amino(peptide: list, spect: list, cyclic: bool):
    """Scores a cyclic peptide against a spectrum

    A higher score indicates more similarity due to more shared massses

    :param peptide: the peptide (with 1-char amino acid codes in order)
    :type peptide: str
    :param spect: a cyclic spectrum to score against
    :type spect: list (of ints)
    :param cyclic: whether to deal with peptides as cyclic or not
    :type cyclic: bool
    :returns: the number of matching masses between peptide & spect
    :rtype: int
    """

    masses = [amino_to_weight[amino] for amino in peptide]
    return score(masses, spect, cyclic)

def score(peptide: list, spect: list, cyclic: bool):
    """Scores a cyclic peptide against a spectrum

    A higher score indicates more similarity due to more shared massses

    :param peptide: the peptide to score as a list of masses
    :type peptide: list (of ints)
    :param spect: a cyclic spectrum to score against
    :type spect: list (of ints)
    :param cyclic: whether to deal with peptides as cyclic or not
    :type cyclic: bool
    :returns: the number of matching masses between peptide & spect
    :rtype: int
    """

    matches = Counter(ideal_spectrum(peptide, cyclic)) & Counter(spect)
    return sum(matches.values())

def trim(leaderboard: list, spect: list, keep: int) -> list:
    """Trims a leaderboard of peptides to top-keep-plus-ties

    :param leaderboard: the peptide leaderboard to trim, with peptides
                        as lists of masses
    :type leaderboard: list (of lists (of ints))
    :param keep: the minimum number of peptides to keep (if no ties)
    :param spect: a cyclic spectrum to score against
    :type spect: list (of ints)
    :type keep: int
    :returns: the trimmed leaderboard
    :rtype: list (of lists (of ints))
    """

    if len(leaderboard) <= keep:
        return leaderboard

    score_count = {}
    all_scores = []
    for peptide in leaderboard:
        cur_score = score(peptide, spect, False)
        if cur_score in score_count:
            score_count[cur_score].append(peptide)
        else:
            score_count[cur_score] = [peptide]
            all_scores.append(cur_score)

    all_scores.sort(reverse=True)
    trimmed = []
    for cur_score in all_scores:
        trimmed += score_count[cur_score]
        if len(trimmed) >= keep:
            break
    return trimmed

def top_diffs(spect: list, num_acids: int) -> list:
    spect.sort()
    spect_len = len(spect)
    diffs = [spect[i] - spect[j] for i in range(1, spect_len)
             for j in range(i - 1, -1, -1)]
    diff_count = Counter(diffs)
    if 0 in diff_count:
        del diff_count[0]
    tops = diff_count.most_common()
    acids = []
    last_count = 0
    for mass, count in tops:
        if len(acids) >= num_acids and count < last_count:
            break
        if 57 <= mass <= 200:
            acids.append(mass)
            last_count = count
    return acids
    
def leaderboard_sequence(spect: list, keep: int,
                         masses: list=amino_masses) -> list:
    """Determines the most-matching peptide given a spectrum

    :param spect: an observed, perhaps non-ideal, spectrum
    :type spect: list (of ints)
    :param keep: the number of peptides to keep after each trim
    :type keep: int
    :param masses: the list of masses to expand by (default amino_masses)
    :type masses: list (of ints)
    :returns: the most-matching peptide
    :rtype: list (of ints)
    """
    print(len(spect))
    leaderboard = [[]]
    best_peptide = []
    best_score = -1
    total_mass = spect[len(spect) - 1]
    while leaderboard:
        leaderboard = expand(leaderboard, masses)
        # iterate over list backwards for easier removal
        for i in range(len(leaderboard) - 1, -1, -1):
            cur_mass = sum(leaderboard[i])
            if cur_mass == total_mass:
                cur_score = score(leaderboard[i], spect, True)
                #cur_score -= score(leaderboard[i], spect, False)
                
                if cur_score >= best_score:
                    print(leaderboard[i], cur_score)
                    best_peptide = leaderboard[i]
                    best_score = cur_score
                del leaderboard[i]
            elif cur_mass > total_mass:
                del leaderboard[i]
        leaderboard = trim(leaderboard, spect, keep)
    return best_peptide

def convolution_sequence(spect: list, keep: int, num_acids: int) -> list:
    """Determines the most-matching peptide given a spectrum

    Amino acid masses are taken from a convolution

    :param spect: an observed, perhaps non-ideal, spectrum
    :type spect: list (of ints)
    :param keep: the number of peptides to keep after each trim
    :type keep: int
    :param num_acids: the number of acids to keep from the convolution
    :returns: the most-matching peptide
    :rtype: list (of ints)
    """

    return leaderboard_sequence(spect, keep, top_diffs(spect, num_acids))

if __name__ == '__main__':
    with open('data.txt') as data:
        num_acids = int(data.readline().rstrip())
        keep = int(data.readline().rstrip())
        spect = [int(x) for x in data.readline().rstrip().split()]
    print(*convolution_sequence(spect, keep, num_acids), sep='-')
