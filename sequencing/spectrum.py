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
        
    return sorted(spect)

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

def expand(peptides: list) -> list:
    """Expands each peptide in in a list by all possible next-masses

    :param peptides: the peptides to expand
    :type peptides: list (of lists (of ints))
    :returns: the expanded list
    :rtype: list (of lists (of ints))
    """
    
    expanded = []
    for peptide in peptides:
        for mass in amino_masses:
            expanded.append(peptide + [mass])
    return expanded

def sum_list(nums: list) -> int:
    """Sums a list of ints

    :param nums: the list to sum elements of
    :type nums: list (of ints)
    :returns: the sum of all of nums' elements
    :rtype: int
    """
    
    count = 0
    for num in nums:
        count += num
    return count

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

    :param spect: the idea cyclic spectrum of the peptide
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
            if sum_list(all_peptides[i]) == total_mass:
                # only check cyclic spectrum if peptide is possible answer
                if ideal_spectrum(all_peptides[i], True) == spect:
                    good_peptides.append(all_peptides[i])
                # no matter what, stop considering this peptide
                del all_peptides[i]
            # remove intermediate peptides if not consistant with spect
            elif not contains(spect, ideal_spectrum(all_peptides[i], False)):
                del all_peptides[i]
    return good_peptides

if __name__ == '__main__':
    with open('data.txt') as data:
        spect = [int(x) for x in data.readline().rstrip().split()]
    for peptide in sequence(spect):
        print(*peptide, sep='-', end=' ')
