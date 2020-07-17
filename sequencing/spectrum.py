from reference import *

def ideal_spectrum_amino(peptide: str, cyclic: bool) -> list:
    """Finds the ideal spectrum of a peptide

    :param peptide: the peptide (with 1-char amino acid codes in order)
    :type peptide: str
    :param cyclic: whether to include cyclic sub-peptides
    :type cyclic: bool
    :returns: all observed weights in the ideal spectrum, least->greatest
    :rtype: list
    """

    masses = [amino_to_weight[amino] for amino in peptide]
    return ideal_spectrum(masses, cyclic)

def ideal_spectrum(peptide: list, cyclic: bool) -> list:
    """Finds the ideal spectrum of a peptide

    :param peptide: the peptide (with amino masses in order)
    :type peptide: list
    :param cyclic: whether to include cyclic sub-peptides
    :type cyclic: bool
    :returns: all observed weights in the ideal spectrum, least->greatest
    :rtype: list
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
    :type pre_computed: dict
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

def sequence(spect: list) -> list:
    """Determines possible amino acid compositions given an ideal spectrum

    :param spect: the idea cyclic spectrum of the peptide
    :type spect: list
    :returns: all possible peptides as a lists of masses
    :rtype: list
    """
    

if __name__ == '__main__':
    print(count_peptides_with_mass(1347))
