BASES = ('A', 'C', 'G', 'T')

def ham_dist(one: str, two: str) -> int:
    """Calculates HAMming DISTance between two strings

    'Hamming distance' is the number of substituted chars

    :param one: a string
    :type one: str
    :param two: a string to compare to one
    :type two: str
    :returns: the substitution distance between one and two
    :rtype: int
    """
    
    if not one or not two:
        raise ValueError('Cannot find distance between empty strings')
    if len(one) != len(two):
        raise ValueError('Strings must be of the same length')
    dist = 0
    for i in range(len(one)):
        if one[i] != two[i]:
            dist += 1
    return dist

def get_neighbors(pat: str, dist: int) -> list:
    """Finds all neighbors of a DNA string

    'Neighbors' are strings of the same length with at most dist
    base substitutions. This includes the original string

    :param pat: the DNA string to find neighbors of
    :type pat: str
    :param dist: the maximum number of changes to allow
    :type dist: int
    :returns: all neighbors of pat
    :rtype: list
    """
    
    if not pat:
        raise ValueError('Cannot find neighbors of empty string')
    if dist < 0:
        raise ValueError('Cannot have a negative distance')
    # dist of 0 means no changes
    if dist == 0:
        return [pat]
    if len(pat) == 1:
        return list(BASES)
    # lists of strings currently 0, 1, 2, .... d differences from pat
    neighbors = [[] for _ in range(dist + 1)]
    # build off empty string
    neighbors[0] = ['']
    # add on new chars to each existing string one at a time
    for i in range(len(pat)):
        old = list(neighbors)
        neighbors = [[] for _ in range(dist + 1)]
        correct = pat[i]
        try:
            print(BASES.index(correct))
        except ValueError:
            raise ValueError('Non-DNA base "' + correct + '" in given string')
        # loop over all dist-lists
        for cur_dist in range(dist + 1):
            for cur_pat in old[cur_dist]:
                # using the correct char does not move back a list
                neighbors[cur_dist].append(cur_pat + correct)
                if cur_dist < dist:
                    for base in BASES:
                        if base != correct:
                            # using incorrect char does move back a list
                            neighbors[cur_dist + 1].append(cur_pat + base)
    # group together in one list
    all_neighbors = []
    for group in neighbors:
        for neighbor in group:
            all_neighbors.append(neighbor)
    return all_neighbors

def brute_finder(DNAs: list, pat_len: int, dist: int) -> set:
    """Brute-forces all possible motifs

    Trys every neighbor of every substring of length pat_len in the
    first DNA string
    
    :param DNAs: DNA strings with a shared motif
    :type DNAs: list
    :param pat_len: the length of motifs to search for
    :type pat_len: int
    :param dist: the maximum substitutions for each instance of the motif
    :type dist: int
    :returns: all motifs conserved between all DNA strings
    :rtype: set
    """

    if not DNAs:
        raise ValueError('Cannot find motifs between non-existant strings')
    if len(DNAs) < 2:
        raise ValueError('Must compare at lest 2 strings to find motifs')
    for DNA in DNAs:
        if not DNA:
            raise ValueError('Cannot use empty string as a motif location')
    if pat_len < 1:
        raise ValueError('Motifs must be at least 1 base long')
    motifs = set()
    end_indexes = [(len(DNA) - pat_len + 1) for DNA in DNAs]
    for i in range(end_indexes[0]):
        cur_pat = DNAs[0][i:i + pat_len]
        if not cur_pat in motifs:
            for neighbor in get_neighbors(cur_pat, dist):
                if not neighbor in motifs:
                    found_all = True
                    for j in range(1, len(DNAs)):
                        found = False
                        for k in range(end_indexes[j]):
                            compare_pat = DNAs[j][k:k + pat_len]
                            if ham_dist(neighbor, compare_pat) <= dist:
                                found = True
                                break
                        if not found:
                            found_all = False
                            break
                    if found_all:
                        motifs.add(neighbor)
    return motifs
