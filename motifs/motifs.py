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
            BASES.index(correct)
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
        for base in DNA:
            try:
                BASES.index(base)
            except:
                raise ValueError('Non-DNA base "' + base + '" in given string')
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

def all_DNA_strings(pat_len: int) -> list:
    """Produce all DNA strings of a certain length

    :param pat_len: the length of all strings to produce
    :type pat_len: int
    :returns: the DNA strings of length pat_len
    :rtype: list
    """
    
    if pat_len < 0:
        raise ValueError('Cannot find negative-length strings')
    if pat_len == 0:
        return ['']
    all_DNA = []
    # recursion; adding on all bases to each string of length less than
    for short in all_DNA_strings(pat_len - 1):
        for base in BASES:
            all_DNA.append(short + base)
    return all_DNA

def median_string(DNAs: list, pat_len: int) -> str:
    """Finds an optimal median string between all DNA strings

    A 'median string' is of length pat_len and appears in all DNA strings
    passed in with a minimal number of substitutions

    :param DNAs: DNA strings
    :type DNAs: list
    :param pat_len: the length of median strings to search for
    :type pat_len: int
    :returns: the best median string
    :rtype: list
    """
    
    if not DNAs:
        raise ValueError('Cannot find medians between non-existant strings')
    if len(DNAs) < 2:
        raise ValueError('Must compare at lest 2 strings to find medians')
    for DNA in DNAs:
        if not DNA:
            raise ValueError('Cannot use empty string as a median location')
    if pat_len < 1:
        raise ValueError('Medians must be at least 1 base long')
    min_dist = len(DNAs) * pat_len + 1
    best_str = ''
    for cur_str in all_DNA_strings(pat_len):
        # distance for cur_str
        cur_dist = 0
        for DNA in DNAs:
            # distance for this DNA string
            this_dist = pat_len + 1
            for i in range(len(DNA) - pat_len + 1):
                # distance for this substring
                dist = ham_dist(cur_str, DNA[i:i + pat_len])
                if dist < this_dist:
                    this_dist = dist
            cur_dist += this_dist
        if cur_dist < min_dist:
            best_str = cur_str
            min_dist = cur_dist
    return best_str
        
if __name__ == '__main__':
    with open('data.txt') as data:
        pat_len = int(data.readline().rstrip())
        DNAs = []
        for line in data:
            DNAs.append(line.rstrip())
    print(median_string(DNAs, pat_len))
