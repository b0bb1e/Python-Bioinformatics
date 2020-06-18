from random import choice

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
    :rtype: list (of strs)
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

def ensure_validity(DNAs: list, pat_len: int):
    """For all motif finder methods, checks params for validity

    Raises appropriate errors if params are invalid
    
    :param DNAs: DNA strings with a shared motif
    :type DNAs: list (of strs)
    :param pat_len: the length of motifs to search for
    :type pat_len: int
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

def brute_finder(DNAs: list, pat_len: int, dist: int) -> set:
    """Brute-forces all possible motifs

    Trys every neighbor of every substring of length pat_len in the
    first DNA string
    
    :param DNAs: DNA strings with a shared motif
    :type DNAs: list (of strs)
    :param pat_len: the length of motifs to search for
    :type pat_len: int
    :param dist: the maximum substitutions for each instance of the motif
    :type dist: int
    :returns: all motifs conserved between all DNA strings
    :rtype: set (of strs)
    """

    ensure_validity(DNAs, pat_len)
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
    :rtype: list (of strs)
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
    :type DNAs: list (of strs)
    :param pat_len: the length of median strings to search for
    :type pat_len: int
    :returns: the best median string
    :rtype: list
    """
    
    ensure_validity(DNAs, pat_len)
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
            best_str, min_dist = cur_str, cur_dist
    return best_str

def calc_prob(pat: str, profile: list) -> float:
    """Calculates probability of a string given a profile

    :param pat: the string to find probability of
    :type pat: str
    :param profile: a probability profile (4 rows)
    :type profile: list (of lists (of floats))
    :returns: the probability of pat given profile
    :rtype: float
    """
    
    if not pat:
        raise ValueError('Cannot calculate probability of empty string')
    if not profile:
        raise ValueError('Cannot calculate probability without a profile')
    if len(profile) != 4:
        raise ValueError('Profiles must have 4 rows, one for each base')
    # profiles must be exactly as long (column wise) as the patterns
    pat_len = len(profile[0])
    if (pat_len != len(profile[1]) or pat_len != len(profile[2])
        or pat_len != len(profile[3])):
        raise ValueError('All rows of profile must be same length')
    if pat_len != len(pat):
        raise ValueError('Profile is the wrong length for this pattern')
    prob = 1
    for i in range(pat_len):
        try:
            prob *= profile[BASES.index(pat[i])][i]
        except:
            raise ValueError('Non-DNA base "' + pat[i] + '" found')
    return prob

def best_by_profile(DNA: str, profile: list) -> str:
    """Finds most probable substring given a profile

    :param DNA: the string to search in
    :type DNA: str
    :param profile: a probability profile, 4 rows & pat_len columns
    :type profile: list (of lists (of floats))
    :returns: the most probable substring of DNA
    :rtype: str
    """
    
    if not DNA:
        raise ValueError('Cannot search in empty string')
    if not profile:
        raise ValueError('Cannot search without a profile')
    if len(profile) != 4:
        raise ValueError('Profiles must have 4 rows, one for each base')
    pat_len = len(profile[0])
    if (pat_len != len(profile[1]) or pat_len != len(profile[2])
        or pat_len != len(profile[3])):
        raise ValueError('All rows of profile must be same length')
    if pat_len > len(DNA):
        raise ValueError('DNA string is not long enough for this profile')
    best_prob = -1.
    best_str = ''
    # try all possible subbstrings
    for i in range(len(DNA) - pat_len + 1):
        cur_str = DNA[i:i + pat_len]
        cur_prob = calc_prob(cur_str, profile)
        if cur_prob > best_prob:
            best_prob, best_str = cur_prob, cur_str
    return best_str

def get_profile(motifs: list) -> list:
    """Calculate a profile (with pseudocounts!) for some motifs

    :param motifs: DNA strings to build a profile off of
    :type motfis: list (of strs)
    :returns: a completed probability profile
    :rtype: list (of lists (of floats))
    """
    
    motif_len = len(motifs[0])
    profile = [[1 for _ in range(motif_len)] for __ in range(4)]
    for motif in motifs:
        for i in range(motif_len):
            try:
                profile[BASES.index(motif[i])][i] += 1
            except ValueError:
                raise ValueError('Non-DNA base "' + pat[i] + '" found')
    num_motifs = len(motifs) + 4
    # normalize column probabilities
    return [[(val / num_motifs) for val in row] for row in profile]

def consensus_string(profile: list) -> str:
    """Determines a consensus string from a profile

    :param profile: a probability profile, 4 rows
    :type profile: list (of lists (of floats))
    :returns: the most-probable string for this profile
    :rtype: str
    """
    
    con = ''
    for col in range(len(profile[0])):
        # assume A, try all others to prove better
        best_base = 'A'
        best_prob = profile[0][col]
        for row in range(1, 4):
            cur_prob = profile[row][col]
            if cur_prob > best_prob:
                best_base, best_prob = BASES[row], cur_prob
        con += best_base
    return con
        
def score_motifs(motifs: list, profile=None) -> int:
    """Scores motifs by their similarities to each other

    Lower scores = more similar

    :param motifs: the DNA motifs to score
    :type motifs: list (of strs)
    :param profile: a probability profile, 4 rows (default None)
    :type profile: list (of lists (of floats))
    :returns: a similarity score
    :rtype: int
    """
    
    score = 0
    if not profile:
        profile = get_profile(motifs)
    consensus = consensus_string(profile)
    for motif in motifs:
        # score is sum of differences from consensus
        score += ham_dist(motif, consensus)
    return score

def greedy_finder(DNAs: list, pat_len: int) -> list:
    """Use a greedy algorithm to find good motifs

    Moves from motifs -> median -> motifs, saving if better than last

    :param DNAs: DNA strings to search for shared motifs
    :type DNAs: list (of strs)
    :param pat_len: the length of motifs to search for
    :type pat_len: int
    :returns: each string's version of a motif
    :rtype: list (of strs)
    """

    ensure_validity(DNAs, pat_len)
    num_DNAs = len(DNAs)
    # assume first substrings of each string is best
    best_motifs = [DNA[0:pat_len] for DNA in DNAs]
    best_score = score_motifs(best_motifs)
    # try starting with each substring of the first DNA
    for i in range(len(DNAs[0]) - pat_len + 1):
        cur_motifs = [DNAs[0][i:i + pat_len]]
        # add on new motifs 1 at a time, by most prob
        for j in range(1, num_DNAs):
            cur_motifs.append(best_by_profile(DNAs[j], get_profile(cur_motifs)))
        cur_score = score_motifs(cur_motifs)
        # update if necessary
        if cur_score < best_score:
            best_motifs, best_score = cur_motifs, cur_score
    return best_motifs

def one_random_finder(DNAs: list, pat_len: int) -> list:
    """Run a randomized algorithm once to find decent motifs

    Moves from motifs -> median -> motifs, saving if better than last
    and leaving if not

    :param DNAs: DNA strings to search for shared motifs
    :type DNAs: list (of strs)
    :param pat_len: the length of motifs to search for
    :type pat_len: int
    :returns: each string's version of a motif
    :rtype: list (of strs)
    """
    
    num_DNAs = len(DNAs)
    best_motifs = []
    for i in range(num_DNAs):
        start = choice(range(len(DNAs[i]) - pat_len + 1))
        best_motifs.append(DNAs[i][start:start + pat_len])
    best_score = score_motifs(best_motifs)
    while True:
        profile = get_profile(best_motifs)
        cur_motifs = [best_by_profile(DNA, profile) for DNA in DNAs]
        cur_score = score_motifs(cur_motifs, profile)
        if cur_score < best_score:
            best_motifs, best_score = cur_motifs, cur_score
        else:
            return best_motifs, best_score

def random_finder(DNAs: list, pat_len: int) -> list:
    """Use a randomized algorithm to find good motifs

    Runs one_random_finder 10000 times, returning best result

    :param DNAs: DNA strings to search for shared motifs
    :type DNAs: list (of strs)
    :param pat_len: the length of motifs to search for
    :type pat_len: int
    :returns: each string's version of a motif
    :rtype: list (of strs)
    """

    ensure_validity(DNAs, pat_len)
    best_motifs, best_score = one_random_finder(DNAs, pat_len)
    for i in range(999):
        cur_motifs, cur_score = one_random_finder(DNAs, pat_len)
        if cur_score < best_score:
            best_motifs, best_score = cur_motifs, cur_score
    return best_motifs

if __name__ == '__main__':
    with open('data.txt') as data:
        pat_len, num_DNAs = [int(x) for x in data.readline().split()]
        DNAs = []
        for line in data:
            DNAs.append(line.rstrip())
    for motif in random_finder(DNAs, pat_len):
        print(motif)
