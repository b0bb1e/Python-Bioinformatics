BASES = ('A', 'C', 'G', 'T')
# matches bases with their complements
COMP = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def DNA_to_num(DNA: str) -> int:
    """Converts a DNA string to a number

    For all strings of a certain length, lower numbers will correspond
    to an earlier point in the alphabet. Numbers are re-used in-between
    lengths.

    :param DNA: the DNA string of ACGT to convert
    :type DNA: str
    :returns: a number corresponding to the DNA string
    :rtype: int
    """

    if not DNA:
        raise ValueError('Cannot convert empty string to number')
    num = 0
    chars = len(DNA)
    for i in range(chars):
        try:
            # earlier chars are given more weight/value in the number
            num += BASES.index(DNA[i]) * (4 ** (chars - i - 1))
        except ValueError:
            raise ValueError('Non-DNA base "' + DNA[i] + '" in given string')
    return num

def num_to_DNA(num: int, chars: int) -> str:
    """Converts a number to a DNA string

    :param num: the number to convert
    :type num: int
    :param chars: the length of the eventual string
    :type chars: int
    :returns: the DNA string of length chars equivalent to num
    :rtype: str
    """

    if num < 0:
        raise ValueError('DNA numbers must be non-negative')
    if chars < 1:
        raise ValueError('DNA strings must be at least length 1')
    DNA = ''
    for i in range(chars):
        # the smaller numbers, pulled off first, go in the back
        DNA = BASES[int(num % 4)] + DNA
        num = num // 4
    if num != 0:
        raise ValueError('Number is too large for length given')
    return DNA

def find_starts(DNA: str, pat: str) -> list:
    """Find all start indexes of a substring

    :param DNA: the longer string to search in
    :type DNA: str
    :param pattern: the substring to search for
    :type pattern: str
    :returns: all indexes where pattern starts in DNA
    :rtype: list
    """
    
    if not DNA:
        raise ValueError('Cannot search in empty string')
    if not pat:
        raise ValueError('Cannot search for empty string')
    pat_len = len(pat)
    DNA_len = len(DNA)
    if pat_len > DNA_len:
        return []
    starts = []
    # only search possible starting positions
    for i in range(DNA_len - pat_len + 1):
        if DNA[i:i + pat_len] == pat:
            starts.append(i)
    return starts

def rev_comp(pat: str) -> str:
    """Finds the reverse complement of a DNA string

    :param pat: the DNA string to REVerse COMPlement
    :type pat: str
    :returns: the reverse complement of pat
    :rtype: str
    """

    if not pat:
        raise ValueError('Cannot reverse-complement empty string')
    rev = ''
    for base in pat:
        try:
            rev = COMP.get(base) + rev
        except TypeError:
            raise ValueError('Non-DNA base "' + base + '" in given string')
    return rev

def min_skew(DNA: str) -> list:
    """Finds all indexes of minimum skew in a DNA string

    'Minimum skew' is when #G - #C is at a minimum. Other bases are ignored

    :param DNA: the DNA string to calculate skew for
    :type DNA: str
    :returns: all indexes where #G - #C is minimized
    :rtype: list:
    """
    
    if not DNA:
        raise ValueError('Cannot search in empty string')
    mins = []
    min_skew = 0
    skew = 0
    for i in range(len(DNA)):
        if DNA[i] == 'C':
            skew -= 1
        elif DNA[i] == 'G':
            skew += 1
        elif DNA[i] != 'A' and DNA[i] != 'T':
            raise ValueError('Non-DNA base "' + DNA[i] + '" in given string')
        if skew < min_skew:
            min_skew = skew
            mins = [i]
        elif skew == min_skew:
            mins.append(i)
    return mins

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

def calc_freq(DNA: str, pat_len: int) -> list:
    """Calculates the number of times each string of a certain length appears

    :param DNA: the longer string to search in
    :type DNA: str
    :param pat_len: the length of substrings to count
    :type pat_len: int
    :returns: a list where [i] is the times the string num_to_DNA(i, pat_len)
              appears
    :rtype: list
    """
    
    if not DNA:
        raise ValueError('Cannot search in empty string')
    if pat_len < 1:
        raise ValueError('Patterns must have length at least 1')
    freq = [0 for _ in range(4 ** pat_len)]
    # search all starting positions
    for i in range(len(DNA) - pat_len + 1):
        freq[DNA_to_num(DNA[i:i + pat_len])] += 1
    return freq

def find_most_freq(DNA: str, pat_len: int, min_times: int = 2,
                   all_above: bool = False) -> list:
    """Find substrings which appear with enough frequency

    :param DNA: the longer string to search in
    :type DNA: str
    :param pat_len: the length of substrings to find
    :type pat_len: int
    :param min_times: the minimum number of times the substring must
                      appear (default 2)
    :type min_times: int
    :param all_above: wether to save all substrings which appear at
                      least min_times or only the most frequent ones
    :type all_above: bool
    :returns: all most frequent substrings with rules specified above
    :rtype: list
    """
    
    if min_times < 2:
        raise ValueError('Minimum number of appearances must be >=2')
    freq = calc_freq(DNA, pat_len)
    most_freq = []
    for i in range(len(freq)):
        if freq[i] >= min_times:
            # only clear most_freq & update min_times if allowed to
            if not all_above and freq[i] > min_times:
                most_freq = []
                min_times = freq[i]
            most_freq.append(num_to_DNA(i, pat_len))
    return most_freq

def find_clumps(DNA: str, pat_len: int, min_times: int,
                window_len: int) -> list:
    """Find substrings which appear with enough frequency in a restricted window

    :param DNA: the longer string to search in
    :type DNA: str
    :param pat_len: the length of substrings to find
    :type pat_len: int
    :param min_times: the minimum number of times the substring must
                      appear (default 2)
    :type min_times: int
    :param window_len: the char-length of the maximal window which enough
                          substring instances must be contained in
    :type window_len: int
    :returns: all most frequent substrings with rules specified above
    :rtype: list
    """
    
    if window_len > len(DNA):
        raise ValueError('Window length cannot be longer than DNA length')
    if window_len < pat_len + min_times:
        raise ValueError('Window is not long enough to accomadate the '
                         + 'minimum number of patterns')
    freq_pats = find_most_freq(DNA, pat_len, min_times, True)
    clumped = []
    for pat in freq_pats:
        starts = find_starts(DNA, pat)
        for i in range(len(starts) - min_times + 1):
            # check distance of patterns from each other
            if (starts[i + min_times - 1] - starts[i]
                < window_len - pat_len + 1):
                clumped.append(pat)
                break
    return clumped

def find_approx_starts(DNA: str, pat: str, dist: int) -> list:
    """Find all start indexes of an approximate substring

    :param DNA: the longer string to search in
    :type DNA: str
    :param pattern: the substring to search for
    :type pattern: str
    :param dist: the maximum number of substitutions to allow
    :type dist: int
    :returns: all indexes where pattern approximately starts in DNA
    :rtype: list
    """
    
    if dist == 0:
        return find_starts(DNA, pat)
    if dist < 0:
        raise ValueError('Cannot use a negative distance')
    if not DNA:
        raise ValueError('Cannot search in empty string')
    if not pat:
        raise ValueError('Cannot search for empty string')
    starts = []
    pat_len = len(pat)
    # only search possible starting positions
    for i in range(len(DNA) - len(pat) + 1):
        if ham_dist(DNA[i:i + pat_len], pat) <= dist:
            starts.append(i)
    return starts

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
    # set up lists with the single-base strings
    correct = pat[0]
    try:
        BASES.index(correct)
    except ValueError:
        raise ValueError('Non-DNA base "' + correct + '" in given string')
    for base in BASES:
        if base == correct:
            neighbors[0] = [base]
        else:
            neighbors[1].append(base)
    # add on new chars to each existing string one at a time
    for i in range(1, len(pat)):
        old = list(neighbors)
        neighbors = [[] for _ in range(dist + 1)]
        correct = pat[i]
        try:
            BASES.index(correct)
        except ValueError:
            raise ValueError('Non-DNA base "' + correct + '" in given string')
        # loop over all dist-lists
        for cur_dist in range(0, dist + 1):
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

def find_most_approx_freq(DNA: str, pat_len: int, dist: int,
                          allow_rev_comp: bool = False) -> list:
    """Find substrings which appear with enough frequency, with mismatches

    :param DNA: the longer string to search in
    :type DNA: str
    :param pat_len: the length of substrings to find
    :type pat_len: int
    :param dist: the maximum number of substitutions to allow
    :type dist: int
    :param allow_rev_comp: whether to count reverse complements together
                           (default of False)
    :type allow_rev_comp: bool
    :returns: all most frequent substrings with rules specified above
    :rtype: list
    """
    
    if dist == 0:
        return find_most_freq(DNA, pat_len)
    if pat_len < 1:
        raise ValueError('Patterns must have length at least 1')
    if not DNA:
        raise ValueError('Cannot search in empty string')
    if len(DNA) < pat_len:
        raise ValueError('Cannot search for patterns longer than the string')
    # build up frequency map
    freq = {}
    for i in range(len(DNA) - pat_len + 1):
        pat = DNA[i: i + pat_len]
        # add ALL neighbors up 1
        for neighbor in get_neighbors(pat, dist):
            try:
                freq[neighbor] += 1
            except KeyError:
                freq[neighbor] = 1
    freq_pats = []
    best_freq = 0
    for pat, num in freq.items():
        # add rev-comp if allowed and re-comp is different
        if allow_rev_comp:
            try:
                rev = rev_comp(pat)
                if rev != pat:
                    num += freq[rev_comp(pat)]
            except KeyError:
                pass
        if num > best_freq:
            freq_pats = [pat]
            best_freq = num
        elif num == best_freq:
            freq_pats.append(pat)
    return freq_pats

if __name__ == '__main__':
    # read in parameters
    with open('data.txt') as data:
        DNA = data.readline().rstrip()
        pat_len, dist = [int(x) for x in data.readline().split()]
    # print all results
    for pat in find_most_approx_freq(DNA, pat_len, dist, True):
        print(pat, end=' ')
