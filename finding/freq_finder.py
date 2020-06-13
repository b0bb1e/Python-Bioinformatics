BASES = ('A', 'C', 'G', 'T')

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
    
    num = 0
    chars = len(DNA)
    for i in range(chars):
        try:
            # earlier chars are given more weight/value in the number
            num += BASES.index(DNA[i]) * (4 ** (chars - i - 1))
        except ValueError:
            print('Invalid base "' + DNA[i] + '" is ignored')
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
    
    DNA = ''
    for i in range(chars):
        # the smaller numbers, pulled off first, go in the back
        DNA = BASES[int(num % 4)] + DNA
        num = num // 4
    if num != 0:
        raise ValueError('Number is too large for length given')
    return DNA

def find_starts(DNA: str, pattern: str) -> list:
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
    if not pattern:
        raise ValueError('Cannot search for empty string')
    pat_len = len(pattern)
    DNA_len = len(DNA)
    if pat_len > DNA_len:
        return []
    starts = []
    # only search possible starting positions
    for i in range(DNA_len - pat_len + 1):
        if DNA[i:i + pat_len] == pattern:
            starts.append(i)
    return starts

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

def find_freq(DNA: str, pat_len: int, min_times: int = 2,
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
                window_length: int) -> list:
    """Find substrings which appear with enough frequency in a restricted window

    :param DNA: the longer string to search in
    :type DNA: str
    :param pat_len: the length of substrings to find
    :type pat_len: int
    :param min_times: the minimum number of times the substring must
                      appear (default 2)
    :type min_times: int
    :param window_length: the char-length of the maximal window which enough
                          substring instances must be contained in
    :type window_length: int
    :returns: all most frequent substrings with rules specified above
    :rtype: list
    """
    
    if window_length > len(DNA):
        raise ValueError('Window length cannot be longer than DNA length')
    if window_length < pat_len + min_times:
        raise ValueError('Window is not long enough to accomadate the '
                         + 'minimum number of patterns')
    freq_pats = find_freq(DNA, pat_len, min_times, True)
    clumped = []
    for pat in freq_pats:
        starts = find_starts(DNA, pat)
        for i in range(len(starts) - min_times + 1):
            # check distance of patterns from each other
            if starts[i + min_times - 1] - starts[i] < window_length:
                clumped.append(pat)
                break
    return clumped
