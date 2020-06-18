def comp(DNA: str, pat_len: int) -> list:
    """Sort all substrings of pat_len length

    :param DNA: the string to pull substrings from
    :type DNA: str
    :param pat_len: the length of substrings to pull
    :type pat_len: int
    :returns: all substrings, sorted
    :rtype: list (of strs)
    """

    if not DNA:
        raise ValueError('Cannot pull substrings from empty string')
    if pat_len < 1:
        raise ValueError('Substrings must be at least length 1')
    DNA_len = len(DNA)
    if DNA_len < pat_len:
        raise ValueError('No substrings of that length')

    return sorted([DNA[i:i + pat_len] for i in range(DNA_len - pat_len + 1)])

if __name__ == '__main__':
    with open('data.txt') as data:
        pat_len = int(data.readline().rstrip())
        DNA = data.readline().rstrip()
    with open('output.txt', mode='w') as output:
        for sub_str in comp(DNA, pat_len):
            output.write(sub_str + '\n')
