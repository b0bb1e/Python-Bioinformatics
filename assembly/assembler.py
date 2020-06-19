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

def path_to_DNA(path: list) -> str:
    """Convert a path of substrings to a condensed string

    :param path: the ordered path of substrings
    :type path: list (of strs)
    :returns: the overall string
    :rtype: str
    """

    if not path:
        raise ValueError('Cannot convert nonexistant path')
    sub_len = len(path[0])
    for sub_str in path:
        if len(sub_str) != sub_len:
            raise ValueError('All substrings in path must be the same length')
    DNA = path[0][:-1]
    for sub_str in path:
        DNA += sub_str[-1]
    return DNA

if __name__ == '__main__':
    with open('data.txt') as data:
        path = []
        for line in data:
            path.append(line.rstrip())
    print(path_to_DNA(path))
