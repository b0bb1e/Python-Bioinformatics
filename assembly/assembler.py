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
    pat_len = len(path[0])
    for pat in path:
        if len(pat) != pat_len:
            raise ValueError('All substrings in path must be the same length')
    DNA = path[0][:-1]
    for pat in path:
        DNA += pat[-1]
    return DNA

def overlap_graph(pats: list) -> dict:
    """Constructs an overlap graph for a group of patterns

    Patterns directionally overlap if one's headless string matches
    the other's tailless one

    :param pats: a group of substrings
    :type pats: list (of strs)
    :returns: the overlap graph in dictionary form
    :rtype: dict (strs: lists (of strs))
    """

    if not pats:
        raise ValueError('Cannot convert nonexistant pattern')
    pat_len = len(pats[0])
    for pat in pats:
        if len(pat) != pat_len:
            raise ValueError('All patterns must be the same length')
    num_pats = len(pats)

    overlap = {}
    for start in range(num_pats):
        start_pat = pats[start]
        for end in range(num_pats):
            if start != end and start_pat[1:] == pats[end][:-1]:
                if not start_pat in overlap:
                    overlap[start_pat] = []
                overlap[start_pat].append(pats[end])
    return overlap

if __name__ == '__main__':
    with open('data.txt') as data:
        pats = []
        for line in data:
            pats.append(line.rstrip())
    overlap = overlap_graph(pats)
    with open('output.txt', mode='w') as output:
        for start in overlap:
            for end in overlap[start]:
                output.write('{0} -> {1}\n'.format(start, end))
