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

def verify_pats(pats: list):
    """Verify that a group of pattens is valid; raises errors

    :param pats: a group of patterns
    :type pats: list (of strs)
    :raises: ValueError (if the patterns are invalid)
    """

    if not pats:
        raise ValueError('Cannot convert nonexistant pattern')
    pat_len = len(pats[0])
    for pat in pats:
        if len(pat) != pat_len:
            raise ValueError('All patterns must be the same length')

def path_to_DNA(path: list) -> str:
    """Convert a path of substrings to a condensed string

    :param path: the ordered path of substrings
    :type path: list (of strs)
    :returns: the overall string
    :rtype: str
    """

    verify_pats(path)
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

    verify_pats(pats)
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

def self_overlap_graph(pats: list) -> dict:
    """Construct a self-overlap graph for a group of patterns

    The tailess bit of a sting self-overlaps with the headless bit
    
    :param pats: a group of substrings
    :type pats: list (of strs)
    :returns: the overlap graph in dictionary form
    :rtype: dict (strs: lists (of strs))
    """

    verify_pats(pats)
    graph = {}
    for pat in pats:
        start = pat[:-1]
        if not start in graph:
            graph[start] = []
        graph[start].append(pat[1:])
    return graph

def new_cycle(graph: dict, old_cycle: list=None) -> list:
    """Find a new cycle within a graph

    :param graph: a self-overlap graph with left-over edges
    :type graph: dict (strs: lists (of strs))
    :param old_cycle: a previously-constructed cycle (default None)
    :type old_cycle: list (of strs)
    :returns: a new cycle guarenteed to be longer than old_cycle
    :rtype: list (of strs)
    """
    
    if old_cycle:
        # find a node with unused edges
        for node in range(len(old_cycle)):
            if old_cycle[node] in graph:
                end = node
                break
        # rotate cycle to start & end with this node
        try:
            cycle = old_cycle[end:] + old_cycle[1:end + 1]
            cur = cycle[0]
        except NameError:
            raise ValueError('Graph has a closed cycle')
    else:
        cur = next(iter(graph.keys()))
        cycle = []

    # while there are unused edges from the current node
    while cur in graph:
        # use an edge
        cycle.append(graph[cur].pop(-1))
        # delete nodes from graph if no edges are left
        if not graph[cur]:
            del graph[cur]
        cur = cycle[-1]
    # force a cycle if there wasn't one before
    if not old_cycle:
        cycle.append(cycle[0])
    return cycle

def graph_to_cycle(graph: dict) -> list:
    """Constructs a Eulerian cycle from a self-overlap graph

    :param graph: a self-overlap graph
    :type graph: dict (strs: lists (of strs))
    :returns: a cycle using all paths in graph
    :rtype: list
    """

    cycle = new_cycle(graph)
    while graph:
        cycle = new_cycle(graph, cycle)
    return cycle

def graph_to_path(graph: dict) -> list:
    """Constructs a Eulerian path from a self-overlap graph

    :param graph: a self-overlap graph
    :type graph: dict (strs: lists (of strs))
    :returns: a path using all paths in graph
    :rtype: list
    """

    # used to find the start (1 extra out) and end (1 extra in) nodes
    net_outs = {}
    for start in graph:
        if start in net_outs:
            net_outs[start] += len(graph[start])
        else:
            net_outs[start] = len(graph[start])
        for out in graph[start]:
            if out in net_outs:
                net_outs[out] -= 1
            else:
                net_outs[out] = -1
    for node in net_outs:
        if net_outs[node] == -1:
            end = node
        elif net_outs[node] == 1:
            start = node
        elif net_outs[node] != 0:
            raise ValueError('Given graph is not a Eulerian path')
    # connect end to start to make a cycle
    if end in graph:
        graph[end].append(start)
    else:
        graph[end] = [start]
    cycle = graph_to_cycle(graph)
    start_indexes = []
    for node in range(len(cycle)):
        if cycle[node] == start:
            start_indexes.append(node)
    # rotate cycle so that start & end appear in proper spots
    for index in start_indexes:
        path = cycle[index:] + cycle[1:index + 1]
        if path[-2] == end:
            # get rid of the edge between end and start
            del path[-1]
            break
    return path

def assemble(pats: list) -> str:
    """Assemble a string from substrings

    :param pats: a group of substrings
    :type pats: list (of strs)
    :returns: the assembled string
    :rtype: str
    """

    verify_pats(pats)
    graph = self_overlap_graph(pats)
    path = graph_to_path(graph)
    return path_to_DNA(path)

def all_binary(bin_len: int) -> list:
    """Find all binary strings of a certain length

    :param bin_len: the length the binary strings should be
    :type bin_len: int
    :returns: all the possible binarys strings of this length
    :rtype: list (of strs)
    """
    
    if bin_len == 0:
        return ['']
    elif bin_len < 0:
        raise ValueError('Binary strings cannot have negative length')
    all_bin = []
    # extend all shorter strings by 0 and 1
    for short_bin in all_binary(bin_len - 1):
        all_bin += [short_bin + '0', short_bin + '1']
    return all_bin

def universal_binary(bin_len: int) -> str:
    """Constuct a ciruclar universal binary string

    :param bin_len: the length of the universal binary patterns
    :type bin_len: int
    :returns: a string containing all binary strings of legnth bin_len
    :rtype: str
    """
    
    graph = self_overlap_graph(all_binary(bin_len))
    cycle = graph_to_cycle(graph)
    # delete last bin_len - 1 chars since this is CIRCULAR
    cycle = cycle[:-(bin_len - 1)]
    return path_to_DNA(cycle)

if __name__ == '__main__':
    print(universal_binary(9))
