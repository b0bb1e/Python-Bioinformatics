def _calc_grid(one: str, two: str, indel_penalty: int,
              v_taxi_start: bool, h_taxi_start: bool,
              v_taxi_end: bool, h_taxi_end: bool,
              score_matrix: dict, match: int, no_match: int) -> list:
    """Calculate optimal values and backtracks for alignment grid

    Either score_matrix should have a value,
    or match and no_match should have values

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :param v_taxi_start: whether vertical taxis from the source are allowed
    :type v_taxi_start: bool
    :param h_taxi_start: whether horizontal taxis from the source are allowed
    :type h_taxi_start: bool
    :param v_taxi_end: whether vertical taxis to the sink are allowed
    :type v_taxi_end: bool
    :param h_taxi_end: whether horizonal taxis to the sink are allowed
    :type h_taxi_end: bool
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param match: the score bonus for matching chars in alignment
    :type match: int (positive)
    :param no_match: the socre deduction for non-matching chars in alignment
    :type no_match: int (negative)
    :returns: the alignment grid
    :rtype: list (of lists (of tuples (int, str)))
    """
    
    one_len, two_len = len(one), len(two)
    # set up the first column
    if v_taxi_start:
        # if taxis are allowed, then all optimally backtrack to source
        grid = [[(0, 's')] for i in range(one_len + 1)]
    else:
        grid = [[(indel_penalty * i, 'v')] for i in range(one_len + 1)]
    # set up the first row
    if h_taxi_start:
        # if taxis are allowed, then all optimally backtrack to source
        grid[0] = [(0, 's') for i in range(two_len + 1)]
    else:
        grid[0] = [(indel_penalty * i, 'h') for i in range(two_len + 1)]
    # if taxis to sink are allowed, set up trackers for where to taxi from
    if v_taxi_end or h_taxi_end:
        max_val, max_row, max_col = 0, 0, 0
    can_taxi = v_taxi_start or h_taxi_start or v_taxi_end or h_taxi_end
    # set up the source node
    grid[0][0]= (0, None)
    
    for one_i in range(one_len):
        for two_i in range(two_len):
            # assume horizontal move is the best
            best_val = grid[one_i + 1][two_i][0] + indel_penalty
            backtrack = 'h'

            # try vertical move
            vert_val = grid[one_i][two_i + 1][0] + indel_penalty
            if vert_val > best_val:
                best_val, backtrack = vert_val, 'v'

            # try diagonal move
            diag_val = grid[one_i][two_i][0]
            # use path weight of score matrix or matchness, as appropriate
            if score_matrix:
                diag_val += score_matrix[one[one_i]][two[two_i]]
            else:
                if one[one_i] == two[two_i]:
                    diag_val += match
                else:
                    diag_val += no_match
            if diag_val > best_val:
                best_val, backtrack = diag_val, 'd'

            if can_taxi:
                # change backtrack to start if better
                if best_val < 0 and (v_taxi_start and h_taxi_start):
                    best_val, backtrack = 0, 's'
                # if can taxi to sink, and is best so far, update trackers
                if (best_val > max_val
                    and ((v_taxi_end and h_taxi_end)
                         or (v_taxi_end and two_i == two_len - 1)
                         or (h_taxi_end and one_i == one_len - 1))):
                    max_val, max_row, max_col = best_val, one_i + 1, two_i + 1
            grid[one_i + 1].append((best_val, backtrack))
    # taxi to end if allowed and better
    if (v_taxi_end or h_taxi_end) and grid[one_len][two_len][0] < max_val:
        # special backtrack to exact cell
        grid[one_len][two_len] = (max_val, (max_row, max_col))
    return grid

def _backtrack_alignment(one: str, two: str, grid: list) -> str:
    """Backtrack through alignment grid to determine optimal alignment

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param grid: the alignment grid
    :type grid: list (of lists (of tuples (int, str)))
    :returns: the optimal alignment, with one + newline + two
    :rtype: str
    """
    
    one_align, two_align = '', ''
    one_len, two_len = len(one), len(two)
    # if sink has special taxi backtrack, start at its backtrack
    if len(grid[one_len][two_len][1]) == 2:
        cur_row, cur_col = grid[one_len][two_len][1]
    # otherwise start at sink
    else:
        cur_row, cur_col = one_len, two_len
    
    # while not yet backtracked to source
    while not (cur_row == 0 and cur_col == 0):
        backtrack = grid[cur_row][cur_col][1]
        # vertical backtrack
        if backtrack == 'v':
            cur_row -= 1
            one_align = one[cur_row] + one_align
            two_align = '-' + two_align
        # horizontal backtrack
        elif backtrack == 'h':
            cur_col -= 1
            one_align = '-' + one_align
            two_align = two[cur_col] + two_align
        # diagonal backtrack
        elif backtrack == 'd':
            cur_row -= 1
            cur_col -= 1
            one_align = one[cur_row] + one_align
            two_align = two[cur_col] + two_align
        # backtrack to source
        elif backtrack == 's':
            cur_row = 0
            cur_col = 0
    return one_align + '\n' + two_align

def _align(one: str, two: str, indel_penalty: int,
           v_taxi_start: bool, h_taxi_start: bool,
           v_taxi_end: bool, h_taxi_end: bool, score_matrix: dict=None,
           match: int=None, no_match: int=None) -> (int, str):
    """Optimally align two strings

    Either score_matrix should have a value,
    or match and no_match should have values

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :param v_taxi_start: whether vertical taxis from the source are allowed
    :type v_taxi_start: bool
    :param h_taxi_start: whether horizontal taxis from the source are allowed
    :type h_taxi_start: bool
    :param v_taxi_end: whether vertical taxis to the sink are allowed
    :type v_taxi_end: bool
    :param h_taxi_end: whether horizonal taxis to the sink are allowed
    :type h_taxi_end: bool
    -keyword params below this point-
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param match: the score bonus for matching chars in alignment
    :type match: int (positive)
    :param no_match: the socre deduction for non-matching chars in alignment
    :type no_match: int (negative)
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    grid = _calc_grid(one, two, indel_penalty, v_taxi_start, h_taxi_start,
                      v_taxi_end, h_taxi_end, score_matrix,
                      match, no_match)
    score = grid[len(one)][len(two)][0]
    return score, _backtrack_alignment(one, two, grid)

def global_align(one: str, two: str, score_matrix: dict) -> (int, str):
    """Globally align two strings

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    return _align(one, two, -5, False, False, False, False, score_matrix)

def local_align(one: str, two: str, score_matrix: dict) -> (int, str):
    """Locally align two strings

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    return _align(one, two, -5, True, True, True, True, score_matrix)

def find_edit_distance(one: str, two: str) -> int:
    """Find the edit distance between two strings

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :returns: the optimal score and alignment
    :rtype: int
    """
    
    grid = _calc_grid(one, two, -1, False, False, False, False, None,
                     match=0, no_match=-1)
    return -grid[len(one)][len(two)][0]

def fitting_align(long: str, short: str) -> (int, str):
    """Find the fitting alignment of two strings

    :param long: the longer string to use a section of
    :type long: str
    :param short: the shorter string to use the whole of
    :type short: str
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    return _align(long, short, -1, True, False, True, False,
                 match=1, no_match=-1)

def overlap_align(before: str, after: str) -> (int, str):
    """Find the overlap alignment of two strings

    :param before: the string to use a suffix of
    :type before: str
    :param after: the string to use a prefix of
    :type after: str
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    return _align(before, after, -2, True, False, False, True,
                 match=1, no_match=-2)

def read_score_matrix(file_name: str) -> dict:
    """Read a scoring matrix from a file

    :param file_name: the file with the matrix
    :type fil_name: str
    :returns: a scoring matrix for proteins
    :rtype: dict (strs : dicts (strs: ints))
    """
    
    with open(file_name) as score_file:
        acids = score_file.readline().split()
        scores = {acid: {} for acid in acids}
        for line in score_file:
            line_scores = [int(x) for x in line[1:].split()]
            scores[line[0]] = dict(zip(acids, line_scores))
    return scores

if __name__ == '__main__':
    with open('data.txt') as data:
        before = data.readline().rstrip()
        after = data.readline().rstrip()
    print(*overlap_align(before, after), sep='\n')
