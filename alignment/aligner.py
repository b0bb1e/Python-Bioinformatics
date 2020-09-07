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
            # current cell is at grid[one_i + 1][two_i + 1]
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
    :returns: the edit distance
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

def _calc_affine_grids(one: str, two: str, gap_open: int, gap_ext: int,
                       score_matrix: dict) -> (list, list, list):
    """Calcualte the three levels of an affine alignment grid

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param gap_open: the penalty for the first indel in a gap
    :type gap_open: int (negative)
    :param gap_ext: the penalty for each indel after the first
    :type gap_ext: int (negative, but less so than gap_open)
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :returns: the optimal score and alignment
    :rtype: tuple (list(of lists (of tuples (int, str)))), same, same)
    """
    
    one_len, two_len = len(one), len(two)
    # set up left column (horiz doesn't have)
    vert = [[(gap_open + (gap_ext * (i - 1)), 'v')]
               for i in range(1, one_len + 1)]
    vert[0][0] = (gap_open, 'd')
    diag = [[(gap_open + (gap_ext * (i - 1)), 'v')]
            for i in range(one_len + 1)]
    # set up horiz's rows
    horiz = [[] for i in range(one_len + 1)]
    # set up top row (vert doesn't have)
    horiz[0] = [(gap_open + (gap_ext * (i - 1)), 'h')
                 for i in range(1, two_len + 1)]
    horiz[0][0] = (gap_open, 'd')
    diag[0] = [(gap_open + (gap_ext * (i - 1)), 'h')
                for i in range(two_len + 1)]
    # re-set source values
    diag[0][0] = (0, None)
    
    for one_i in range(one_len):
        for two_i in range(two_len):
            # current cell to calculate is at vert[one_i][two_i + 1],
            # horiz[one_i + 1][two_i], and diag[one_i + 1][two_i + 1]

            # assume vert backtracks to diag
            vert[one_i].append((diag[one_i][two_i + 1][0] + gap_open, 'd'))
            if one_i > 0:
                if_ext = vert[one_i - 1][two_i + 1][0] + gap_ext
                # chagne to backtrack to vert if better
                if if_ext > vert[one_i][two_i + 1][0]:
                    vert[one_i][two_i + 1] = (if_ext, 'v')

            # similar logic for horiz
            horiz[one_i + 1].append((diag[one_i + 1][two_i][0] + gap_open, 'd'))
            if two_i > 0:
                if_ext = horiz[one_i + 1][two_i - 1][0] + gap_ext
                if if_ext > horiz[one_i + 1][two_i][0]:
                    horiz[one_i + 1][two_i] = (if_ext, 'h')

            # assume diag backtracks to vert
            diag[one_i + 1].append((vert[one_i][two_i + 1][0], 'v'))
            # change to backtrack to horiz if better
            if horiz[one_i + 1][two_i][0] > diag[one_i + 1][two_i + 1][0]:
                diag[one_i + 1][two_i + 1] = (horiz[one_i + 1][two_i][0], 'h')
            if_diag = (diag[one_i][two_i][0]
                       + score_matrix[one[one_i]][two[two_i]])
            # change to backtrack to diag if better
            if if_diag > diag[one_i + 1][two_i + 1][0]:
                diag[one_i + 1][two_i + 1] = (if_diag, 'd')
                
    return vert, horiz, diag

def _backtrack_affine_alignment(one: str, two: str, vert: list, horiz: list,
                                diag: list) -> str:
    """Backtrack through and affine alignment grid
        to determine optimal alignment

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param vert: the vertical-moves (indels for two) grid
    :type grid: list (of lists (of tuples (int, str)))
    :param horiz: the horizonal-moves (indels for one) grid
    :type grid: list (of lists (of tuples (int, str)))
    :param horiz: the diagonal-moves (no indels) grid
    :type grid: list (of lists (of tuples (int, str)))
    :returns: the optimal alignment, with one + newline + two
    :rtype: str
    """

    # start at sink in diagonal level
    cur_row, cur_col = len(one), len(two)
    level = 'd'
    one_align, two_align = '', ''
    
    # backtrack until hit the source node at (0, 0)
    while not (cur_row == 0 and cur_col == 0):
        if level == 'd':
            backtrack = diag[cur_row][cur_col][1]
            # only d-d moves will update position and alignment
            if backtrack == 'd':
                cur_row -= 1
                cur_col -= 1
                one_align = one[cur_row] + one_align
                two_align = two[cur_col] + two_align
        elif level == 'v':
            backtrack = vert[cur_row - 1][cur_col][1]
            # all v- moves go back one row
            cur_row -= 1
            one_align = one[cur_row - 1] + one_align
            two_align = '-' + two_align
        elif level == 'h':
            backtrack = horiz[cur_row][cur_col - 1][1]
            # all h- mvoes go back one row
            cur_col -= 1
            one_align = '-' + one_align
            two_align = two[cur_col - 1] + two_align
        level = backtrack
    return one_align + '\n' + two_align

def affine_align(one: str, two: str, gap_open: int, gap_ext: int,
                 score_matrix: dict) -> (int, str):
    """Optimally align two strings with affine gap penalities

    Affine gap penalties means the opening of a gap is penalized much
    more than the extenion of a gap

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param gap_open: the penalty for the first indel in a gap
    :type gap_open: int (negative)
    :param gap_ext: the penalty for each indel after the first
    :type gap_ext: int (negative, but less so than gap_open)
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """

    vert, horiz, diag = _calc_affine_grids(one, two, gap_open, gap_ext,
                                           score_matrix)
    score = diag[len(one)][len(two)][0]
    return score, _backtrack_affine_alignment(one, two, vert, horiz, diag)

def find_middle_edge(one: str, two: str, score_matrix: dict,
                     indel_penalty: int) -> ((int, int), str):
    """Find the middle edge in an alignment graph

    Uses linear space and quadratic time

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :returns: the middle point & its optimal path forward
    :rtype: tuple (tuple (int, int), str)
    """

    return _find_middle_edge(one, two, score_matrix, -5,
                             0, 0, len(one), len(two))

def _find_middle_edge(one: str, two: str, score_matrix: dict,
                      indel_penalty: int, min_row: int, min_col: int,
                      max_row: int, max_col: int) -> ((int, int), str):
    """Find the middle edge in a certain area of an alignment graph

    Uses linear space and quadratic time

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :param min_row: the lowest row number (top=0) of the area to search
    :type min_row: int
    :param min_col: the lowest column number (left=0) of the area to search
    :type min_col: int
    :param max_row: the highest row number (top=0) of the area to search
    :type max_row: int
    :param max_col: the highest column number (top=0) of the area to search
    :type max_col: int
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :returns: the middle point & its optimal path forward
    :rtype: tuple (tuple (int, int), str)
    """

    middle = int((min_col + max_col) / 2) + 1
    left_of_middle = _calc_col_from_edge(one, two, score_matrix, -5, True,
                                         min_row, min_col, max_row, middle - 1)
    middle_from_right = _calc_col_from_edge(one, two, score_matrix, -5, False,
                                            max_row, max_col, min_row, middle)
    backtracks = ['h']
    middle_from_left = [left_of_middle[0] + indel_penalty]
    for row in range(min_row + 1, max_row + 1):
        cur_index = row - min_row
        best_val = middle_from_left[cur_index - 1] + indel_penalty
        best_backtrack = 'v'

        horiz_val = left_of_middle[cur_index] + indel_penalty
        if horiz_val > best_val:
            best_val, best_backtrack = horiz_val, 'h'

        diag_val = (left_of_middle[cur_index - 1]
                    + score_matrix[one[row - 1]][two[middle - 1]])
        if diag_val > best_val:
            best_val, best_backtrack = diag_val, 'd'

        middle_from_left.append(best_val)
        backtracks.append(best_backtrack)

    max_val = float('-inf')
    max_row = 0
    for row in range(len(middle_from_left)):
        cur_val = middle_from_left[row] + middle_from_right[row]
        if cur_val > max_val:
            max_val, max_row = cur_val, row
    
    if backtracks[max_row] != 'v':
        mid_row = max_row - 1
    else:
        mid_row = max_row
    return ((mid_row, middle - 1), backtracks[max_row])

def _calc_col_from_edge(one: str, two: str, score_matrix: dict,
                        indel_penalty: int, going_right: bool, start_row: int,
                        start_col: int, end_row: int, end_col: int) -> list:
    """Calculate the optimal values for a column in the alignment grid

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :param going_right: whether the calculator should sweep left-to-right
    :type going_right: bool
    :param start_row: the starting row of the sweep
    :type start_row: int
    :param start_col: the starting column of the sweep
    :type start_col: int
    :param end_row: the ending row of the sweep
    :type end_row: int
    :param end_col: the ending column of the sweep
    :type end_col: int
    :returns: the optimal values for the column at the end of the sweep
    :rtype: list (of ints)
    """

    if going_right:
        change = 1
    else:
        change = -1

    cur_col = [indel_penalty * i for i in range(len(one) + 1)]
    for col in range(start_col + change, end_col + change, change):
        for row in range(end_row, start_row, change):
            cur_index = (end_row - row) * change
            cur_col[cur_index] = max(cur_col[cur_index] + indel_penalty,
                                     cur_col[cur_index - 1] + indel_penalty,
                                     (cur_col[cur_index - 1]
                                      + score_matrix[one[row - 1]][two[col - 1]]
                                      ))
        cur_col[0] = cur_col[0] + indel_penalty
    if not going_right:
        cur_col.reverse()
    return cur_col
    
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
        one = data.readline().rstrip()
        two = data.readline().rstrip()
    score_matrix = read_score_matrix('blossom.txt')
    edge = find_middle_edge(one, two, score_matrix, -5)
