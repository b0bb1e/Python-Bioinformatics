def calc_grid(one: str, two: str, score_matrix: dict,
              indel_penalty: int) -> list:
    """Calculate optimal values and backtracks for alignment grid

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :returns: the alignment grid, where the i-th row and j-th column holds
              the optimal value after using i letters from one and j from two
              in [0] and the optimal backtrack ('h'orizontal, 'v'ertical,
              'd'iagonal) in [1]
    :rtype: list (of lists (of tuples (int, str)))
    """
    
    one_len = len(one)
    two_len = len(two)
    # set up the first column
    grid = [[(indel_penalty * i, (i - 1, 0))] for i in range(one_len + 1)]
    # set up the top row
    grid[0] = [(indel_penalty * i, (0, i - 1)) for i in range(two_len + 1)]
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
            diag_val = (grid[one_i][two_i][0]
                        + score_matrix[one[one_i]][two[two_i]])
            if diag_val > best_val:
                best_val, backtrack = diag_val, 'd'
            grid[one_i + 1].append((best_val, backtrack))
    return grid

def backtrack_alignment(one: str, two: str, grid: list) -> str:
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
    
    one_align = ''
    two_align = ''
    cur_row = len(one)
    cur_col = len(two)
    
    while cur_row != 0 or cur_col != 0:
        backtrack = grid[cur_row][cur_col][1]
        if backtrack == 'v':
            cur_row -= 1
            one_align = one[cur_row] + one_align
            two_align = '-' + two_align
        elif backtrack == 'h':
            cur_col -= 1
            one_align = '-' + one_align
            two_align = two[cur_col] + two_align
        elif backtrack == 'd':
            cur_row -= 1
            cur_col -= 1
            one_align = one[cur_row] + one_align
            two_align = two[cur_col] + two_align
    return one_align + '\n' + two_align

def align(one: str, two: str, score_matrix: dict,
          indel_penalty: int) -> (int, str):
    """Optimally align two strings

    :param one: the string along the side of the grid
    :type one: str
    :param two: the string along the top of the grid
    :type two: str
    :param score_matrix: a scoring matrix for proteins
    :type score_matrix: dict (strs : dicts (strs: ints))
    :param indel_penalty: the score deduction for using an indel in alignment
    :type indel_penalty: int (negative)
    :returns: the optimal score and alignment
    :rtype: tuple (int, str)
    """
    
    grid = calc_grid(one, two, score_matrix, indel_penalty)
    score = grid[len(one)][len(two)][0]
    return score, backtrack_alignment(one, two, grid)

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
    print(align(one, two, score_matrix, -5), sep='\n')
