def _calc_col_from_edge_(one: str, two: str, score_matrix: dict,
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

def _calc_col_from_left(one: str, two: str, score_matrix: dict,
                        indel_penalty: int, start_row: int, start_col: int,
                        end_row: int, end_col: int) -> list:
    print("going left col #", start_col, " -> col#", end_col)

    # set up initial column (all indels)
    num_rows = end_row - start_row + 1
    last_col = [indel_penalty * i for i in range(num_rows)]
    cur_col = list(last_col)
    # sweep from one right of start column to end column
    for col in range(start_col + 1, end_col + 1):
        print(cur_col)
        # top cell must be an indel
        cur_col[0] = last_col[0] + indel_penalty
        # sweep from one after top row to bottom row
        for row in range(start_row + 1, end_row + 1):
            cur_index = row - start_row
            # maximize this cell's value
            cur_col[cur_index] = max(last_col[cur_index] + indel_penalty,
                                     cur_col[cur_index - 1] + indel_penalty,
                                     (last_col[cur_index - 1] +
                                      score_matrix[one[row - 1]][two[col - 1]]
                                      ))
    print(cur_col)
    return cur_col


def _calc_col_from_right(one: str, two: str, score_matrix: dict,
                         indel_penalty: int, start_row: int, start_col: int,
                         end_row: int, end_col: int) -> list:
    print("going right col #", start_col, " -> col#", end_col)

    # set up initial column (all indels)
    num_rows = (start_row - end_row) + 1
    last_col = [indel_penalty * i for i in range(num_rows)]
    last_col.reverse()
    cur_col = list(last_col)
    # sweep from one left of start column to end column
    for col in range(start_col - 1, end_col - 1, -1):
        print(cur_col)
        # bottom cell must be an indel
        cur_col[num_rows - 1] = last_col[num_rows - 1] + indel_penalty
        # sweep from bottom row to one before top row
        for row in range(start_row - 1, end_row - 1, -1):
            cur_index = end_row - row
            # maximize this cell's value
            cur_col[cur_index] = max(last_col[cur_index] + indel_penalty,
                                     cur_col[cur_index + 1] + indel_penalty,
                                     (last_col[cur_index + 1]
                                      + score_matrix[one[row]][two[col]]
                                      ))        
    print(cur_col)
    return cur_col
