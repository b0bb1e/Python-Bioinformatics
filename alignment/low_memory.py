class Low_Memory_Aligner:
    """An aligner which can align two letter strings optimally in linear memory

    A score matrix which has all possible match-ups of letters must be given.

    find_midddle_edge: find the middle edge in the current alignment graph
    """
    
    def __init__(self, one, two, score_matrix, indel_penalty):
        """Initialize all variables needed to construct the alignment graph

        :param one: the string along the side of the grid
        :type one: str
        :param two: the string along the top of the grid
        :type two: str
        :param score_matrix: a scoring matrix for proteins
        :type score_matrix: dict (strs : dicts (strs: ints))
        :param indel_penalty: the penalty for using an indel in an alignment
        :type indel_penalty: int
        """

        self.one = one
        self.two = two
        self.s_m= score_matrix
        self.i_p = indel_penalty

        # the penalty can't be positive
        if self.i_p > 0:
            self.i_p = -indel_pentaly

    def find_middle_edge(self) -> ((int, int), str):
        """Find the middle edge in the alignment graph

        Uses linear space and quadratic time

        :returns: the middle point & its optimal path forward
        :rtype: tuple (tuple (int, int), str)
        """

        return self._find_middle_edge(0, 0, len(self.one), len(self.two))

    def _find_middle_edge(self, min_row: int, min_col: int, max_row: int,
                          max_col: int) -> ((int, int), str):
        """Find the middle edge in a certain area of the alignment graph

        Uses linear space and quadratic time

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
                  in the format ((row, col), direction)
        :rtype: tuple (tuple (int, int), str)
        """

        if min_col == max_col:
            return ((min_row, min_col), 'v')

        mid_col = int((min_col + max_col) / 2)
        # sweep from left edge to middle column
        mid_from_left = self._calc_col_from_left(min_row, min_col,
                                                 max_row, mid_col)
        # sweep from right edge to column one right of middle column
        right_of_mid = self._calc_col_from_right(min_row, mid_col + 1,
                                                 max_row, max_col)
        
        backtracks = ['h']
        # optimal values for column one right of middle column from left
        one_past_mid = [mid_from_left[0] + self.i_p]
        for row in range(min_row + 1, max_row + 1):
            cur_index = row - min_row
            best_val = one_past_mid[cur_index - 1] + self.i_p
            best_backtrack = 'v'

            horiz_val = mid_from_left[cur_index] + self.i_p
            if horiz_val >= best_val:
                best_val, best_backtrack = horiz_val, 'h'

            diag_val = (mid_from_left[cur_index - 1]
                        + self.s_m[self.one[row - 1]][self.two[mid_col]])
            if diag_val >= best_val:
                best_val, best_backtrack = diag_val, 'd'

            one_past_mid.append(best_val)
            backtracks.append(best_backtrack)

        max_val = float('-inf')
        max_row = 0
        # find optimal node to go to in one right of middle column
        for row in range(len(one_past_mid)):
            cur_val = one_past_mid[row] + right_of_mid[row]
            if cur_val > max_val:
                max_val, max_row = cur_val, row

        # figure out coordinates of middle node
        if backtracks[max_row] == 'd':
            mid_row = max_row - 1
        else:
            mid_row = max_row
        return ((mid_row, mid_col), backtracks[max_row])

    def _calc_col_from_left(self, top_row: int, left_col: int,
                            bottom_row: int, right_col: int) -> list:
        """Calculate from left values for a column in the alignment grid

        :param top_row: the top row of the sweep
        :type top_row: int
        :param left_col: the left column of the sweep
        :type left_col: int
        :param bottom_row: the bottom row of the sweep
        :type bottom_row: int
        :param right_col: the right column of the sweep
        :type right_col: int
        :returns: the optimal values for the column at the end of the sweep
        :rtype: list (of ints)
        """

        # set up initial column (all indels)
        num_rows = (bottom_row - top_row) + 1
        last_col = [self.i_p * i for i in range(num_rows)]
        cur_col = list(last_col)
        
        # sweep from one right of start column to end column
        for col in range(left_col + 1, right_col + 1):
            # top cell must be an indel
            cur_col[0] = last_col[0] + self.i_p
            # sweep from one after top row to bottom row
            for row in range(top_row + 1, bottom_row + 1):
                cur_i = row - top_row
                # maximize this cell's value
                cur_col[cur_i] = (
                    max(last_col[cur_i] + self.i_p,
                        cur_col[cur_i - 1] + self.i_p,
                        (last_col[cur_i - 1]
                         + self.s_m[self.one[row - 1]][self.two[col - 1]])))
                
        return cur_col


    def _calc_col_from_right(self, top_row: int, left_col: int,
                             bottom_row: int, right_col: int) -> list:
        """Calculate from right values for a column in the alignment grid

        :param top_row: the top row of the sweep
        :type top_row: int
        :param left_col: the left column of the sweep
        :type left_col: int
        :param bottom_row: the bottom row of the sweep
        :type bottom_row: int
        :param right_col: the right column of the sweep
        :type right_col: int
        :returns: the optimal values for the column at the end of the sweep
        :rtype: list (of ints)
        """

        # set up initial column (all indels)
        num_rows = (bottom_row - top_row) + 1
        last_col = [self.i_p * i for i in range(num_rows)]
        last_col.reverse()
        cur_col = list(last_col)
        
        # sweep from one left of start column to end column
        for col in range(right_col - 1, left_col - 1, -1):
            # bottom cell must be an indel
            cur_col[num_rows - 1] = last_col[num_rows - 1] + self.i_p
            # sweep from on above bottom row to top row
            for row in range(bottom_row - 1, top_row - 1, -1):
                cur_i = top_row - row
                # maximize this cell's value
                cur_col[cur_i] = max(last_col[cur_i] + self.i_p,
                                     cur_col[cur_i + 1] + self.i_p,
                                     (last_col[cur_i + 1]
                                      + self.s_m[self.one[row]][self.two[col]]))        

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
    lm_aligner = Low_Memory_Aligner(one, two, score_matrix, -5)
    print(lm_aligner.find_middle_edge())
