def calc_grid(one: str, two: str, score_matrix: dict,
                   indel_penalty: int) -> list:
    grid = [(0, None)]
    one_len = len(one)
    two_len = len(two)
    cur = 0
    # change to list of lists of tuples, instead of list of tuples
    # tuples can be (val, row, col) instead of (val, index)
    for one_i in range(-1, one_len):
        for two_i in range(-1, two_len):
            grid.append((float('-inf'), None))
            on_left_edge = cur % (one_len + 1) == 0
            on_top_row = cur < one_len + 1
            
            if not on_left_edge:
                val_if_used = grid[cur - 1][0] + indel_penalty
                if val_if_used > grid[cur][0]:
                    grid[cur] = (val_if_used, cur - 1)
            if not on_top_row:
                val_if_used = grid[cur - (one_len + 1)][0] + indel_penalty
                if val_if_used > grid[cur][0]:
                    grid[cur] = (val_if_used, cur - (one_len + 1))
            if (not on_left_edge) and (not on_top_row):
                val_if_used = (grid[cur - 1 - (one_len + 1)][0] +
                               score_matrix[one[one_i]][two[two_i]])
                if val_if_used > grid[cur][0]:
                    grid[cur] = (val_if_used, cur - 1 - (one_len + 1))
            cur += 1
    print(grid)
    return grid

def backtrack_alignment(one: str, two: str, grid: list) -> list:
    alignment = ['', '']
    one_len = len(one)
    two_len = len(two)
    one_i = one_len - 1
    two_i = two_len - 1
    cur = (one_len + 1) * (two_len + 1) - 1
    while cur != 0:
        diff = cur - grid[cur][1]
        if diff == 1:
            alignment[0] += one[one_i]
            alignment[1] += '-'
            one_i -= 1
        elif diff == one_len + 1:
            alignment[0] += '-'
            alignment[1] += two[two_i]
            two_i -= 1
        elif diff == one_len + 2:
            alignment[0] += one[one_i]
            alignment[1] += two[two_i]
            one_i -= 1
            two_i -= 1
        cur = grid[cur][1]
    return alignment

def align(one: str, two: str, score_matrix: dict,
          indel_penalty: int) -> (int, list):
    grid = calc_grid(one, two, score_matrix, indel_penalty)
    score = grid[(len(one) + 1) * (len(two) + 1) - 1][0]
    return score, backtrack_alignment(one, two, grid)

def read_score_matrix(file_name: str) -> dict:
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
    aligner = ProteinAligner(one, two, score_matrix, -5)
    score, alignment = aligner.align()
    print(score, *alignment, sep='\n')
