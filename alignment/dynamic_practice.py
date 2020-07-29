def min_coins(value: int, denoms: list, pre_calc: dict=None) -> int:
    """Finds the minimum number of coins needed for change

    :param value: the value to change
    :type value: int
    :param denoms: the denominations of coins available
    :type denoms: list (of ints)
    :param pre_calc: dynamic programming dictionary, with {value: return}
                     defaults to None (if no values have been calculated)
    :type pre_calc: dict (int: int)
    :returns: the minimum number of coins
    :rtype: int
    """

    if pre_calc is None:
        pre_calc = {}
    if value < 0:
        raise ValueError('Cannot change a negative amount')
    elif value == 0:
        return 0
    elif value in pre_calc:
        return pre_calc[value]

    cur_min = value
    for denom in denoms:
        value_if_used = value - denom
        if value_if_used >= 0:
            coins_if_used = 1 + min_coins(value_if_used, denoms, pre_calc)
            if cur_min > coins_if_used:
                cur_min = coins_if_used
    return cur_min

def max_length(down: list, right: list) -> int:
    """Finds the maximum length of a R/D path given edge weights

    :param down: a matrix of down-path weights
    :type down: list (of lists (of ints))
    :param right: a matrix of right-path weights
    :type right: list (of lists (of ints))
    :returns: the maximal path weight through the grid
    :rtype: int
    """

    max_weights = [[0]]
    width = len(down[0])
    height = len(right)
    
    for i in range(width - 1):
        max_weights[0].append(right[0][i] + max_weights[0][i])
    for i in range(height - 1):
        max_weights.append([down[i][0] + max_weights[i][0]])

    for row in range(1, height):
        for col in range(1, width):
            max_weights[row].append(max(
                max_weights[row - 1][col] + down[row - 1][col],
                max_weights[row][col - 1] + right[row][col - 1]))

    return max_weights[height - 1][width - 1]

def build_string(backtrack: list, one: str) -> str:
    """Builds a string based on a backtarack matrix

    :param backtrack: the backtrack matrix
    :type backtrack: list (of lists (of strs))
    :param one: the string along the top of the matrix
    :type one: str
    :returns: the string as backtracked
    :rtype: str
    """

    row, col = len(backtrack) - 1, len(backtrack[0]) - 1
    built = ""
    while row != 0 or col != 0:
        if backtrack[row][col] == 'd':
            row -= 1
        elif backtrack[row][col] == 'r':
            col -= 1
        elif backtrack[row][col] == 'b':
            built = one[col - 1] + built
            row -= 1
            col -= 1
    return built
    

def lcs(one: str, two: str) -> str:
    """Finds the longest common subsequence of two strings

    :param one: one of the strings
    :type one: str
    :param two: the other string
    :type two: str
    :returns: the longest common subsequence
    :rtype: str
    """

    width = len(one) + 1
    height = len(two) + 1
    weights = [[0] for i in range(height)]
    backtrack = [['d'] for i in range(height)]
    weights[0] += [0 for i in range(width - 1)]
    backtrack[0] += ['r' for i in range(width - 1)]

    for row in range(1, height):
        for col in range(1, width):
            weights[row].append(weights[row - 1][col])
            backtrack[row].append('d')
            if_right = weights[row][col - 1]
            if_diag = weights[row - 1][col - 1] + (one[col - 1] == two[row - 1])
            if if_right > weights[row][col]:
                weights[row][col] = if_right
                backtrack[row][col] = 'r'
            elif if_diag > weights[row][col]:
                weights[row][col] = if_diag
                backtrack[row][col] = 'b'
    return build_string(backtrack, one)
                     
if __name__ == '__main__':
    with open('data.txt') as data:
        one = data.readline().rstrip()
        two = data.readline().rstrip()
    print(lcs(one, two))
