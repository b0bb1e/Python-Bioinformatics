import dag

class ProteinAligner(dag.DAG):
    def __init__(self, one: str, two: str, score_matrix: dict,
                 indel_penalty: int):
        """Initialize paths using scoring matrix

        :param one: the protein that runs along the 'top'
        :type one: str
        :param two: the protein that runs along the 'side'
        :type two: str
        :param score_matrix: a protein scoring matrix
        :type score_matrix: dict (strs: dicts (strs: ints))
        :param indel_penalty: the penalty for using an indel in alignment
        :type indel_penalty: int (negative)
        """

        one_len = len(one)
        two_len = len(two)

        def node_id(one_i: int, two_i: int) -> int:
            return (one_i + 1) + ((two_i + 1) * (one_len + 1))
        paths = {0: {}}
        for i in range(one_len):
            cur_id = node_id(i, two_len - 1)
            paths[cur_id] = {cur_id + 1: indel_penalty}
        for i in range(two_len):
            cur_id = node_id(one_len - 1, i)
            paths[cur_id] = {cur_id + (one_len + 1): indel_penalty}
        for one_i in range(-1, one_len - 1):
            for two_i in range(-1, two_len -1):
                cur_id = node_id(one_i, two_i)
                paths[cur_id] = {}
                paths[cur_id][cur_id +1] = indel_penalty
                paths[cur_id][cur_id + (one_len + 1)] = indel_penalty
                paths[cur_id][cur_id + 1 + (one_len + 1)] = (
                    score_matrix[one[one_i + 1]][two[two_i + 1]])
        self._one = one
        self._two = two
        super().__init__(0, ((one_len + 1) * (two_len + 1)) - 1, paths)
        print('finished set-up')
    
    def align(self) -> (int, list):
        score, path = super().longest_path()
        print('finished path-finding')
        alignment = ['', '']
        last = path[0]
        one_i = 0
        two_i = 0
        one_len = len(self._one)
        for node in path:
            diff = node - last
            if diff == 1:
                alignment[0] += self._one[one_i]
                alignment[1] += '-'
                one_i += 1
            elif diff == one_len + 1:
                alignment[0] += '-'
                alignment[1] += self._two[two_i]
                two_i += 1
            elif diff == one_len + 2:
                alignment[0] += self._one[one_i]
                alignment[1] += self._two[two_i]
                one_i += 1
                two_i += 1
            last = node
        return score, alignment

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
