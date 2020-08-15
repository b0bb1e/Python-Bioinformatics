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
        
        paths = {0: {}}
        one_len = len(one)
        two_len = len(two)
        for i in range(one_len):
            paths[i + 1] = {i: indel_penalty}
        for i range(two_len):
            paths[(i + 1) * (one_len + 1)] = {i * (one_len + 1): indel_penalty}
        for one_i in range(one_len):
            for two_i in range(two_len):
                cur_id = (one_i + 1) + ((two_i + 1) * (one_len + 1))
                paths[cur_id] = {}
                paths[cur_id][cur_id - 1] = indel_penalty
                paths[cur_id][cur_id - (one_len - 1)] = indel_penalty
                pahts[cur_id][cur_id - 1 - (one_len - 1)] = (
                    score_matrix[one[one_i]][two[two_i]])
        super().__init__(0, ((one_len + 1) * (two_len + 1)) - 1, paths)
