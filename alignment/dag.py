from node import Node

class DAG:
    def __init__(self, source: int, sink: int, paths: dict):
        self._nodes = {}
        print(paths)
        for start in paths:
            if not start in self._nodes:
                self._nodes[start] = Node(start)
            for end in paths[start]:
                if not end in self._nodes:
                    self._nodes[end] = Node(end)
                self._nodes[end].add_path(self._nodes[start], paths[start][end])
        self._source = self._nodes[source]
        self._sink = self._nodes[sink]
        for id in self._nodes:
            print(id, str(self._nodes[id]))

    def longest_path(self) -> (int, list):
        for id in range(self._source.id, self._sink.id + 1):
            self._nodes[id].maximize_value()
        return self._sink.value, self.backtrack_path()

    def backtrack_path(self) -> list:
        cur_node = self._sink
        path = []
        while cur_node is not self._source:
            path.append(cur_node.id)
            if cur_node.backtrack is None:
                raise RuntimeError('Backtrack has not been calcualted for this'
                                   ' node (id#' + str(cur_node.id) + ')')
            else:
                cur_node = cur_node.backtrack
        path.append(self._source.id)
        path.reverse()
        return path
