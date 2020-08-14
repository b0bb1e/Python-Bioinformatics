from node import Node

class DAG:
    """A generalized Directed Acyclic Graph
    Assumes that increasing id#s is a valid topological ordering
    Subclasses only have to override the consructor, and ensure nodes
    are topologically ordered

    longest_path: calculate the longest path in the graph
    backtrack_path: find longest path. Call only after longest_path
                    does the actual calculation
    """
    
    def __init__(self, source: int, sink: int, paths: dict):
        """Initialize all the Nodes in this graph

        :param source: the id# of the source node
        :type source: int
        :param sink: the id# of the sink node
        :type sink: int
        :param paths: a lookup table of all paths from a certain node
                      {start node id#: {end node id#: weight of path}}
        :type paths: dict (of ints: dicts (of ints: ints))
        """
        
        self._nodes = {}
        for start in paths:
            if start not in self._nodes:
                self._nodes[start] = Node(start)
            for end in paths[start]:
                if end not in self._nodes:
                    self._nodes[end] = Node(end)
                self._nodes[end].add_path(self._nodes[start], paths[start][end])
        self._source = self._nodes[source]
        self._sink = self._nodes[sink]

    def longest_path(self) -> (int, list):
        """Calculate the longest path and its weight

        :returns: the path's weight and its nodes in order
        :rtype: tuple (int, list (of ints))
        """

        for id in range(0, self._source.id):
            # protection against skipped id#s
            if id in self._nodes:
                # any nodes before the source are unreachable
                self._nodes[id].minimize_value()
        for id in range(self._source.id + 1, self._sink.id + 1):
            if id in self._nodes:
                self._nodes[id].maximize_value()
        return self._sink.value, self.backtrack_path()

    def backtrack_path(self) -> list:
        """Backtrack to find the longest path

        :returns: the path's nodes in order
        :rtype: list (of ints)
        """
        
        cur_node = self._sink
        path = []
        while cur_node is not self._source:
            path.append(cur_node.id)
            if cur_node.backtrack is None:
                raise RuntimeError('Backtrack has not been calcualted for this'
                                   ' node (id#' + str(cur_node.id) + ')')
            else:
                # move backwards in optimal path
                cur_node = cur_node.backtrack
        # source was not added within the loop
        path.append(self._source.id)
        # backtracking added nodes in reverse order
        path.reverse()
        return path

    def __str__(self) -> str:
        """Get a string representation of the graph

        :returns: information about all nodes in the graph
        :rtype: string
        """
        
        for id in self._nodes:
            print(self._nodes[id])

if __name__ == '__main__':
    with open('data.txt') as data:
        source = int(data.readline().rstrip())
        sink = int(data.readline().rstrip())
        
        paths = {}
        for line in data:
            # extract information from the line
            start = int(line.split('->')[0])
            end = int(line.split('->')[1].split(':')[0])
            value = int(line.split(':')[1].rstrip())
            if not start in paths:
                paths[start] = {}
            paths[start][end] = value
    graph = DAG(source, sink, paths)
    path_len, path = graph.longest_path()
    print(path_len)
    print(*path, sep='->')
