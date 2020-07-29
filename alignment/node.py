class Node:
    def __init__(self, id: int, paths: dict=None):
        self._id = id
        if paths is None:
            paths = {}
        self._paths = paths
        self._value = None
        self._backtrack = None

    def add_path(self, node, value: int):
        self._paths[node] = value

    @property
    def id(self) -> int:
        return self._id
    
    @property
    def value(self) -> int:
        return self._value

    @property
    def backtrack(self):
        return self._backtrack

    def maximize_value(self):
        if not self._paths:
            self._value = 0
        else:
            self._value = float('-inf')
            for back_node in self._paths:
                if back_node.value is None:
                    raise RuntimeError('A back-node (id#' + str(back_node.id) +
                                       ') does not have its value yet')
                else:
                    if_used = self._paths[back_node] + back_node.value
                    if if_used > self._value:
                        self._value = if_used
                        backtrack = back_node

    def __str__(self):
        ret = 'id#' + str(self._id) + ', paths:'
        for end in self._paths:
            ret += ' id#' + str(end.id) + ':' + str(self._paths[end])
        return ret
