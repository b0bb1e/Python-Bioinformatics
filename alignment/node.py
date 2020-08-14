class Node:
    """A node in a directed, weighted graph

    add_path: add a path leading into the node
    minimize_value: set the value of this node to a minimum
    maximize_value: choose an optimal path entering the node,
                    set value and backtrack accordingly

    read-only attributes: id, value, backtrack
    """
    
    def __init__(self, id: int, paths: dict=None):
        """Initialize all attributes to defaults/passed in values"""
        self._id = id
        if paths is None:
            paths = {}
        self._paths = paths
        self._value = 0
        self._backtrack = None

    def add_path(self, node, value: int):
        """Add a weighted path leading in to this node

        :param node: the node at the start of this path
        :type node: Node
        :param value: the weight of this path
        :type value: int
        """
        
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

    def minimize_value(self):
        """Set personal value to a minimum"""
        self._value = float('-inf')

    def maximize_value(self):
        """Calculate optimal path/backtrack to get highest value"""
        self.minimize_value()
        for back_node in self._paths:
            # calculate weight if this path is used
            if_used = self._paths[back_node] + back_node.value
            if if_used > self._value:
                self._value = if_used
                self._backtrack = back_node

    def __str__(self) -> str:
        """Get a string representation of the node

        :returns: information about this node
        :rtype: string
        """
        
        ret = 'id#' + str(self._id) + ', paths:'
        for end in self._paths:
            ret += ' id#' + str(end.id) + ':' + str(self._paths[end])
        return ret
