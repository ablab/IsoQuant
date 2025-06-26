import numpy as np
from collections import defaultdict, deque

class CellTypeTree:

    def __init__(self, cell_types: set):

        self.children = defaultdict(set)
        self.parents = defaultdict(str)
        self.cell_types = None
        self.n_of_cell_types = 0

        self.primary_cell_types = {}

        self._construct_tree(cell_types)
        self._construct_cell_types()

    def _construct_tree(self, cell_types: set) -> None:
        for cell_type in cell_types:
            if cell_type != 'NA':
                self._add_cell_type(cell_type)
        #print(self.children)
        #print(self.parents)

    def _add_cell_type(self, cell_type: str) -> None:
        current = "ROOT"
        for sub_cell_type in cell_type.split(":"):
            self.children[current].add(sub_cell_type)
            self.parents[sub_cell_type] = current
            current = sub_cell_type
    
    def _construct_cell_types(self) -> None:
        self.cell_types = {}
        index = 0
        queue = deque(["ROOT"])
        while len(queue) > 0:
            current = queue.popleft()
            for cell_type in self.children[current]:
                queue.append(cell_type)
            self.cell_types[current] = index
            index += 1
        self.n_of_cell_types = index
        #print(self.cell_types)

    def _add_cell_count(self, cell_type: str, count: int, cell_type_counts: np.array) -> None:
        while cell_type != "":
            cell_type_counts[self.cell_types[cell_type]] += count
            cell_type = self.parents[cell_type]

    def get_leaf_types(self):
        return [cell_type for cell_type in self.cell_types if len(self.children[cell_type]) == 0]

    def transform_counts(self, counts: defaultdict) -> list:
        full_cell_type_counts = np.zeros(shape = self.n_of_cell_types, dtype = int)
        for cell_type in self.cell_types:
            if cell_type in counts: self._add_cell_count(cell_type, counts[cell_type], full_cell_type_counts)
        return full_cell_type_counts
    
    def get_cell_type_index(self, cell_type: str) -> int:
        return self.cell_types[cell_type]
    
    def get_leaft_counts(self) -> dict:
        return {cell_type: count for cell_type, count in self.cell_types.items() if self.children[cell_type] == {}}
    
    def is_leaf(self, cell_type: str):
        return cell_type not in self.parents 
    
    def get_child_leafs(self, cell_type: str):
        if self.children[cell_type] == []:
            return [cell_type]
        leafs = []
        for child in self.children[cell_type]:
            leafs += self.get_child_leafs(cell_type = child)
        return leafs
    
    def get_cell_types(self):
        return list(self.cell_types.keys())
    
    def transform_counts_to_arrays(self, counts: dict) -> dict:
        count_array = {}
        for (u, v, cell_type), k in counts.items():
            if (u, v) not in count_array:
                count_array[(u, v)] = np.zeros(self.n_of_cell_types)
            group = cell_type
            while group != "ROOT":
                count_array[(u, v)][self.get_cell_type_index(group)] += k
                group = self.parents[group]
            count_array[(u,v)][self.get_cell_type_index(group)] += k
        return count_array
    

    def transform_count_to_arrays(self, counts: dict) -> np.array:
        results = np.zeros(shape=(max([a for a, _ in counts]) + 1, (self.n_of_cell_types)))
        for (i, cell_type), k in counts.items():
            group = cell_type
            while group != "ROOT":
                results[i][self.get_cell_type_index(group)] += k
                group = self.parents[group]
            results[i][self.get_cell_type_index(group)] += k
        return results
        


if __name__ == "__main__":
    
    cell_types = [
        "G1:G1.1:G1.1.1",
        "G1:G1.1:G1.1.2",
        "G2:G2.1:G2.1.1",
        "G1:G1.2:G1.2.1",
        "G1:G1.2:G1.2.2",
        "G2:G2.2:G2.2.2",
        "G2:G2.1:G2.1.2",
        "G1:G1.1:G1.1.1"
        ]
    
    from random import randint, random

    cc1 = defaultdict(int)
    cc2 = defaultdict(int)
    cc3 = defaultdict(int)

    for ct in cell_types:
        ct = ct.split(":")[-1]
        if random() > 0.4: cc1[ct] = randint(1, 50)
        if random() > 0.4: cc2[ct] = randint(1, 50)
        if random() > 0.4: cc3[ct] = randint(1, 50)

    c_tree = CellTypeTree(cell_types)
    results = c_tree.transform_counts(cc1)
    print(c_tree.transform_counts(cc2))
    print(c_tree.transform_counts(cc3))
    print(c_tree.parents["ROOT"])

    for i, ct in enumerate(c_tree.cell_types.keys()):
        print(ct, results[i])

    for ct, i in c_tree.get_leaft_counts().items():
        print(ct, i)