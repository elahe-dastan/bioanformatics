import copy


class Tree:
    def __init__(self, nodes):
        self.nodes = nodes


class Node:
    def __init__(self, seq, tag):
        self.seq = seq
        self.tag = tag
        self.s = {'A': 1000000000, 'C': 1000000000, 'G': 1000000000, 'T': 1000000000}
        self.daughter = None
        self.son = None


def small_parsimony(t: Tree, character):
    for _, v in t.nodes.items():
        if v.tag == 1:
            for k in ['A', 'C', 'G', 'T']:
                if v.seq[character] == k:
                    v.s[k] = 0
    while True:
        root = None
        for _, v in t.nodes.items():
            if v.tag == 0 and v.daughter.tag == 1 and v.son.tag == 1:
                root = v
                v.tag = 1
                for k in ['A', 'C', 'G', 'T']:
                    min_i = 10000000000000
                    for i in ['A', 'C', 'G', 'T']:
                        d = 1
                        if i == k:
                            d = 0
                        tmp = v.daughter.s[i] + d
                        if tmp < min_i:
                            min_i = tmp
                    min_j = 10000000000000
                    for j in ['A', 'C', 'G', 'T']:
                        d = 1
                        if j == k:
                            d = 0
                        tmp = v.son.s[j] + d
                        if tmp < min_j:
                            min_j = tmp
                    v.s[k] = min_i + min_j
        b = True
        for _, v in t.nodes.items():
            if v.tag == 0:
                b = False
                break
        if b:
            min = 1000000000000
            for k in ['A', 'C', 'G', 'T']:
                if root.s[k] < min:
                    min = root.s[k]

            return min


if __name__ == "__main__":
    # read the number of nodes of the tree
    n = int(input())

    # read each edge of the tree
    nodes_dict = {}
    for _ in range(n - 1):
        nodes = input().split(" ")
        if nodes_dict.get(int(nodes[1]), -1) == -1:
            nodes_dict[int(nodes[1])] = Node("", 0)

        if nodes_dict.get(int(nodes[0]), -1) == -1:
            node = Node("", 0)
            node.daughter = nodes_dict[int(nodes[1])]
            nodes_dict[int(nodes[0])] = node
        elif nodes_dict[int(nodes[0])].daughter is None:
            nodes_dict[int(nodes[0])].daughter = nodes_dict[int(nodes[1])]
        else:
            nodes_dict[int(nodes[0])].son = nodes_dict[int(nodes[1])]

    # read the number of leaves of the tree
    m = int(input())

    leaves = []
    seq = ""
    for _ in range(m):
        leave = input().split(" ")
        node = nodes_dict[int(leave[0])]
        node.seq = leave[1]
        seq = leave[1]
        node.tag = 1

    t = Tree(nodes_dict)

    score = 0

    for i in range(len(seq)):
        score += small_parsimony(copy.deepcopy(t), i)

    print(score)
