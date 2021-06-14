class Edge:
    def __init__(self, i, j):
        self.i = i
        self.j = j


class Leave:
    def __init__(self, num, seq):
        self.num = num
        self.seq = seq


if __name__ == "__main__":
    # read the number of nodes of the tree
    n = int(input())

    # read each edge of the tree
    edges = []
    for _ in range(n - 1):
        nodes = input().split(" ")
        edges.append(Edge(int(nodes[0]), int(nodes[1])))

    # read the number of leaves of the tree
    m = int(input())

    leaves = []
    for _ in range(m):
        leave = input().split(" ")
        leaves.append(Leave(int(leave[0]), leave[1]))
