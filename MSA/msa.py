"""
this module implements multiple sequence alignment (msa).
"""
import functools
import typing


class MSA:
    """
    progressive alignment builds up a final MSA by combining pairwise
    alignments beginning with the most similar pair
    and progressing to the most distantly related.
    all progressive alignment methods require two stages:

    - a first stage in which the relationships between the sequences
    are represented as a tree, called a guide tree
    - a second step in which the MSA is built by adding the sequences
    sequentially to the growing MSA according to the guide tree.
    """

    def __init__(self, sequences):
        """
        initiates all required data structure
        """
        self.base_sequences: typing.List[str] = sequences
        self.sequences: typing.Any = sequences
        self.n = len(self.sequences)
        self.distance_matrix: typing.List[typing.List[float]] = [
            [0.0 for _ in range(self.n)] for _ in range(self.n)
        ]
        self.divergence: typing.List[float] = [0.0 for i in range(self.n)]
        self.new_distance_matrix: typing.List[typing.List[float]] = [
            [0.0 for _ in range(self.n)] for _ in range(self.n)
        ]

    @staticmethod
    def global_align(x, y, s_match, s_mismatch, s_gap):
        """
        global alignment of two sequence based on a code which is availabe
        on Quera for this homework without any modification.
        """
        A = []
        for i in range(len(y) + 1):
            A.append([0] * (len(x) + 1))
        for i in range(len(y) + 1):
            A[i][0] = s_gap * i
        for i in range(len(x) + 1):
            A[0][i] = s_gap * i
        for i in range(1, len(y) + 1):
            for j in range(1, len(x) + 1):

                A[i][j] = max(
                    A[i][j - 1] + s_gap,
                    A[i - 1][j] + s_gap,
                    A[i - 1][j - 1]
                    + (
                        s_match
                        if (y[i - 1] == x[j - 1] and y[i - 1] != "-")
                        else 0
                    )
                    + (
                        s_mismatch
                        if (
                            y[i - 1] != x[j - 1]
                            and y[i - 1] != "-"
                            and x[j - 1] != "-"
                        )
                        else 0
                    )
                    + (s_gap if (y[i - 1] == "-" or x[j - 1] == "-") else 0),
                )

        align_X = ""
        align_Y = ""
        i = len(x)
        j = len(y)

        while i > 0 or j > 0:

            current_score = A[j][i]

            if (
                i > 0
                and j > 0
                and (
                    (
                        (x[i - 1] == y[j - 1] and y[j - 1] != "-")
                        and current_score == A[j - 1][i - 1] + s_match
                    )
                    or (
                        (
                            y[j - 1] != x[i - 1]
                            and y[j - 1] != "-"
                            and x[i - 1] != "-"
                        )
                        and current_score == A[j - 1][i - 1] + s_mismatch
                    )
                    or (
                        (y[j - 1] == "-" or x[i - 1] == "-")
                        and current_score == A[j - 1][i - 1] + s_gap
                    )
                )
            ):
                align_X = x[i - 1] + align_X
                align_Y = y[j - 1] + align_Y
                i = i - 1
                j = j - 1
            elif i > 0 and (current_score == A[j][i - 1] + s_gap):
                align_X = x[i - 1] + align_X
                align_Y = "-" + align_Y
                i = i - 1
            else:
                align_X = "-" + align_X
                align_Y = y[j - 1] + align_Y
                j = j - 1
        return (align_X, align_Y, A[len(y)][len(x)])

    def fill_distance_matrix(self):
        """
        distance matrix contains the score of pairwaise global
        alignment between sequences.
        only in the first time we calculate distance matrix by
        global alignments, each next steps will calculate it by
        neighbor joining algorithm.
        """
        for i in range(self.n):
            for j in range(i + 1, self.n):
                _, _, score = self.global_align(
                    self.sequences[i], self.sequences[j], 1, -1, -2
                )
                self.distance_matrix[i][j] = score
                self.distance_matrix[j][i] = score

    def guide_tree(self):
        """
        generate a guide by repeating the steps that are discussed
        here:

        https://www.deduveinstitute.be/~opperd/private/neighbor.html
        https://en.wikipedia.org/wiki/Neighbor_joining
        """
        while self.n > 2:
            self.calculate_divergence()
            self.build_new_distance_matrix()
            row, column = self.choose_neighbor()
            self.u_distances(row, column)
            self.distances_from_u(row, column)

        root = (self.sequences[0], self.sequences[1])
        _, _, common = self.traverse(root)

        score = 0
        for seq in self.base_sequences:
            _, res, sc = self.global_align(common, seq, 1, -1, -2)
            score += sc
            print(res)
        print(score)

    @staticmethod
    def consensus(*ss: str) -> str:
        n = len(ss[0])

        result = ""
        for i in range(n):
            scores: typing.Dict[str, int] = {}
            for s in ss:
                if s[i] == "-":
                    # continue on gap, we will handle it later
                    pass
                elif s[i] in scores:
                    scores[s[i]] += 1
                else:
                    scores[s[i]] = 1

            max_score = 0
            max_ch = "-"
            for ch in scores:
                if scores[ch] > max_score:
                    max_score = scores[ch]
                    max_ch = ch
                elif scores[ch] == max_score and ch < max_ch:
                    max_score = scores[ch]
                    max_ch = ch
            result += max_ch

        return result

    @staticmethod
    def traverse(root) -> typing.Tuple[str, str, str]:
        """
        calculate node alignments.
        leaves correspond to sequences.
        internal nodes represent alignments.
        """
        r1 = root[0]
        r2 = root[1]

        rr1: str = ""
        br11: str = ""
        br12: str = ""

        rr2: str = ""
        br21: str = ""
        br22: str = ""

        if isinstance(r1, tuple):
            # intermediate node
            br11, br12, rr1 = MSA.traverse(r1)
        else:
            # leaf node
            rr1 = r1

        if isinstance(r2, tuple):
            # intermediate node
            br21, br22, rr2 = MSA.traverse(r2)
        else:
            # leaf node
            rr2 = r2

        ai, aj, _ = MSA.global_align(rr1, rr2, 1, -1, -2)

        if isinstance(r1, tuple) and isinstance(r2, tuple):
            return ai, aj, MSA.consensus(ai, aj)

        if isinstance(r1, tuple):
            return ai, aj, MSA.consensus(ai, aj, br11, br12)

        if isinstance(r2, tuple):
            return ai, aj, MSA.consensus(ai, aj, br21, br22)

        return ai, aj, MSA.consensus(ai, aj)

    def calculate_divergence(self):
        """
        calculate the divergence based on the distance matrix.
        """
        for i in range(self.n):
            self.divergence[i] = (
                functools.reduce(lambda x, y: x + y, self.distance_matrix[i])
                - self.distance_matrix[i][i]
            )

    def build_new_distance_matrix(self):
        """
        build a new distance matrix whic is called
        Q-matrix in wikipedia article.
        """
        for i in range(self.n):
            for j in range(i + 1, self.n):
                distance = self.distance_matrix[i][j] - (
                    self.divergence[i] + self.divergence[j]
                ) / (self.n - 2)
                self.new_distance_matrix[i][j] = distance
                self.new_distance_matrix[j][i] = distance

    def choose_neighbor(self):
        """
        choose two neighbors that has the minimum score.
        here we only check the strict lower than and not the equal
        and lower than to be compatible with the Quera description.
        """
        row, column = 0, 0
        min_score = self.new_distance_matrix[row][column]
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if self.new_distance_matrix[i][j] < min_score:
                    min_score = self.new_distance_matrix[i][j]
                    row, column = i, j

        return row, column

    def u_distances(self, seq_a, seq_b):
        """
        referes to

        https://en.wikipedia.org/wiki/Neighbor_joining#Distance_from_the_pair_members_to_the_new_node
        """
        u_seq_a_seq_b = self.distance_matrix[seq_a][seq_b] / 2 + (
            self.divergence[seq_a] - self.divergence[seq_b]
        ) / (2 * (self.n - 2))

        u_seq_b_seq_a = self.distance_matrix[seq_a][seq_b] - u_seq_a_seq_b

        return u_seq_a_seq_b, u_seq_b_seq_a

    def distances_from_u(self, seq_a, seq_b):
        """
        generate a new distance matrix by merging
        the seq_a and seq_b sequences.
        """
        # others is a set of remaining sequences
        others = set(range(self.n)) - {seq_a, seq_b}

        new_matrix: typing.List[typing.List[float]] = [
            [0 for _ in range(self.n - 1)] for _ in range(self.n - 1)
        ]

        # create a map between old and new indecies
        new_indecies = {}
        for new_index, old_index in enumerate(others):
            new_indecies[new_index + 1] = old_index

        for i in range(1, self.n - 1):
            new_matrix[i][0] = (
                self.distance_matrix[seq_a][new_indecies[i]]
                + self.distance_matrix[seq_b][new_indecies[i]]
                - self.distance_matrix[seq_a][seq_b] / 2
            )
            new_matrix[0][i] = new_matrix[i][0]

        for i in range(1, self.n - 1):
            for j in range(1, self.n - 1):
                new_matrix[i][j] = self.distance_matrix[new_indecies[i]][
                    new_indecies[j]
                ]
                new_matrix[j][i] = new_matrix[i][j]

        # we are going to store guide tree into sequence list
        self.sequences = [(self.sequences[seq_a], self.sequences[seq_b])] + [
            self.sequences[new_indecies[i]] for i in range(1, self.n - 1)
        ]

        self.distance_matrix = new_matrix
        self.n = self.n - 1

    def do(self):
        """
        do the msa algorithm step by step
        """
        self.fill_distance_matrix()
        self.guide_tree()


if __name__ == "__main__":
    # read the number of sequences
    n = int(input())

    # read each sequence one by one
    sequences = []
    for _ in range(n):
        sequence = input()
        sequences.append(sequence)

    msa = MSA(sequences)
    msa.do()
