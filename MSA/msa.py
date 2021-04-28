import typing
import functools


class MSA:
    def __init__(self, sequences):
        self.sequences = sequences
        self.n = len(self.sequences)
        self.distance_matrix: typing.List[typing.List[int]] = [
            [0 for _ in range(self.n)] for _ in range(self.n)
        ]
        self.divergence = [0 for i in range(self.n)]
        self.new_distance_matrix = [
            [0 for _ in range(self.n)] for _ in range(self.n)
        ]

    @staticmethod
    def global_align(x, y, s_match, s_mismatch, s_gap):
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
        """
        for i in range(self.n):
            for j in range(i + 1, self.n):
                _, _, score = self.global_align(
                    self.sequences[i], self.sequences[j], 1, -1, -2
                )
                self.distance_matrix[i][j] = score
                self.distance_matrix[j][i] = score
        print(self.distance_matrix)

    def guide_tree(self):
        self.calculate_divergence()
        self.build_new_distance_matrix()
        row, column = self.choose_neighbor()
        self.u_distances(row, column)

    def calculate_divergence(self):
        for i in range(self.n):
            self.divergence[i] = (
                functools.reduce(lambda x, y: x + y, self.distance_matrix[i])
                - self.distance_matrix[i][i]
            )

    def build_new_distance_matrix(self):
        for i in range(self.n):
            for j in range(i + 1, self.n):
                temp = self.distance_matrix[i][j] - (
                    self.divergence[i] + self.divergence[j]
                ) / (self.n - 2)
                self.new[i][j] = temp
                self.new[j][i] = temp

    def choose_neighbor(self):
        min = self.new[0][0]
        row, column = 0, 0
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if self.new[i][j] < min:
                    min = self.new[i][j]
                    row, column = i, j

        return row, column

    def u_distances(self, row, column):
        row_u = self.distance_matrix[row][column] / 2 + (
            self.divergence[row] - self.divergence[column]
        ) / (2 * (self.n - 2))

        column_u = self.distance_matrix[row][column] - row_u

    def distances_from_u(self, row, column):
        new_matrix = np.delete(self.distance_matrix, [row, column], 0)
        new_matrix = np.delete(new_matrix, [row, column], 1)
        new_matrix = np.insert(new_matrix, 0, np.zeros(self.n - 2), axis=0)
        new_matrix = np.insert(new_matrix, 0, np.zeros(self.n - 2), axis=1)
        others = np.delete(np.delete(np.arange(self.n), row), column)
        for i in others:
            new_matrix[i][0] = (
                self.distance_matrix[row][i]
                + self.distance_matrix[column][i]
                - self.distance_matrix[row][column] / 2
            )
            new_matrix[0][i] = new_matrix[i][0]

    def do(self):
        self.fill_distance_matrix()
        self.guide_tree()


if __name__ == "__main__":
    n = int(input("enter n"))

    sequences = []
    for _ in range(n):
        sequence = input("enter sequence")
        sequences.append(sequence)

    msa = MSA(sequences)
    msa.do()
