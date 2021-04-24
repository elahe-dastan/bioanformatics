from Bio import pairwise2
import numpy as np


class MSA:
    def __init__(self, sequences):
        self.sequences = sequences
        self.n = len(self.sequences)
        self.distance_matrix = np.zeros(shape=(len(self.sequences), len(self.sequences)))
        self.divergence = np.zeros(len(sequences))
        self.new = np.zeros(shape=(self.n, self.n))

    def fill_distance_matrix(self):
        for i in range(self.n):
            for j in range(i + 1, self.n):
                temp = pairwise2.align.globalms(self.sequences[i], self.sequences[j], 3, -1, -1.5, -1.5)
                # print(temp[0].score)  # prints the score of the first alignment
                # print(temp[0].seqA)  # prints the first sequence of the first alignment
                # print(temp[0].seqB)  # prints the second sequence of the first alignment
                self.distance_matrix[i][j] = temp[0].score
                self.distance_matrix[j][i] = temp[0].score
        print(self.distance_matrix)

    def guide_tree(self):
        self.new = np.zeros(shape=(self.n, self.n))
        self.calculate_divergence()
        self.new_distance_matrix()
        row, column = self.choose_neighbor()
        self.u_distances(row, column)

    def calculate_divergence(self):
        for i in range(self.n):
            self.divergence[i] = np.sum(self.distance_matrix[i]) - self.distance_matrix[i][i]

    def new_distance_matrix(self):
        for i in range(self.n):
            for j in range(i + 1, self.n):
                temp = self.distance_matrix[i][j] - (self.divergence[i] + self.divergence[j]) / (self.n - 2)
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
        row_u = self.distance_matrix[row][column] / 2 + (self.divergence[row] - self.divergence[column]) /\
                (2 * (self.n - 2))

        column_u = self.distance_matrix[row][column] - row_u

    def distances_from_u(self, row, column):
        new_matrix = np.delete(self.distance_matrix, [row, column], 0)
        new_matrix = np.delete(new_matrix, [row, column], 1)
        new_matrix = np.insert(new_matrix, 0, np.zeros(self.n - 2), axis=0)
        new_matrix = np.insert(new_matrix, 0, np.zeros(self.n - 2), axis=1)
        others = np.delete(np.delete(np.arange(self.n), row), column)
        for i in others:
            new_matrix[i][0] = self.distance_matrix[row][i] + self.distance_matrix[column][i] - self.distance_matrix[row][column] / 2
            new_matrix[0][i] = new_matrix[i][0]

    def do(self):
        self.fill_distance_matrix()
        self.guide_tree()


if __name__ == "__main__":
    n = int(input("enter n"))

    sequences = np.zeros([])
    for i in range(n):
        sequence = input("enter sequence")
        sequences = np.append(sequences, sequence)

    msa = MSA(sequences)
    msa.do()
