from Bio import pairwise2
import numpy as np


class MSA:
    def __init__(self, sequences):
        self.sequences = sequences
        self.n = len(self.sequences)
        self.distance_matrix = np.zeros(shape=(len(self.sequences), len(self.sequences)))
        self.divergence = np.zeros(len(sequences))

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
        for i in range(self.sequences):
            self.divergence[i] = np.sum(self.sequences[i])

    def new_distance_matrix(self):
        new = np.zeros(shape=(self.n, self.n))
        for i in range(self.n):
            for j in range(i + 1, self.n):
                temp = self.distance_matrix[i][j] - (self.divergence[i] + self.divergence[j]) / (self.n - 2)
                new[i][j] = temp
                new[j][i] = temp

    def choose_neighbor(self):
        min = self.distance_matrix[0][0]
        row, column = 0, 0
        for i in range(self.n):
            for j in range(i+1, self.n):
                if self.distance_matrix[i][j] < min:
                    min = self.distance_matrix[i][j]
                    row, column = i, j

        return row, column

    def u_distances(self, row, column):
        row_u =

if __name__ == "__main__":
    n = int(input("enter n"))

    sequences = np.zeros([])
    for i in range(n):
        sequence = input("enter sequence")
        sequences = np.append(sequences, sequence)
