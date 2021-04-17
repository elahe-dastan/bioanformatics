from Bio import pairwise2
import numpy as np

class MSA:
    def __init__(self, sequences):
        self.sequences = sequences
        self.distance_matrix = np.zeros(shape=(len(self.sequences), len(self.sequences)))

    def fill_distance_matrix(self):
        for i in range(len(self.sequences)):
            for j in range(i+1, len(self.sequences)):
                temp = pairwise2.align.globalms(self.sequences[i], self.sequences[j], 3, -1, -1.5, -1.5)
                # print(temp[0].score)  # prints the score of the first alignment
                # print(temp[0].seqA)  # prints the first sequence of the first alignment
                # print(temp[0].seqB)  # prints the second sequence of the first alignment
                self.distance_matrix[j][i] = temp[0].score
        print(self.distance_matrix)

    def guide_tree(self):



if __name__ == "__main__":
    n = int(input("enter n"))

    sequences = np.zeros([])
    for i in range(n):
        sequence = input("enter sequence")
        sequences = np.append(sequences, sequence)
