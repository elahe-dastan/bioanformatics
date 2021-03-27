import numpy as np
from scoring_matrix import *


class Alignment:
    def __init__(self, gap_penalty: int):
        self.scoring_matrix = PAM250
        self.gap_penalty = gap_penalty

    def local_alignment(self, a: str, b: str):
        rows, cols = (len(a) + 1, len(b) + 1)
        matrix = np.zeros((rows, cols))
        matrix = self.initialization(matrix)
        matrix = self.fill_matrix(matrix)
        print(matrix)

    def initialization(self, matrix: np.array) -> np.array:
        matrix[0][0] = 0
        for i in range(1, len(matrix)):
            matrix[i][0] = -5 * i
        for i in range(1, len(matrix[0])):
            matrix[0][i] = -5 * i
        return matrix

    def fill_matrix(self, matrix: np.array) -> np.array:
        for i in range(1, len(matrix)):
            for j in range(1, len(matrix[0])):
                matrix[i][j] = score
        return matrix
