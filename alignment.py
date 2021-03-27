import numpy as np
from scoring_matrix import *


class Alignment:
    def __init__(self, gap_penalty: int, a: str, b: str):
        self.scoring_matrix = PAM250
        self.gap_penalty = gap_penalty
        self.a = a
        self.b = b
        rows, cols = (len(a) + 1, len(b) + 1)
        self.matrix = np.zeros((rows, cols))

    def local_alignment(self):
        self.initialization()
        self.fill_matrix()
        print(self.matrix)

    def initialization(self):
        self.matrix[0][0] = 0
        for i in range(1, len(self.matrix)):
            self.matrix[i][0] = -5 * i
        for i in range(1, len(self.matrix[0])):
            self.matrix[0][i] = -5 * i

    def fill_matrix(self):
        for i in range(1, len(self.matrix)):
            for j in range(1, len(self.matrix[0])):
                self.matrix[i][j] = self.score(i, j)

    def score(self, i: int, j: int) -> int:
        s = self.scoring_matrix[self.a[i - 1]][self.b[j - 1]] # match/mismatch score
        m = self.matrix[i - 1][j - 1] + s # m stands for max
        gap = self.matrix[i-1][j] + self.gap_penalty
        if gap > m:
            m = gap
        gap = self.matrix[i][j-1] + self.gap_penalty
        if gap > m:
            m = gap
        return m
