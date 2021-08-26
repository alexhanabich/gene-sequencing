    def populate_tables(self):
        for i in range(1, self.num_row):
            for j in range(1, self.num_col):
                # compute left top and diagonal costs
                left = math.inf if j < 0 or self.table[i][j - 1] == None else self.table[i][j - 1] + INDEL
                top = math.inf if j < 0 or i < 0 or self.table[i - 1][j + 1] == None else self.table[i - 1][j + 1] + INDEL
                diagonal = self.table[i - 1][j - 1]
                if self.seq1[i - 1] == self.seq2[j - 1]:
                    diagonal += MATCH
                else:
                    diagonal += SUB
                self.fill_cells(left, top, diagonal, i, j)