from GeneSequencing import MAXINDELS
import collections
import math

INDEL = 5
MATCH = -3
SUB = 1
MAXINDELS = 3
class Foo:
    def create_tables(self):
        table = [[None for i in range(self.num_col)] for j in range(self.num_row)]
        prev = [[None for i in range(self.num_col)] for j in range(self.num_row)]
        # fill base values horizontally for table and prev
        cost = 0
        for i in range(self.num_col):
            table[0][i] = cost
            prev[0][i] = 0
            cost += 5
        # fill base values vertically for table and prev
        cost = 0
        for i in range(self.num_row):
            table[i][0] = cost
            prev[i][0] = 1
            cost += 5
        # set prev[0][0] to None so signal the end of traversal
        prev[0][0] = None
        return table, prev


    def create_tables_banded(self):
        self.num_col = MAXINDELS * 2 + 1
        table = [[None for i in range(self.num_col)] for j in range(self.num_row)]
        prev = [[None for i in range(self.num_col)] for j in range(self.num_row)]
        cost = 0
        for i in range(MAXINDELS, (MAXINDELS * 2 + 1)):
            table[0][i] = cost
            prev[0][i] = 0
            cost += 5
        # fill base values vertically for table and prev
        cost = 0
        for i in range(MAXINDELS + 1):
            table[i][MAXINDELS - i] = cost
            prev[i][MAXINDELS - i] = 1
            cost += 5
        prev[0][0] = None
        return table, prev


    def fill_cells(self, left, top, diagonal, i, j):
        if left <= top and left <= diagonal:
            self.table[i][j] = left
            self.prev[i][j] = 0
        elif top <= diagonal:
            self.table[i][j] = top
            self.prev[i][j] = 1
        else:
            self.table[i][j] = diagonal
            self.prev[i][j] = 2


    def populate_tables(self):
        for i in range(1, self.num_row):
            for j in range(1, self.num_col):
                # compute left top and diagonal costs
                left = self.table[i][j - 1] + INDEL
                top = self.table[i - 1][j] + INDEL
                diagonal = self.table[i - 1][j - 1]
                if self.seq1[i - 1] == self.seq2[j - 1]:
                    diagonal += MATCH
                else:
                    diagonal += SUB
                self.fill_cells(left, top, diagonal, i, j)


     # visual representation of the table, o represents the diagonal val if it was a n by m table
        #       o***
        #      *o***
        #     **o***
        #    ***o***
        #    ***o***
        #    ***o***
        #    ***o***
        #    ***o***
        #    ***o***
        #    ***o**
        #    ***o*
        #    ***o
    def populate_tables_banded(self):
        for i in range(1, self.num_row):
            first_idx = (MAXINDELS + 1) - i if (MAXINDELS + 1) - i > 0 else 0
            rev_i = (self.num_row - 1) - i
            last_idx = MAXINDELS + rev_i + len(self.seq2) - len(self.seq1) if MAXINDELS + rev_i < 6 else 6
            print('rev', rev_i)
            print(last_idx)
            for j in range(first_idx, last_idx + 1):
                left = math.inf if j < 0 or self.table[i][j - 1] == None else self.table[i][j - 1] + INDEL
                top = math.inf if j > (self.num_col - 2) or i < 0 or self.table[i - 1][j + 1] == None else self.table[i - 1][j + 1] + INDEL
                diagonal = math.inf if j <= (MAXINDELS - i) else self.table[i - 1][j]
                print('diag', diagonal)
                if self.seq1[i - 1] == self.seq2[j + i - MAXINDELS - 1]:
                    diagonal += MATCH
                else:
                    diagonal += SUB
                print('diag after', diagonal)
                self.fill_cells(left, top, diagonal, i, j)



    def extract_solution(self):
        i = self.num_row - 1
        j = self.num_col - 1
        score = self.table[i][j]
        alignment1 = collections.deque([])
        alignment2 = collections.deque([])
        # seq1 and seq2 indexes are 1 shorter than the table
        while self.prev[i][j] != None:
            if self.prev[i][j] == 0:
                alignment1.appendleft('-')
                alignment2.appendleft(self.seq2[j - 1])
                j -= 1
            elif self.prev[i][j] == 1:
                alignment1.appendleft(self.seq1[i - 1])
                alignment2.appendleft('-')
                i -= 1
            else:
                alignment1.appendleft(self.seq1[i - 1])
                alignment2.appendleft(self.seq2[j - 1])
                i -= 1
                j -= 1
        return score, ''.join(alignment1), ''.join(alignment2)


    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.seq1 = seq1
        self.seq2 = seq2
        self.num_row = len(seq1) + 1 if len(seq1) + 1 < align_length else align_length + 1
        self.num_col = len(seq2) + 1 if len(seq2) + 1 < align_length else align_length + 1
        if banded:
            self.table, self.prev = self.create_tables_banded()
            self.populate_tables_banded()
        else:
            self.table, self.prev = self.create_tables()
            self.populate_tables()
        # score, alignment1, alignment2 = self.extract_solution()	
        score = self.table[self.num_row - 1][MAXINDELS]
        print(score)
        print(self.table)

seq1 = ['p','o','l','y','n','o','m','i','a','l']
seq2 = ['e','x','p','o','n','e','n','t','i','a','l']
foo = Foo()
foo.align(seq1, seq2, True, 1000)