from GeneSequencing import MAXINDELS
import collections
import math

INDEL = 5
MATCH = -3
SUB = 1
MAXINDELS = 3
LEFT = 0
TOP = 1
DIAGONAL = 2
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
		for i in range(MAXINDELS + 1):
			table[0][i] = cost
			prev[0][i] = LEFT
			cost += 5
		# fill base values vertically for table and prev
		cost = 0
		for i in range(MAXINDELS + 1):
			table[i][0] = cost
			prev[i][0] = TOP
			cost += 5
		# set prev[0][0] to None so signal the end of traversal
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
		# init loop
		for i in range(1, MAXINDELS + 1):
			temp_num_col = i + MAXINDELS + 1
			for j in range(1, temp_num_col):
				# compute left top and diagonal costs
				left = self.table[i][j - 1] + INDEL
				top = math.inf if self.table[i - 1][j] == None else self.table[i - 1][j] + INDEL
				diagonal = self.table[i - 1][j - 1]
				if self.seq1[i - 1] == self.seq2[j - 1]:
					diagonal += MATCH
				else:
					diagonal += SUB
				self.fill_cells(left, top, diagonal, i, j)
		# second loop	
		for i in range(MAXINDELS + 1, self.num_row):
			temp_i = i - MAXINDELS + 1
			temp_num_col = min(self.num_col, len(self.seq1) - temp_i)
			for j in range(temp_num_col):
				# compute left top and diagonal costs
				left = self.table[i][j - 1] + INDEL if  j > 0 else math.inf
				top = self.table[i - 1][j + 1] + INDEL if j + 1 < temp_num_col else math.inf
				diagonal = self.table[i - 1][j]
				print(i,j, left, top, diagonal)
				if self.seq1[i - 1] == self.seq2[j - 1]:
					diagonal += MATCH
				else:
					diagonal += SUB
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

	def extract_solution_banded(self):
		i = self.num_row - 1
		j = self.num_col - 1
		loopcnt = 0
		while True:
			j -= 1
			if self.table[i][j] != None:
				break
		score = self.table[i][j]
		alignment1 = 'testesttest'
		alignment2 = 'testtesttest'
		return score, alignment1, alignment2


	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.seq1 = seq1
		self.seq2 = seq2
		self.num_row = len(seq1) + 1 if len(seq1) + 1 < align_length else align_length + 1
		self.num_col = len(seq2) + 1 if len(seq2) + 1 < align_length else align_length + 1
		score = 0
		align1 = ''
		align2 = ''
		if banded:
			self.table, self.prev = self.create_tables_banded()
			self.populate_tables_banded()
			score, align1, align2 = self.extract_solution_banded()
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