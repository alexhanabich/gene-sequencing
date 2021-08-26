#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random
import collections

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment
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
				# find the most optimal of the 3
				if left <= top and left <= diagonal:
					self.table[i][j] = left
					self.prev[i][j] = 0
				elif top <= diagonal:
					self.table[i][j] = top
					self.prev[i][j] = 1
				else:
					self.table[i][j] = diagonal
					self.prev[i][j] = 2


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
		return score, alignment1, alignment2


	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.seq1 = seq1
		self.seq2 = seq2
		self.num_row = len(seq1) + 1 if len(seq1) + 1 < align_length else align_length + 1
		self.num_col = len(seq2) + 1 if len(seq2) + 1 < align_length else align_length + 1
		self.table, self.prev = self.create_tables()
		self.populate_tables()
		score, alignment1, alignment2 = self.extract_solution()
		

		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')				
		
		return {'align_cost':score, 'seqi_first100':''.join(alignment1), 'seqj_first100':''.join(alignment2)}


