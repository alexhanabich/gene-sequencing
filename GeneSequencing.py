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

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		# self.MaxCharactersToAlign = align_length
		num_row = len(seq1) + 1 if len(seq1) + 1 < align_length else align_length + 1
		num_col = len(seq2) + 1 if len(seq2) + 1 < align_length else align_length + 1
		table = [[None for i in range(num_col)] for j in range(num_row)]
		prev = [[None for i in range(num_col)] for j in range(num_row)]
		# fill base values horizontally for table and prev
		cost = 0
		for i in range(num_col):
			table[0][i] = cost
			prev[0][i] = 0
			cost += 5
		# fill base values vertically for table and prev
		cost = 0
		for i in range(num_row):
			table[i][0] = cost
			prev[i][0] = 1
			cost += 5

		prev[0][0] = None


		for i in range(1, num_row):
			for j in range(1, num_col):
				# compute left top and diagonal costs
				left = table[i][j - 1] + INDEL
				top = table[i - 1][j] + INDEL
				diagonal = table[i - 1][j - 1]
				if seq1[i - 1] == seq2[j - 1]:
					diagonal += MATCH
				else:
					diagonal += SUB
				# find the most optimal of the 3
				if left <= top and left <= diagonal:
					table[i][j] = left
					prev[i][j] = 0
				elif top <= diagonal:
					table[i][j] = top
					prev[i][j] = 1
				else:
					table[i][j] = diagonal
					prev[i][j] = 2

		i = num_row - 1
		j = num_col - 1
		score = table[i][j]
		alignment1 = collections.deque([])
		alignment2 = collections.deque([])
		# seq1 and seq2 indexes are 1 shorter than the table
		while prev[i][j] != None:
			if prev[i][j] == 0:
				alignment1.appendleft('-')
				alignment2.appendleft(seq2[j - 1])
				j -= 1
			elif prev[i][j] == 1:
				alignment1.appendleft(seq1[i - 1])
				alignment2.appendleft('-')
				i -= 1
			else:
				alignment1.appendleft(seq1[i - 1])
				alignment2.appendleft(seq2[j - 1])
				i -= 1
				j -= 1

		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')				
		
		return {'align_cost':score, 'seqi_first100':''.join(alignment1), 'seqj_first100':''.join(alignment2)}


