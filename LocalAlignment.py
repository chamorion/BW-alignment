#!/usr/bin/python3
from decimal import *


class Node:
	def __init__(self, x, y):
		self.x = x
		self.y = y
		self.value = -float('inf')
		self.left = None
		self.top_left = None
		self.top = None
		self.length = 0
	
	def get_x(self):
		return self.x
	
	def get_y(self):
		return self.y
	
	def get_value(self):
		return self.value

	def get_left(self):
		return self.left
	
	def get_top_left(self):
		return self.top_left
	
	def get_top(self):
		return self.top

	def set_value(self, value):
		self.value = value

	def set_left(self, left):
		self.left = left
	
	def set_top_left(self, top_left):
		self.top_left = top_left

	def set_top(self, top):
		self.top = top


def do_alignment(seq1, seq2, previous_score, previous_length):
    gap = -1
    mismatch = -1
    match = 1
    value_matrix = []

    # initialize value matrix
    for i in range(len(seq1) + 1):
        value_matrix.append([Node(i, j) for j in range(len(seq2) + 1)])

    # calculate value matrix
    # initialize row value
    for i in range(len(seq1) + 1):
        value_matrix[i][0].set_value(i * gap)
        value_matrix[i][0].length = i
        value_matrix[i][0].set_top(value_matrix[i-1][0] if i > 0 else None)

    # initialize column value
    for j in range(len(seq2) + 1):
        value_matrix[0][j].set_value(j * gap)
        value_matrix[0][j].length = j
        value_matrix[0][j].set_left(value_matrix[0][j-1] if j > 0 else None)

    # calculate value of every item
    highest_x = 0
    highest_y = 0
    highest_value = 0
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            index_dic = {0:value_matrix[i - 1][j], 1:value_matrix[i - 1][j - 1], 2:value_matrix[i][j - 1]}
            value_vec = [value_matrix[i - 1][j].get_value() + gap, value_matrix[i - 1][j - 1].get_value() + s,
                         value_matrix[i][j - 1].get_value() + gap]
            value_max = max(value_vec)
            length_vec = [index_dic[index].length for index,value in enumerate(value_vec) if value == value_max]
            length_min = min(length_vec)
            value_matrix[i][j].set_value(value_max)
            value_matrix[i][j].length = length_min + 1
            if Decimal(value_matrix[i][j].value + previous_score) / Decimal(value_matrix[i][j].length + previous_length) >= Decimal('0.8') and\
                     i ** 2 + j ** 2 > highest_x ** 2 + highest_y ** 2:
                     highest_x = i
                     highest_y = j
                     highest_value = value_matrix[i][j].value
            
            if value_matrix[i - 1][j].get_value() + gap == value_max and value_matrix[i - 1][j].length == length_min:
                value_matrix[i][j].set_top(value_matrix[i - 1][j])
            if value_matrix[i - 1][j - 1].get_value() + s == value_max and value_matrix[i - 1][j - 1].length == length_min:
                value_matrix[i][j].set_top_left(value_matrix[i - 1][j - 1])
            if value_matrix[i][j - 1].get_value() + gap == value_max and value_matrix[i][j - 1].length == length_min:
                value_matrix[i][j].set_left(value_matrix[i][j - 1])

    # trace backward
    line1 = ""
    bars = ""
    line2 = ""
    line1_count = 0
    line2_count = 0
    x = highest_x
    y = highest_y
    while x != 0 or y != 0:
        if value_matrix[x][y].get_top_left() is not None:
            line1 = seq1[x-1] + line1
            line2 = seq2[y-1] + line2
            bars = ("|" if seq1[x-1] == seq2[y-1] else " ") + bars
            x = x - 1
            y = y - 1
        elif value_matrix[x][y].get_top() is not None:
            line1 = seq1[x-1] + line1
            line2 = "-" + line2
            line2_count += 1
            bars = " " + bars
            x = x - 1
        elif value_matrix[x][y].get_left() is not None:
            line1 = "-" + line1
            line2 = seq2[y-1] + line2
            line1_count += 1
            bars = " " + bars
            y = y - 1

    return line1, bars, line2, highest_value, [len(line1) - line1_count, len(line2) - line2_count]
