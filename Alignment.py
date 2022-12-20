#!/usr/bin/python3
from decimal import *


class Alignment:
    is_hbb_comp = False
    chr, seq_query, seq_ref, seq_pipe = '', '', '', ''
    getcontext().prec = 16
    LAMBDA, K, e_value, similarity= Decimal(0), Decimal(0), Decimal(0), Decimal(0)
    M, N, raw_score = 0, 0, 0
    section_query, section_ref = [], []

    def __init__(self, chr, LAMBDA, K, M, N, section_query, section_ref, is_hbb_comp):
        self.chr = chr
        self.LAMBDA = LAMBDA
        self.K = K
        self.M = M
        self.N = N
        self.section_query = section_query
        self.section_ref = section_ref
        self.is_hbb_comp = is_hbb_comp

    def __lt__(self, algn) -> bool:
        return self.section_query[1] - self.section_query[0] < algn.section_query[1] - algn.section_query[0]

    def __gt__(self, algn) -> bool:
        return self.e_value > algn.e_value
    
    def __str__(self) -> str:
        print('{}, length: {}, position({}): {}, score: {}, similarity: {}'.format(self.chr, len(self.seq_query),\
             'b' if self.is_hbb_comp else 'f', str(self.section_ref[0])+'-'+str(self.section_ref[1]),\
                 self.raw_score, self.similarity.quantize(Decimal('.0001'), ROUND_DOWN)))
        print('K:{}, lambda:{}, M:{}, N:{}, e_value: {}'.format(self.K.quantize(Decimal('.0001'), ROUND_DOWN),\
             self.LAMBDA.quantize(Decimal('.0001'), ROUND_DOWN), self.M, self.N, self.e_value))
        self.print_alignment() 
        return ''

    def cal_e_value(self, raw_score):
        self.raw_score = raw_score
        self.e_value = self.K * self.M * self.N * Decimal(-1 * self.LAMBDA * self.raw_score).exp()

    def print_alignment(self):
        align_vec = [self.seq_query, self.seq_pipe, self.seq_ref]
        bin_length = 100
        length = len(align_vec[0])
        start = 0
        while start < length:
            end = (start + bin_length) if (start + bin_length) < length else length
            print(align_vec[0][start:end])
            print(align_vec[1][start:end])
            print(align_vec[2][start:end])
            print()
            start = end
    
    def cal_similarity(self):
        count = 0
        for i in range(len(self.seq_query)):
            count += 1 if self.seq_query[i] == self.seq_ref[i] else 0
        self.similarity = Decimal(count) / Decimal(len(self.seq_query))

    def is_disjoint_with(self, algn):
        return self.chr != algn.chr or self.is_hbb_comp ^ algn.is_hbb_comp or\
             self.section_query[1] <= algn.section_query[0] or algn.section_query[1] <= self.section_query[0] or\
                 self.section_ref[1] <= algn.section_ref[0] or algn.section_ref[1] <= self.section_ref[0]

    def can_merge_with(self, algn):
        return not self.is_disjoint_with(algn) and self.section_query[0] - algn.section_query[0] == self.section_ref[0] - algn.section_ref[0]

    def merge_with(self, algn):
        return Alignment(self.chr, self.LAMBDA, self.K, self.M, self.N,\
             [min(self.section_query[0], algn.section_query[0]), max(self.section_query[1], algn.section_query[1])],\
                 [min(self.section_ref[0], algn.section_ref[0]), max(self.section_ref[1], algn.section_ref[1])], self.is_hbb_comp)
