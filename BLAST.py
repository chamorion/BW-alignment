#!/usr/bin/python3
import time
from decimal import * 
from BWT import *
from MinPQ import *
from MaxPQ import *
from Alignment import *
from LA import *
from AlgnTree import *


def read_sa(path):
    sa = []
    with open(path, 'r') as f:
        for line in f:
            sa += line.strip().split(',')
    return [int(i) for i in sa]


def calculate_prp(seq):
    getcontext().prec = 16
    prp = [seq.count('A') + seq.count('T'), seq.count('C') + seq.count('G'),\
     seq.count('A') + seq.count('T'), seq.count('C') + seq.count('G'), seq.count('N')]
    prp = [Decimal(i)/Decimal(len(seq)) for i in prp]
    return prp


def lambda_and_K(prp1, prp2):
    getcontext().prec = 16
    p = Decimal(0)
    q = Decimal(0)
    for i in range(len(prp1)):
        p += prp1[i]*prp2[i]
        for j in range(len(prp2)):
            q += prp1[i]*prp2[j] if i != j else Decimal(0)
    return Decimal.ln(q/p), (q - p)**2/q


def cal_least_length(K, M, N, LAMBDA, least_score):
    getcontext().prec = 16
    least_length = (Decimal.ln(K*M*N) - Decimal.ln(least_score))/LAMBDA
    return int(least_length.quantize(Decimal('1.'), ROUND_CEILING))


def exch_ch(sq, ch1, ch2):
    sq = sq.replace(ch1, "_")
    sq = sq.replace(ch2, ch1)
    sq = sq.replace("_", ch2)
    return sq


def complementary_strand(seq):
    seq = exch_ch(seq, 'A', 'T')
    seq = exch_ch(seq, 'C', 'G')
    return ''.join(reversed(list(seq)))


def gen_alignment_list(chr, LAMBDA, K, M, N, bw, seq, is_hbb_comp, least_length):
    start = 0
    root = None
    while start + least_length <= M:
        start_refs = bw.exact_match(seq[start:start + least_length])
        for i in start_refs:
            if root is None:
                root = AlgnTree(Alignment(chr, LAMBDA, K, M, N, [start, start + least_length], [i, i + least_length], is_hbb_comp))
            else:
                root.insert(AlgnTree(Alignment(chr, LAMBDA, K, M, N, [start, start + least_length], [i, i + least_length], is_hbb_comp)))
        start += 1
    return root.pop_all() if root is not None else []
    

def algn_top_n(minpq, maxpq, n, hbb, hbb_comp, genome):
    for i in range(n):
        if maxpq.size > 0:
            algn = maxpq.pop_max()
            seq = hbb_comp if algn.is_hbb_comp else hbb
            raw_score = algn.section_query[1] - algn.section_query[0]
            shortest_length = algn.section_query[1] - algn.section_query[0]
            length = min([2 * algn.section_query[0], algn.section_ref[0]])
            left_query, left_ref = ''.join(reversed(list(seq[0:algn.section_query[0]]))),\
                ''.join(reversed(list(genome[algn.section_ref[0] - length:algn.section_ref[0]])))
            line1, bars, line2, highest_value, left_length = do_alignment(left_query, left_ref, raw_score, shortest_length)
            raw_score += highest_value
            shortest_length += len(line1)
            algn.seq_query, algn.seq_pipe, algn.seq_ref = ''.join(reversed(list(line1))),\
                ''.join(reversed(list(bars))), ''.join(reversed(list(line2)))
            algn.seq_query += seq[algn.section_query[0]:algn.section_query[1]]
            algn.seq_pipe += '|' * (algn.section_query[1] - algn.section_query[0])
            algn.seq_ref += genome[algn.section_ref[0]:algn.section_ref[1]]
            length = min([2 * (len(seq) - algn.section_query[1]), len(genome) - algn.section_ref[1]])
            right_query, right_ref = seq[algn.section_query[1]:], genome[algn.section_ref[1]:algn.section_ref[1] + length]
            line1, bars, line2, highest_value, right_length = do_alignment(right_query, right_ref, raw_score, shortest_length)
            raw_score += highest_value
            algn.seq_query += line1
            algn.seq_pipe += bars
            algn.seq_ref += line2
            algn.section_ref[0] -= left_length[0] - 1
            algn.section_query[0] -= left_length[1] - 1
            algn.section_ref[1] += right_length[0]
            algn.section_query[1] += right_length[1]
            algn.cal_e_value(raw_score)
            algn.cal_similarity()
            minpq.insert(algn)


def do_blast(chr, genome, hbb, hbb_comp, minpq):
    getcontext().prec = 16
    print(chr)
    bw = BWT(chr)
    print('alphabet of {0}: {1}'.format(chr, bw.alphabet))
    prp_hbb = calculate_prp(hbb)
    prp_genome = calculate_prp(genome)
    M = len(hbb)
    N = len(genome)
    LAMBDA, K = lambda_and_K(prp_hbb, prp_genome)
    print('length of HBB: {0}, length of {1}: {2}, lambda: {3}, K: {4}'.format(M, chr, N, LAMBDA, K))
    # maxpq = MaxPQ()
    least_length = cal_least_length(K, Decimal(M), Decimal(N), LAMBDA, Decimal(3))
    print('least length: {}'.format(least_length))
    algn_list_f = gen_alignment_list(chr, LAMBDA, K, M, N, bw, hbb, False, least_length)
    algn_list_b = gen_alignment_list(chr, LAMBDA, K, M, N, bw, hbb_comp, True, least_length)
    algn_top_n(minpq, MaxPQ(algn_list_f + algn_list_b), 10, hbb, hbb_comp, genome)
    print('size of minpq: {}'.format(minpq.size))


t_0 = time.time()
hbb = ''
with open('/home/zhangtk2022/homework/L7/HBB.fa', 'r') as f:
    for line in f:
        hbb += line.strip()
hbb_comp = complementary_strand(hbb)
genome = ''
chr_list = ['chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',\
    'chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY']
minpq = MinPQ()
count = 0
chr_dic = {}
with open("/data/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa","r") as f:
    for line in f:
        if count < len(chr_list) and chr_list[count] in line:
            if len(genome) != 0:
                do_blast(chr, genome, hbb, hbb_comp, minpq)
            chr = chr_list[count]
            chr_dic[chr] = {True:set(), False:set()}
            count += 1
            genome = ''
        elif '>K' in line:
            break
        else:
            genome += line.strip()
if len(genome) != 0:
    do_blast(chr, genome, hbb, hbb_comp, minpq)
count = 0
while minpq.size > 0 and count < 10:
    algn = minpq.pop_min()
    interval = set(range(algn.section_ref[0] - 1, algn.section_ref[1]))
    if interval.isdisjoint(chr_dic[algn.chr][algn.is_hbb_comp]):
        chr_dic[algn.chr][algn.is_hbb_comp] = chr_dic[algn.chr][algn.is_hbb_comp] | interval
        print(algn)
        count += 1
t_end = time.time()
print('time: {}'.format(t_end - t_0))
