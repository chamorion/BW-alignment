#!/usr/bin/python3
import time


def gen_t(s):
    t = [''] * len(s)
    t[-1] = 'S'
    i = len(s) - 2
    while i >= 0:
        if s[i] < s[i + 1]:
            t[i] = 'S'
        elif s[i] > s[i + 1]:
            t[i] = 'L'
        else:
            t[i] = t[i + 1]
        i -= 1
    return t


def gen_dic(s, k):
    if s == str(s):
        str_set = [None] * 256
        for ch in s:
            str_set[ord(ch)] = ch
    else:
        str_set = [None] * k
        for i in s:
            str_set[i] = i
    str_set = [i for i in str_set if i is not None]
    index_dic = {ch : index for index, ch in enumerate(str_set)}
    return index_dic


def is_lms(t, i):
    return i > 0 and t[i - 1] == 'L' and t[i] == 'S'


def induced_sort(s, p, index_dic, bucket, k, t):
    sa = [-1] * len(s)

    # step1: insert every LMS to the end of bucket in sa
    count = [0] * k
    for lms in p:
        ch = s[lms]
        index = index_dic[ch]
        sa[bucket[index] - 1 - count[index]] = lms
        count[index] += 1

    # step2: scan sa from left to right and insert L type to the head of the bucket in sa
    count = [0] * k
    for i in sa:
        if i <= 0 or t[i - 1] == 'S':
            continue
        ch = s[i - 1]
        index = index_dic[ch]
        sa[bucket[index - 1] + count[index]] = i - 1
        count[index] += 1

    # step3: scan sa from right to left and insert S type to the end of the bucket in sa
    count = [0] * k
    for i in reversed(sa):
        if i <= 0 or t[i - 1] == 'L':
            continue
        ch = s[i - 1]
        index = index_dic[ch]
        sa[bucket[index] - 1 - count[index]] = i - 1
        count[index] += 1

    return sa


def gen_sa(s, k):
    if len(s) < 2:
        return s
    
    # generate type array and position array of left-most-S substrings
    t = gen_t(s)
    p = [i for i in range(len(s)) if is_lms(t, i)]

    # generate dictionary of {character : index} and bucket of every character
    index_dic = gen_dic(s, k)
    bucket = [0] * k
    for i in s:
        bucket[index_dic[i]] += 1
    for i in range(k - 1):
        bucket[i + 1] += bucket[i]

    # sort left-most-S substrings
    sorted_lms = induced_sort(s, p, index_dic, bucket, k, t)

    # generate s1 by naming each left-most-S substrings by its bucket index
    sorted_lms = [i for i in sorted_lms if is_lms(t, i)]
    s1 = [-1] * len(s)
    num = 0
    s1[sorted_lms[0]] = num
    for m in range(len(sorted_lms) - 1):
        i = sorted_lms[m]
        j = sorted_lms[m + 1]
        dif = False
        for n in range(len(s)):
            if i + n >= len(s) or j + n >= len(s):
                break
            if s[i + n] != s[j + n] or is_lms(t, i + n) != is_lms(t, j + n):
                dif = True
                break
            elif n > 0 and (is_lms(t, i + n) or is_lms(t, j + n)):
                break
        if dif:
            num += 1
        s1[j] = num

    s1 = [i for i in s1 if i != -1]

    # generate sa1
    if num + 1 < len(s1):
        sa1 = gen_sa(s1, num + 1)
    else:
        # bucket sort s1 to get sa1
        sa1 = [0] * (num + 1)
        for i, name in enumerate(s1):
            sa1[name] = i

    return induced_sort(s, [p[i] for i in reversed(sa1)], index_dic, bucket, k, t)


def suffix_array(s):
    alphabet = sorted(list(set(s)))
    print('Alphabet: {}'.format(alphabet))
    return gen_sa(s, len(alphabet))


def generate_sa_CSV(sa, chr):
    BIN_LENGTH = 128
    sa = [str(i) for i in sa]
    index = 0
    output = []
    while index + BIN_LENGTH <= len(sa):
        output.append(','.join(sa[index:index + BIN_LENGTH]))
        output[-1] += '\n'
        index += BIN_LENGTH
    output.append(','.join(sa[index:]))
    with open('/data/student2022/zhangtk2022/Genome/{}.txt'.format(chr),'w') as f:
        f.writelines(output)


def gen_genome():
    log_record = [] # chromosome,length,time consumed
    genome = ''
    with open("/data/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa","r") as f:
        for line in f:
            if '>chr' in line:
                if len(genome) != 0:
                    genome += '$'
                    print(len(genome))
                    start = time.time()
                    sa = suffix_array(genome)
                    end = time.time()
                    generate_sa_CSV(sa,title)
                    log_record.append('{0},{1},{2}\n'.format(title,len(genome)-1,end-start))
                title = line.strip()
                genome = ''
            elif '>K' in line:
                break
            else:
                genome += line.strip()
    if len(genome) != 0:
        genome += '$'
        print(len(genome))
        start = time.time()
        sa = suffix_array(genome)
        end = time.time()
        generate_sa_CSV(sa,title)
        log_record.append('{0},{1},{2}'.format(title,len(genome)-1,end-start))
    # generate log file while nohup the exe
    with open('/home/zhangtk2022/homework/L7/log.txt','w') as f:
        f.writelines(log_record)
        
