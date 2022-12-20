#!/usr/bin/python3
class BWT:
    GAP_LENGTH = 32
    def __init__(self, *args):
        if len(args) == 1:
            path = args[0]
            self.__decompression(path)
            self.alphabet = sorted(list(set(self.bw)))
        elif len(args) == 2:
            self.s = args[0]
            self.suffix_array = args[1]
            self.alphabet = sorted(list(set(self.s)))
            self.__BWT_encoder()
        self.__gen_bucket_fm()

    # encode BW string
    def __BWT_encoder(self):
        bw = []
        for i in self.suffix_array:
            if i == 0:
                bw.append("$")
            else:
                bw.append(self.s[i - 1])
        self.bw =  ''.join(bw)
    
    # gnerate bucket and fm array (the occurence)
    def __gen_bucket_fm(self):
        self.bucket = {key:[0,0] for key in self.alphabet}
        self.occrs = {key:[] for key in self.alphabet}
        index = 0
        while index < len(self.bw):
            self.bucket[self.bw[index]][1] += 1
            if index % self.GAP_LENGTH == 0:
                for ch in self.alphabet:
                    self.occrs[ch].append(self.bucket[ch][1])
            index += 1
        for ch in self.alphabet:
            if ch != '$':
                self.bucket[ch][1] += self.bucket[last_ch][1]
                self.bucket[ch][0] = self.bucket[last_ch][1]
            last_ch = ch        

    def __sa_at(self, index):
        count = 0
        while self.suffix_array[index] == -1:
            ch = self.bw[index]
            index = self.bucket[ch][0] - 1 + self.__count_occrs(ch, index)
            count += 1
        return self.suffix_array[index] + count

    # count occurence
    def __count_occrs(self, ch, index):
        if index < 0:
            return 0
        n = len(self.bw) // self.GAP_LENGTH 
        start = index // self.GAP_LENGTH
        rm = index % self.GAP_LENGTH
        if rm == 0:
            return self.occrs[ch][start]
        count = 0
        if rm >= self.GAP_LENGTH // 2 and start < n - 1:
            start += 1 
            for i in range(self.GAP_LENGTH - rm): 
                if self.bw[self.GAP_LENGTH*start - i] == ch: 
                    count -= 1
        else:
            for i in range(1, rm + 1):
                if self.bw[self.GAP_LENGTH*start + i] == ch:
                    count += 1
        return self.occrs[ch][start] + count

    # return matched substring
    def exact_match(self, query):
        index = len(query) - 1 
        ch = query[index] 
        lo = self.bucket[ch][0] 
        hi = self.bucket[ch][1] - 1 
        while lo <= hi and index > 0:
            ch = query[index - 1] 
            lo = self.bucket[ch][0] + self.__count_occrs(ch, lo - 1)
            hi = self.bucket[ch][0] - 1 + self.__count_occrs(ch, hi)
            index -= 1 
        return [self.__sa_at(i) for i in range(lo, hi + 1)] if lo <= hi else []

    def compression(self, chr):
        rl = self.__RL_encoder()
        BIN_LENGTH = 128
        rl = [str(i) for i in rl]
        index = 0
        output = []
        while index + BIN_LENGTH <= len(rl):
            output.append(','.join(rl[index:index + BIN_LENGTH]))
            output[-1] += '\n'
            index += BIN_LENGTH
        output.append(','.join(rl[index:]))
        with open('/home/zhangtk2022/homework/L7/Genome/{}/bw.txt'.format(chr),'w') as f:
            f.writelines(output)
    
    def __decompression(self, chr):
        csv = BWT.read_csv('/home/zhangtk2022/homework/L7/Genome/{}/bw.txt'.format(chr))
        self.bw = BWT.MTF_decoder(BWT.RL_decoder(csv), ['$', 'A', 'C', 'G', 'N', 'T'])
        sa = BWT.read_csv('/home/zhangtk2022/homework/L7/Genome/{}/sa.txt'.format(chr))
        self.suffix_array = [-1] * len(self.bw)
        for i in range(0, len(sa), 2):
            self.suffix_array[sa[i + 1]] = sa[i]
    
        # space compression
    def __MTF_encoder(self):
        mtf_list = []
        alphabet = self.alphabet
        for i in range(len(self.bw)):
            ch = self.bw[i]
            for j in range(len(alphabet)):
                if alphabet[j] == ch:
                    mtf_list.append(j)
                    index = j
                    while index > 0:
                        alphabet[index] = alphabet[index - 1]
                        index -= 1
                    alphabet[0] = ch
        return mtf_list

    def __RL_encoder(self):
        mtf_list = self.__MTF_encoder()
        rl_list = []
        i = 0
        while i < len(mtf_list):
            symbol = mtf_list[i]
            j = 0
            while i + j < len(mtf_list) and mtf_list[i+j] == symbol:
                j += 1
            rl_list += [symbol, j]
            i += j
        return rl_list

    @staticmethod
    def read_csv(path):
        csv = []
        with open(path, 'r') as f:
            for line in f:
                csv += line.strip().split(',')
        return [int(i) for i in csv]

    @staticmethod
    def RL_decoder(rl):
        decoded_rl = []
        for i in range(0,len(rl),2):
            decoded_rl += [rl[i]] * rl[i + 1]
        return decoded_rl

    @staticmethod
    def MTF_decoder(mtf, alphabet):
        s = []
        for i in mtf:
            ch = alphabet[i]
            s.append(ch)
            j = i
            while j > 0:
                alphabet[j] = alphabet[j - 1]
                j -= 1
            alphabet[0] = ch
        return ''.join(s)
