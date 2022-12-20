#!/usr/bin/python3
class MaxPQ:
    def __init__(self, *args):
        if len(args) == 0:
            self.pq = [None]
            self.size = 0
        elif len(args) == 1:
            self.pq = [None] + args[0]
            self.size = len(self.pq) - 1
            for i in range(self.size // 2, 0, -1):
                self.__sink(i)


    def __exch(self, i, j):
        node = self.pq[i]
        self.pq[i] = self.pq[j]
        self.pq[j] = node

    def __swim(self, i):
        while i > 1 and self.pq[i // 2] < self.pq[i]:
            self.__exch(i // 2, i)
            i = i // 2

    def __sink(self, i):
        while 2 * i <= self.size:
            j = 2 * i
            if j < self.size and self.pq[j] < self.pq[j + 1]:
                j += 1
            if not self.pq[i] < self.pq[j]:
                break
            self.__exch(i, j)
            i = j

    def is_empty(self):
        return self.size == 0

    def insert(self, node):
        self.size += 1
        self.pq.append(node)
        self.__swim(self.size)

    def pop_max(self):
        node = self.pq[1]
        self.__exch(1, self.size)
        self.pq[self.size] = None
        self.size -= 1
        self.__sink(1)
        return node
