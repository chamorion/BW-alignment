#!/usr/bin/python3
class AlgnTree:
    def __init__(self, algn):
        self.algn = algn
        self.le = None
        self.ri = None
    
    def insert(self, algn_tree):
        if algn_tree.algn.section_ref[0] < self.algn.section_ref[0]:
            if self.le is None:
                self.le = algn_tree
                return
            self.le.insert(algn_tree)
        else:
            if self.algn.can_merge_with(algn_tree.algn):
                self.algn = self.algn.merge_with(algn_tree.algn)
                return
            if self.ri is None:
                self.ri = algn_tree
                return
            if self.algn.section_ref[0] < algn_tree.algn.section_ref[0]:
                self.ri.insert(algn_tree)
            else:
                algn_tree.ri = self.ri
                self.ri = algn_tree
    
    def __pop_leaves(self, algn_list):
        if self.le is not None:
            algn_list = self.le.__pop_leaves(algn_list)
        if self.ri is not None:
            algn_list = self.ri.__pop_leaves(algn_list)
        algn_list.append(self.algn)
        return algn_list

    def pop_all(self):
        return self.__pop_leaves([])
