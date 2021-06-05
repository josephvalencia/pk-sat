from pysmt.shortcuts import Symbol, And, Not, Plus, Equals, Real,GE,LE,LT, Implies, get_model, is_sat, get_formula_size
from pysmt.typing import BV8, REAL, INT, ArrayType
import time


# every position can be in at most one bp over all pages
# eq 5 in IPKnot paper
def one_pairing(pages):
    conditions = []
    m = len(pages)
    for p in range(m): # for each page
        page = pages[p]
        (rows, cols) = page.shape()
        for i in range(rows):
            for j in range(i + 1, cols): # for every cell in adjacency matrix
                or_cond = False
                for q in range(p + 1, m): # for every other page
                    qage = pages[q]
                    or_cond = Or(or_con, qage[i][j]) # build up a bil ol' disjuction
                conditions.append(Implies(page[i][j], Not(or_cond)))
                # if [i][j] is true in page p, then it is not true in any other page
    return conditions


# no internal pseudoknots
# eq 6
def no_internal_pks(pages):
    conditions = []
    for page in pages:
        for i in range(len(seq)):
            for k in range(len(seq)):
                for j in range(len(seq)):
                    for l in range(len(seq)):
                        if i < k and k < j and j < l:
                            conditions.append(Not(And(page[i][j],page[k][l])))
    return conditions


def yes_external_pks(pages):
    conditions = []
    m = len(pages)
    for p in range(m): # for all pages p
        page = pages[p]
        (rows, cols) = page.shape()
        for q in range(0, p): # for all pages lower than p
            qage = pages[q]
            # ensure pseudoknots (not sure about this logic)
            for i in range(rows): # for all i
                for k in range(i + 1, rows): # for all k > i
                    for j in range(k + 1, cols): # for all j > k
                        for l in range(j + 1, cols): # for all l > j
                            conditions.append(Implies(page[i][j], qage(k, l))) # THIS IS THE LOGIC I'M UNSURE OF
    return conditions


def structural_constraints(pages):
    return one_pairing(pages) + no_internal_pks(pages) + yes_external_pks(pages)


ig __name__ == '__main__':
    print(one_pairing('CG', ))
