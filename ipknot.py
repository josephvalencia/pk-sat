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
                    or_cond = Or(or_con, qage[i][j])
                conditions.append(Implies(page[i][j], Not(or_cond)))
                # if [i][j] is true in page p, then it is not true in any other page
    return conditions



ig __name__ == '__main__':
    print(one_pairing('CG', ))
