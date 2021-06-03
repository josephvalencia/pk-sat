import numpy.matlib
import numpy as np
import os

def get_matrix(f, l): # give it a filename and the sequence length
    f = open(f, 'r')
    l += 1 # because sequences index at 1, not 0
    m = np.matlib.zeros((l, l))
    for line in f.readlines():
        if len(line.split()): # ignore blank line
            (i, j, p) = [float(s) for s in line.split()] # get vals as floats
            (i, j)    = map(lambda x: int(x), (i, j))
            m[i, j] = p
    f.close()
    return m # this matrix is very sparse, don't be surprised if it's hard to find non-zero values

# give a sequence (as a filename), interface with LinearPartition to generate a readable file for it
def make_matrix_file(f):
    os.system(f'cat {f} | ../LinearPartition/linearpartition --prefix ./prob_matrix-{f}')


# combine above 2 functions (hacky patch with file creation and deletion)
def make_matrix(f, l):
    make_matrix_file(f)
    f = f'prob_matrix-{f}_1'
    m = get_matrix(f, l)
    os.system(f'rm {f}')
    return m


if __name__ == '__main__':
    m = make_matrix('CRW_00774.fa', 989)
    print(m)
