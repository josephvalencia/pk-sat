import numpy.matlib as np
import os
from Bio import SeqIO
from collections import defaultdict

def get_matrix(f, l): # give it a filename and the sequence length
    f = open(f, 'r')
    l += 1 # because sequences index at 1, not 0
    m = np.zeros((l, l))
    for line in f.readlines():
        if len(line.split()): # ignore blank line
            (i, j, p) = [float(s) for s in line.split()] # get vals as floats
            (i, j)    = map(lambda x: int(x), (i, j))
            m[i, j] = p
    f.close()
    return m # this matrix is very sparse, don't be surprised if it's hard to find non-zero values

def get_adjacency_dict(f,l):
    
    k = 3
    MIN_THRESH = 1/(2**k+1) 
    print("MIN_THRESH = {}".format(MIN_THRESH))
    pairs = []
    storage = defaultdict(lambda : defaultdict(int))
    with open(f) as inFile:
        for line in inFile.readlines():
            if len(line.split()): # ignore blank line
                (i, j, p) = [float(s) for s in line.split()] # get vals as floats
                (i, j)    = map(lambda x: int(x), (i, j))
                if p > MIN_THRESH:
                    storage[i][j] = p
                    pairs.append((i,j))
    
    return pairs,storage

# give a sequence (as a filename), interface with LinearPartition to generate a readable file for it
def make_matrix_file(f):
    os.system(f'echo {f} | ../LinearPartition/linearpartition --prefix ./prob_matrix-tmp')

# combine above 2 functions (hacky patch with file creation and deletion)
def make_matrix(f, l):
    make_matrix_file(f)
    f = f'prob_matrix-{f}_1'
    m = get_matrix(f, l)
    os.system(f'rm {f}')
    return m

def make_adjacency_dict(f, l):
    make_matrix_file(f)
    f = f'prob_matrix-tmp_1'
    m = get_adjacency_dict(f, l)
    os.system(f'rm {f}')
    return m

def parse_fasta(f):
    with open(f) as inFile:
        for record in SeqIO.parse(inFile,'fasta'):
            yield record.seq

def to_uppercase(seq):
    return seq.upper()

def to_dot_bracket(struct):
    return struct.replace(':','.')

def find_bps(struct):

    bp_list = []
    left_paren_locs = []
    left_brack_locs = []

    # find left members
    for i in range(len(struct)):
        c = struct[i]
        if c == '(':
            left_paren_locs.append(i)
        elif c == '[':
            left_brack_locs.append(i)
    # find right paren matches
    for p in left_paren_locs:
        count = 1
        for j in range(p+1,len(struct)):
            c = struct[j]
            if c == '(':
                count +=1
            elif c == ')':
                count -=1
            if count == 0:
                # plus one for 1-base
                bp_list.append((p+1,j+1))
                break
    # find right brack matches
    for p in left_brack_locs:
        count = 1
        for j in range(p+1,len(struct)):
            c = struct[j]
            if c == '[':
                count +=1
            elif c == ']':
                count -=1
            if count == 0:
                # plus one for 1-base
                bp_list.append((p+1,j+1))
                break 
    return bp_list

def recall_precision(pred,true):

    pred_set = set(pred)
    true_set = set(true)
    TP = 0
    FP = 0
    FN = 0
    
    for p in pred_set:
        if p in true_set:
            TP+=1
            true_set.remove(p)
        else:
            FP+=1
    
    FN = len(true_set)
    recall = 0 if TP+FN == 0 else TP/(TP+FN)
    precision = 0 if TP+FP == 0 else TP/(TP+FP)
    return recall,precision

# FOR DEBUGGING ONLY
def reverse_engineer_partition(seq,true_struct):
    
    make_matrix_file(seq)
    f = f'prob_matrix-tmp_1'
    pairs = find_bps(true_struct)
    log = "log_tmp.txt" 
    print(pairs)
    for s,e in pairs:
        os.system(f'grep \"{s} {e}\" {f} >> {log}')
    
    print('Extracted partition function')
    os.system('awk \'{s+=$3}END{print s}\''+f' {log}') 
    os.system(f'rm {log}')
    os.system(f'rm {f}')

if __name__ == '__main__':

    pkb147 = "AUAAUAGAAUAGGACGUUUGGUUCUAUUUUGUUGGUUUCUAGGACCAUCGU"
    pkb147_struct = to_dot_bracket(":((((((((((:::[[::[[[[[[)))))))))):::::::]]]]]]:]]:")
    reverse_engineer_partition(pkb147,pkb147_struct)
