import numpy as np
import pandas as pd
import os,subprocess
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

def get_adjacency_dict(f,l,theta):

    pairs = []
    storage = defaultdict(lambda : defaultdict(int))
    with open(f) as inFile:
        for line in inFile.readlines():
            if len(line.split()): # ignore blank line
                (i, j, p) = [float(s) for s in line.split()] # get vals as floats
                (i, j)    = map(lambda x: int(x), (i, j))
                if p > theta:
                    storage[i][j] = p
                    pairs.append((i,j))

    return pairs,storage

# give a sequence, interface with LinearPartition to generate readable files of bp probs and PK-free MEA structure
def make_LinearPartition_files(f,name):
    cmd = f'echo {f} | ../LinearPartition/linearpartition --prefix ./prob_matrix_{name}_tmp -M --mea_prefix ./mea_struct_{name}_tmp'
    print(cmd)
    os.system(cmd)

# combine above 2 functions (hacky patch with file creation and deletion)
def make_matrix(f,name, l):
    f = f'prob_matrix_{name}_tmp_1'
    m = get_matrix(f, l)
    return m

def make_adjacency_dict(f,name, l,theta):
    f = f'prob_matrix_{name}_tmp_1'
    m = get_adjacency_dict(f, l,theta)
    return m

def parse_fasta(f):
    with open(f) as inFile:
        for record in SeqIO.parse(inFile,'fasta'):
            yield record.id,record.seq

def to_uppercase(seq):
    return seq.upper()

def to_dot_bracket(struct):
    return struct.replace(':','.')

def load_ground_truth(f):
    gt = {}
    names = []
    with open(f) as inFile:
        entries = inFile.read().split('>')[1:]
        for e in entries:
            l = e.lstrip().rstrip()
            fields = l.split('\n')
            name = fields[0]
            struct = ''.join(fields[1:])
            gt[name] = struct
            names.append(name)
    return gt

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

def reverse_engineer_mea(seq,name,struct=None):

    if struct is None:
        mea_file = f'mea_struct_{name}_tmp_1'
        with open(mea_file) as inFile:
            lines = inFile.readlines()
            mea_struct = lines[1].rstrip()
    else:
        mea_struct = struct

    probs_file = f'prob_matrix_{name}_tmp_1'
    log_file = f'log_{name}_tmp.txt'

    pairs = find_bps(mea_struct)
    # work backwards from MEA structure and sum all bp probs
    for s,e in pairs:
        os.system(f'grep \"{s} {e}\" {probs_file} >> {log_file}')
    result = subprocess.run(['awk','{s+=$3}END{print s}',f'{log_file}'],stdout=subprocess.PIPE)

    score = float(result.stdout.decode('utf-8'))
    return mea_struct,score

def safe_remove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def cleanup_files(name):
    safe_remove(f'prob_matrix_{name}_tmp_1')
    safe_remove(f'mea_struct_{name}_tmp_1')
    safe_remove(f'log_{name}_tmp.txt')

# create a string repr'ing a latex table from a csv file
def mkLaTexTable(fname, usecols=None):
    df = pd.read_csv(fname)
    return df.to_latex(index=False, columns=['name', 'len', 'pred_score', 'time', 'iterations','precision','recall','precision_improvement','recall_improvement'])

if __name__ == '__main__':
    #print(mkLaTexTable('trials.csv'))
    # print(load_ground_truth('all_PKB_structs.txt'))
    pkb115_seq = 'CGGUAGCGCGAACCUUAUCGCGCA'
    pkb115_struct = '.(((.[[[[[[)))...]]]]]].'
    #make_LinearPartition_files(pkb115_seq,'pkb115')
    #print(reverse_engineer_mea(pkb115_seq,'pkb115'))
    #cleanup_files('pkb115')
    make_LinearPartition_files(pkb115_seq,'pkb115')
    print(reverse_engineer_mea(pkb115_seq,'pkb115',struct=pkb115_struct))
    cleanup_files('pkb115')
