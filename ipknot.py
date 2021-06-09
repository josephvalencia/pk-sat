from pysmt.shortcuts import Symbol, And, Not, Or, Plus, Equals, Real,GE,LE,LT,GT, Implies, get_model, is_sat, get_formula_size, AtMostOne
from pysmt.typing import BV8, REAL, INT, ArrayType
import signal,time
from utils import parse_fasta,load_ground_truth, make_LinearPartition_files, make_matrix, make_adjacency_dict,reverse_engineer_mea,cleanup_files,recall_precision,find_bps,to_dot_bracket
from tqdm import tqdm

class TimeoutException(Exception):
    pass

# every position can be in at most one bp over all pages
# eq 5 in IPKnot paper
def one_pairing(page_set,pair_set,probs,n):
    conditions = []
    for i in range(n):
        unique = []
        for page,pairs in zip(page_set,pair_set): # for each page
            for a,b in pairs:
                # i can be the left or right member
                if b == i or a == i:
                    unique.append(page[a][b])
        # only one bp allowed
        conditions.append(AtMostOne(unique))
    return conditions

# no internal pseudoknots
# eq 6
def no_internal_pks(page_set,pair_set,probs,n):
    conditions = []
    for page,pairs in zip(page_set,pair_set):
        for i,j in pairs:
            for k,l in pairs:
                # this ordering indicates a pk, so both bp cannot be true
                if i < k and k < j and j < l:
                    conditions.append(Not(And(page[i][j],page[k][l])))
    return conditions

# as complex PKs should be rare, bifurcations are not allowed on pg2
def no_bifurcations(page_set,pair_set,probs,n):
    p1 = page_set[0]
    p2 = page_set[1]
    p1_pairs = pair_set[0]
    p2_pairs = pair_set[1]
    conditions = []
    for i,j in p2_pairs:
        for k,l in p2_pairs:
            if j < k:
                conditions.append(Not(And(p2[i][j],p2[k][l])))
    return conditions

# every bp on pg2 must have a pk on pg1
# eq 7
def yes_external_pks(page_set,pair_set,probs,n):
    p1 = page_set[0]
    p2 = page_set[1]
    p1_pairs = pair_set[0]
    p2_pairs = pair_set[1]
    conditions = []
    for k,l in p2_pairs:
        potential_pks = []
        for i,j in p1_pairs:
            # (i,j) can be the leftmost bp or rightmost bp in pk
            if (i < k and k < j and j <l) or (k < i and i < l and l<j):
                potential_pks.append(p1[i][j])
        # one or more bp from other page introduces a pk
        if len(potential_pks) > 0:
            is_pk = Implies(p2[k][l],Or(potential_pks))
            conditions.append(is_pk)
        # no bp from other page -> (i,j) not a bp
        else:
            conditions.append(Not(p2[k][l]))
    return conditions

# eq 8 and 9
def no_isolated_pairs(page_set,pair_set,probs,n):
    conditions = []
    for p in range(len(page_set)):
        page = page_set[p]
        pairs = pair_set[p]
        for i,j in pairs:
            neighbors = []
            #outside neighbor
            if i-1 in page and j+1 in page[i-1]:
                neighbors.append(page[i-1][j+1])
            # inside neigbor
            if i+1 in page and j-1 in page[i+1]:
                neighbors.append(page[i+1][j-1])
            # if (i,j) is a pair, (i-1,j+1) or (i+1,j-1) must be a pair
            if len(neighbors) > 0:
                stack_required = Implies(page[i][j],Or(neighbors))
                conditions.append(stack_required)
            # no neighbor means this is an illegal isolated` bp
            else:
                conditions.append(Not(page[i][j]))
    return conditions

def structural_constraints(pages,pairs,probs,n):
    
    total=one_pairing(pages,pairs,probs,n) 
    total += yes_external_pks(pages,pairs,probs,n)
    total+=no_internal_pks(pages,pairs,probs,n)
    total+= no_isolated_pairs(pages,pairs,probs,n)
    total+=no_bifurcations(pages,pairs,probs,n)
    return total

def score_constraints(page_set,score_set,pair_set,probs,n):
    constraints = []
    # assign scores based on bp stacks
    for p in range(len(page_set)):
        page = page_set[p]
        score = score_set[p]
        pairs = pair_set[p]
        for i,j in pairs:
            # no need to check for neighbors, no_isolated_pairs() guarantees this already
            stack_score = Implies(page[i][j],Equals(score[i][j],Real(probs[i][j])))
            nostack_score = Implies(Not(page[i][j]),Equals(score[i][j],Real(0.0)))
            constraints.append(stack_score)
            constraints.append(nostack_score)
    return constraints

def score_threshold_constraint(score_set,pair_set,threshold):
    total = []
    for p in range(len(score_set)):
        score = score_set[p]
        pairs = pair_set[p]
        for i,j in pairs:
            total.append(score[i][j])
    return GT(Plus(total),Real(threshold))

def predict_structure(seq,name,k1,k2):

    # legal pairs, lookup table
    make_LinearPartition_files(seq,name)

    theta_p1 = 1/(2**k1+1)
    theta_p2 = 1/(2**k2+1) 
    
    p1_pairs,probs = make_adjacency_dict(seq,name, len(seq),theta_p1)
    mea_struct,mea_score = reverse_engineer_mea(seq,name)
    print('Starting threshold = {}'.format(mea_score))
    cleanup_files(name)
    p2_pairs = [(a,b) for a,b in p1_pairs if probs[a][b] > theta_p2]
    
    # intialize pairing and score variables for two pages based on LinearPartition probs
    p1 = {}
    p2 = {}
    e1 = {}
    e2 =  {}
    for i,j in p1_pairs:
        if i not in p1:
            p1[i] = {}
            e1[i] = {}
        x = Symbol('X_{},{}'.format(i,j))
        p1[i][j] = x
        E = Symbol('E_{},{}'.format(i,j),REAL)
        e1[i][j] = E
   
    for i,j in p2_pairs:
        if i not in p2:
            p2[i] = {}
            e2[i] = {}
        y = Symbol('Y_{},{}'.format(i,j))
        p2[i][j] = y
        e = Symbol('e_{},{}'.format(i,j),REAL)
        e2[i][j] = e
  
    all_pages = [p1,p2]
    all_scores = [e1,e2]
    all_pairs = [p1_pairs,p2_pairs]
    structural = structural_constraints(all_pages,all_pairs,probs,len(seq))  
    score = score_constraints(all_pages,all_scores,all_pairs,probs,len(seq))
    
    best_struct = mea_struct
    best_score = 0
    growth_rate = 1.001
    threshold = 0.95*mea_score
    verbose = False
    total_elapsed = 0
    MAX_TIME = 600

    #signal.signal(signal.SIGALRM, timeout_handler)
    num_iter = 0 
    for it in range(50):
        s = time.time()
        thresh = score_threshold_constraint(all_scores,all_pairs,threshold)
        all_constraints = structural+ score + [thresh]
        signal.alarm(MAX_TIME)
        try:
            solution = solve_SMT(all_constraints,all_pages,all_scores,all_pairs,seq)
            signal.alarm(0)
        except TimeoutException :
            print('Something went wrong')
            break

        elapsed = time.time() - s
        total_elapsed += elapsed
        threshold *= growth_rate
        num_iter +=1
        
        if solution is not None:
            structure,mea = solution
            best_struct = structure
            best_score = mea
        else:
            print('No solution found')
            break
        
        if verbose:
            print('_______________\nit # {} @  score thresh = {}'.format(it,threshold))
            print('Predicted  = {}'.format(structure))
            print('Expected accuracy score = {}'.format(mea))
            print('Time elapsed = {} seconds'.format(elapsed))
        
    summary = {'name' : name , 'seq' : str(seq) , 'len' : len(seq), 'MEA_score' : mea_score, 'MEA_struct' : mea_struct,'time' : total_elapsed, \
            'iterations' : num_iter ,  'k1' : k1, 'k2' : k2, 'pairs1' : len(p1_pairs) , 'pairs2' : len(p2_pairs),'pred_struct' : best_struct , 'pred_score' : best_score}
    return summary

def timeout_handler(signum,frame):
    print('Exceeded max computation time')
    raise TimeoutException('Timeout')

def solve_SMT(conditions,pages,scores,pairs,seq):

    # build SAT/SMT formula as conjunction of all terms
    formula = conditions[0]
    for c in conditions[1:]:
       formula = And(formula,c)

    # find satisfying assignment if one exists
    print('Finding assignment')
    model = get_model(formula,solver_name="z3")

    if model:
        print('Found!!!!!')
        value = get_structure(model,pages,scores,pairs,seq)
    else:
        value = None

    return value

def get_structure(model,page_set,score_set,pair_set,seq):

    p1 = page_set[0]
    p2 = page_set[1]
    e1 = score_set[0]
    e2 = score_set[1]
    total_energy = 0
    structure = ['.'] * len(seq)
    chars_set = False
    '''
    for i,j in pairs:            
            # sum scores
            energy_a = float(model.get_py_value(e1[i][j]))
            energy_b = float(model.get_py_value(e2[i][j]))
            total_energy += energy_a
            total_energy += energy_b
            # assemble dot-bracket
            bp_a = model.get_py_value(p1[i][j])
            bp_b = model.get_py_value(p2[i][j])
            if bp_a:
                # page variables arbitrarily but left-most output page given rounded brackets
                if not chars_set:
                    page1_chars = ('(',')')
                    page2_chars = ('[',']')
                    chars_set = True
                # subtract 1 for 0-indexing
                structure[i-1] = page1_chars[0]
                structure[j-1] = page1_chars[1]
            elif bp_b:
                if not chars_set:
                    page2_chars = ('(',')')
                    page1_chars = ('[',']')
                    chars_set = True
                # subtract 1 for 0-indexing
                structure[i-1] = page2_chars[0]
                structure[j-1] = page2_chars[1]
    '''

    page_char_set  = [('(',')'),('[',']')]

    for page,score,pairs,chars in zip(page_set,score_set,pair_set,page_char_set): 
        for i,j in pairs:            
            # sum scores
            energy = float(model.get_py_value(score[i][j]))
            total_energy += energy
            # assemble dot-bracket
            is_bp = model.get_py_value(page[i][j])
            if is_bp:
                # subtract 1 for 0-indexing
                structure[i-1] = chars[0]
                structure[j-1] = chars[1]
    
    structure = ''.join(structure)
    return structure,total_energy

def add_evaluation_metrics(summary,truth):
    
    pred = summary['pred_struct']
    linear_partition_pred = summary['MEA_struct']
    
    # score our prediction
    truth = to_dot_bracket(truth)
    true_bp = find_bps(truth)
    pred_bp = find_bps(pred)
    recall,precision = recall_precision(pred_bp,true_bp)
   
    # score the PK-free prediction from LinearPartition
    lp_pred_bp = find_bps(linear_partition_pred)
    lp_recall,lp_precision = recall_precision(lp_pred_bp,true_bp)

    # append metrics to summary dict
    summary['recall'] = recall
    summary['precision'] = precision
    f1 = 0 if precision + recall == 0 else 2*precision*recall/(precision+recall)
    summary['f1'] = f1 
    summary['lp_recall'] = lp_recall
    summary['lp_precision'] = lp_precision
    lp_f1 = 0 if lp_precision + lp_recall == 0 else 2*lp_precision*lp_recall/(lp_precision+lp_recall)
    summary['lp_f1'] = lp_f1
    return summary

if __name__ == '__main__':
 
    gt = load_ground_truth('all_PKB_structs.txt')
    with open('trials_tmp.csv','w') as outFile:
        header = ['name' , 'seq' , 'len', 'MEA_score' , 'MEA_struct' ,'time' , 'k1', 'k2', 'pairs1','pairs2','iterations' ,'pred_struct' ,\
                'pred_score' ,'recall','precision','f1','lp_recall','lp_precision','lp_f1']
        outFile.write(','.join(header)+'\n')
        for k1 in range(2,3):
            for k2 in range(k1,0,-1):
                print(k1,k2)
                for name,seq in parse_fasta('all_PKB.fa'): 
                    truth = gt[name]
                    summary  = predict_structure(seq,name,k1,k2)
                    summary = add_evaluation_metrics(summary,truth)
                    line = [str(summary[x]) for x in header]
                    outFile.write(','.join(line)+'\n')
