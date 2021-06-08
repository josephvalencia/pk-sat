from pysmt.shortcuts import Symbol, And, Not, Or, Plus, Equals, Real,GE,LE,LT,GT, Implies, get_model, is_sat, get_formula_size, AtMostOne
from pysmt.typing import BV8, REAL, INT, ArrayType
import signal,time
from utils import parse_fasta, make_matrix, make_adjacency_dict,load_ground_truth,recall_precision,find_bps
from tqdm import tqdm

# every position can be in at most one bp over all pages
# eq 5 in IPKnot paper
def one_pairing(pages,pairs,probs,n):
    conditions = []
    m = len(pages)
    for i in range(n):
        unique = []
        for page in pages: # for each page
            for a,b in pairs:
                # i can be the left or right member
                if b == i or a == i:
                    unique.append(page[a][b])
        # only one bp allowed
        conditions.append(AtMostOne(unique))
    return conditions

# no internal pseudoknots
# eq 6
def no_internal_pks(pages,pairs,probs,n):
    conditions = []
    for page in pages:
        for i,j in pairs:
            for k,l in pairs:
                # this ordering indicates a pk, so both bp cannot be true
                if i < k and k < j and j < l:
                    conditions.append(Not(And(page[i][j],page[k][l])))
    return conditions

def no_bifurcations(pages,pairs,probs,n):
    p1 = pages[0]
    p2 = pages[1]
    conditions = []
    for i,j in pairs:
        for k,l in pairs:
            if j < k:
                conditions.append(Not(And(p2[i][j],p2[k][l])))
    return conditions

# every bp on pg2 must have a pk on pg1
# eq 7
def yes_external_pks(pages,pairs,probs,n):
    p1 = pages[0]
    p2 = pages[1]
    conditions = []
    for k,l in pairs:
        potential_pks = []
        for i,j in pairs:
            # (i,j) can be the leftmost bp or rightmost bp in pk
            if (i < k and k < j and j <l) or (k < i and i < l and l<j):
                potential_pks.append(p1[i][j])
        # one or more bp from other page introduces a pk
        if len(potential_pks) >0:
            is_pk = Implies(p2[k][l],Or(potential_pks))
            conditions.append(is_pk)
        # no bp from other page -> (i,j) not a bp
        else:
            conditions.append(Not(p2[k][l]))
    return conditions

# eq 8 and 9
def no_isolated_pairs(pages,pairs,probs,n):
    conditions = []
    for page in pages:
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

def score_constraints(pages,scores,pairs,probs,n):
    constraints = []
    # assign scores based on bp stacks
    for page,score in zip(pages,scores): 
        for i,j in pairs:
            # no need to check for neighbors, no_isolated_pairs() guarantees this already
            stack_score = Implies(page[i][j],Equals(score[i][j],Real(probs[i][j])))
            nostack_score = Implies(Not(page[i][j]),Equals(score[i][j],Real(0.0)))
            constraints.append(stack_score)
            constraints.append(nostack_score)
    return constraints

def score_threshold_constraint(scores,pairs,threshold):
    total = []
    for score in scores:
        for i,j in pairs:
            total.append(score[i][j])
    return GT(Plus(total),Real(threshold))

def predict_structure(seq,doubleknots=False):

    # legal pairs, lookup table
    pairs,probs = make_adjacency_dict(seq, len(seq))
    print("# pairs = {}".format(len(pairs)))

    # intialize pairing and score variables for two pages based on LinearPartition probs
    p1 = {}
    p2 = {}
    e1 = {}
    e2 =  {}
    for i,j in pairs:
        if i not in p1:
            p1[i] = {}
            p2[i] = {}
            e1[i] = {}
            e2[i] = {}
        x = Symbol('X_{},{}'.format(i,j))
        y = Symbol('Y_{},{}'.format(i,j))
        p1[i][j] = x
        p2[i][j] = y
        e = Symbol('E_{},{}'.format(i,j),REAL)
        E = Symbol('e_{},{}'.format(i,j),REAL)
        e1[i][j] = e
        e2[i][j] = E
   
    structural = structural_constraints([p1,p2],pairs,probs,len(seq))  
    score = score_constraints([p1,p2],[e1,e2],pairs,probs,len(seq))
    
    best_struct = ''.join(['.']*len(seq))
    best_score = 0
    threshold = 5
    verbose = False
    total_elapsed = 0
    MAX_TIME = 300

    signal.signal(signal.SIGALRM, timeout_handler)
    for it in range(10):
        s = time.time()
        thresh = score_threshold_constraint([e1,e2],pairs,threshold)
        all_constraints = structural+ score + [thresh]
        signal.alarm(MAX_TIME)
        try:
            solution = solve_SMT(all_constraints,[p1,p2],[e1,e2],pairs,seq)
        except Exception:
            break

        if solution is not None:
            structure,mea = solution
            best_struct = structure
            best_score = mea
        else:
            #print('No solution found')
            break
        
        elapsed = time.time() - s
        if verbose:
            print('_______________\nit # {} @  score thresh = {}'.format(it,threshold))
            print('Predicted  = {}'.format(structure))
            print('Expected accuracy score = {}'.format(mea))
            print('Time elapsed = {} seconds'.format(elapsed))
        
        total_elapsed += elapsed
        threshold *= 1.10

    print('***************\nBest Solution')
    print('Predicted = {}'.format(best_struct))
    print('Free Energy = {}'.format(best_score))
    print('Total time elapsed = {} seconds'.format(total_elapsed))
    print('***************')
    return best_struct

def timeout_handler(signum,frame):
    print('Exceeded max computation time')
    raise Exception('Timeout')

def solve_SMT(conditions,pages,scores,pairs,seq):

    # build SAT/SMT formula as conjunction of all terms
    formula = conditions[0]
    for c in conditions[1:]:
       formula = And(formula,c)

    # find satisfying assignment if one exists
    print('Finding assignment')
    model = get_model(formula,solver_name="z3")

    if model:
        value = get_structure(model,pages,scores,pairs,seq)
    else:
        value = None

    return value

def get_structure(model,pages,scores,pairs,seq):

    p1 = pages[0]
    p2 = pages[1]
    e1 = scores[0]
    e2 = scores[1]
    total_energy = 0
    structure = ['.'] * len(seq)
    chars_set = False
    for i,j in pairs:            
            energy_a = float(model.get_py_value(e1[i][j]))
            energy_b = float(model.get_py_value(e2[i][j]))
            if energy_a != 0.0:
                #print('E[{}][{}] = {}'.format(i,j,energy_a))
                pass
            if energy_b != 0.0:
                #print('e[{}][{}] = {}'.format(i,j,energy_b))
                pass
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
                #print('X[{}][{}] = {}'.format(i,j,bp_a))
            elif bp_b:
                if not chars_set:
                    page2_chars = ('(',')')
                    page1_chars = ('[',']')
                    chars_set = True
                # subtract 1 for 0-indexing
                structure[i-1] = page2_chars[0]
                structure[j-1] = page2_chars[1]
                #print('Y[{}][{}] = {}'.format(i,j,bp_b))

    structure = ''.join(structure)
    return structure,total_energy

def evaluate(pred,true):

    true_bp = find_bps(true)
    pred_bp = find_bps(pred)
    recall,precision = recall_precision(pred,true)
    print('Recall = {}, Precision = {}'.format(recall,precision))

if __name__ == '__main__':
 
    gt = load_ground_truth('all_PKB_structs.txt')

    for name,seq in parse_fasta('all_PKB.fa'): 
        pred = predict_structure(seq)
        true = gt[name]
        print(name)
        print(true)
        evaluate(pred,true)
