from pysmt.shortcuts import Symbol, And, Not, Plus, Equals, Real,GE,LE,LT, Implies, get_model, is_sat, get_formula_size
from pysmt.typing import BV8, REAL, INT, ArrayType
import time

def energy_score(a1,a2,b1,b2):

    key = a1+a2+b1+b2

    scores = {'AUAU': -1.1,'AUCG' : -2.1 ,'AUGC' : -2.2 ,'AUGU' : -1.4 , 'AUUG' :  -0.9, 'AUUA' :-0.6 ,\
            'CGAU' : -2.1 , 'CGCG' : -2.4 , 'CGGC' : -3.3 , 'CGGU' : -2.1, 'CGUG' : -2.1 , 'CGUA' : -1.4 ,\
            'GCAU' : -2.2 , 'GCCG' : -3.3 , 'GCGC' : -3.4 , 'GCGU': -2.5 , 'GCUG' : -2.4 , 'GCUA' : -1.5 ,\
            'GUAU' : -1.4 , 'GUCG' : -2.1 , 'GUGC' : -2.5 , 'GUGU' : -1.3, 'GUUG' : -1.3 , 'GUUA' :-0.5 ,\
            'UGAU' : -0.9, 'UGCG' : -2.1 , 'UGGC' : -2.4 , 'UGGU' : -1.3 , 'UGUG' :-1.3 , 'UGUA' : -1.0, \
            'UAAU' : -0.6, 'UACG' : -1.4 , 'UAGC' : -1.5 , 'UAGU' : -0.5 , 'UAUG' : -1.0, 'UAUA' : -0.3}

    if key in scores:
        return scores[key]
    else:
        # 'large' value
        return 10000.0

def predict_structure(seq, doubleknots=False):

    p1 = {}
    p2 = {}
    e1 = {}
    e2 =  {}

    # intialize n*n pairing variables for two pages
    for i in range(len(seq)):
        p1[i] = {}
        p2[i] = {}
        for j in range(len(seq)):
            if i < j :
                x = Symbol('X_{},{}'.format(i,j))
                y = Symbol('Y_{},{}'.format(i,j))
                p1[i][j] = x
                p2[i][j] = y

    # initialize n*n energy scores for two pages
    for i in range(len(seq)):
        e1[i] = {}
        e2[i] = {}
        for j in range(len(seq)):
            if i < j:
                e = Symbol('E_{},{}'.format(i,j),REAL)
                E = Symbol('e_{},{}'.format(i,j),REAL)
                e1[i][j] = e
                e2[i][j] = E

    structural_constraints = build_structural_constraints(seq,p1,p2, doubleknots=doubleknots)
    energy_constraints = build_energy_constraints(seq,p1,p2,e1,e2)
    all_constraints = structural_constraints+energy_constraints

    min_fe = 0
    max_fe = (len(seq) // 2) * -3.4 # best possible configuration is 100% nesting of CG
    print('Max Theoretical FE = {}'.format(max_fe))
    print('# Clauses = {}'.format(len(all_constraints)+1))
    best_struct = ''.join(['.']*len(seq))
    best_fe = 0
    threshold = 0

    total_elapsed = 0
    for it in range(10):
        s = time.time()
        prev_threshold = threshold
        threshold = (min_fe+max_fe) / 2
        step_size = threshold - prev_threshold
        if abs(step_size) < 0.1:
            print('Binary search converged, stopping early!')
            break
        thresh_constraint = energy_threshold_constraint(seq,e1,e2,threshold)
        print('_______________\nit # {} @  FE thresh = {}'.format(it,threshold))
        solution = solve_SMT(all_constraints+[thresh_constraint],p1,p2,e1,e2,seq)
        if solution is not None:
            structure,free_energy = solution
            print('Predicted  = {}'.format(structure))
            print('Free Energy = {}'.format(free_energy))
            min_fe = free_energy
            improvement = best_fe - free_energy
            best_struct = structure
            best_fe = free_energy
        else:
            print('No solution found')
            max_fe = threshold

        elapsed = time.time() - s
        print('Time elapsed = {} seconds'.format(elapsed))
        total_elapsed += elapsed

    print('***************\nBest Solution')
    print('Sequence = {}'.format(seq))
    print('Predicted = {}'.format(best_struct))
    print('Free Energy = {}'.format(best_fe))
    print('Total time elapsed = {} seconds'.format(total_elapsed))
    print('***************')

def build_energy_constraints(seq,p1,p2,e1,e2):

    constraints = []
    MIN_BP_DIST = 3
    # assign scores based on bp stacks
    for i in range(len(seq)):
        for j in range(len(seq)):
            if i < j:
                if j-i >= MIN_BP_DIST:
                    nuc1 = seq[i]
                    nuc2 = seq[j]
                    nuc3 = seq[i+1]
                    nuc4 = seq[j-1]
                    bp_stack_energy = energy_score(nuc1,nuc2,nuc3,nuc4)
                    e1_stack_score = Implies(And(p1[i][j],p1[i+1][j-1]),Equals(e1[i][j],Real(bp_stack_energy)))
                    e2_stack_score = Implies(And(p2[i][j],p2[i+1][j-1]),Equals(e2[i][j],Real(bp_stack_energy)))
                    e1_nostack_score = Implies(Not(And(p1[i][j],p1[i+1][j-1])),Equals(e1[i][j],Real(0.0)))
                    e2_nostack_score = Implies(Not(And(p2[i][j],p2[i+1][j-1])),Equals(e2[i][j],Real(0.0)))
                    new_constraints = [e1_stack_score,e2_stack_score,e1_nostack_score,e2_nostack_score]
                    constraints.extend(new_constraints)
                else:
                    e1_not_beginning = Equals(e1[i][j],Real(0.0))
                    e2_not_beginning = Equals(e2[i][j],Real(0.0))
                    constraints.append(e1_not_beginning)
                    constraints.append(e2_not_beginning)

    # all bp stacks accounted for
    for i in range(1,len(seq)-1):
        for j in range(1,len(seq)-1):
            if i < j and j-i >= MIN_BP_DIST:
                p1_constraint = Not(And(Not(p1[i-1][j+1]),p1[i][j],Not(p1[i+1][j-1])))
                p2_constraint = Not(And(Not(p2[i-1][j+1]),p2[i][j],Not(p2[i+1][j-1])))
                constraints.append(p1_constraint)
                constraints.append(p2_constraint)

    return constraints

def energy_threshold_constraint(seq,e1,e2,threshold):

    # check if total energy threshold is exceeded
    total = []

    for i in range(len(seq)):
        for j in range(len(seq)):
            if i < j:
                total.append(e1[i][j])
                total.append(e2[i][j])

    free_energy = LT(Plus(total),Real(threshold))
    return free_energy

def build_structural_constraints(seq,p1,p2, doubleknots=False):

    conditions = []

    # every position can be in at most one bp
    for i in range(len(seq)):
        for j in range(len(seq)):
            for k in range(len(seq)):
                if i < j:
                    if j < k:
                        p1_internal = Not(And(p1[i][j],p1[j][k]))
                        p2_internal = Not(And(p2[i][j],p2[j][k]))
                        p1_first = Not(And(p1[i][j],p2[j][k]))
                        p2_first = Not(And(p2[i][j],p1[j][k]))
                        new_conditions = [p1_internal,p2_internal,p1_first,p2_first]
                        conditions.extend(new_conditions)
                    elif k < j and i != k:
                        p1_internal = Not(And(p1[i][j],p1[k][j]))
                        p2_internal = Not(And(p2[i][j],p2[k][j]))
                        p1_first = Not(And(p1[i][j],p2[k][j]))
                        p2_first = Not(And(p2[i][j],p1[k][j]))
                        new_conditions = [p1_internal,p2_internal,p1_first,p2_first]
                        conditions.extend(new_conditions)

    # bp must belong to either page1 or page2
    for i in range(len(seq)):
        for j in range(len(seq)):
            if i < j:
                unique_page = Not(And(p1[i][j],p2[i][j]))
                conditions.append(unique_page)

    # no internal pseudoknots
    for i in range(len(seq)):
        for k in range(len(seq)):
            for j in range(len(seq)):
                for l in range(len(seq)):
                    if i < k and k < j and j < l:
                        p1_pk = Not(And(p1[i][j],p1[k][l]))
                        p2_pk = Not(And(p2[i][j],p2[k][l]))
                        conditions.append(p1_pk)
                        conditions.append(p2_pk)
    # only nested bp (no splits) on page2
    for i in range(len(seq)):
        for j in range(len(seq)):
            for k in range(len(seq)):
                for l in range(len(seq)):
                    if i < j and j < k and k < l:
                        p2_split = Not(And(p2[i][j],p2[k][l]))
                        conditions.append(p2_split)

    if doubleknots:
        # no double-crossing pseudoknots
        for i in range(len(seq)):
            for m in range(len(seq)):
                for j in range(len(seq)):
                    for k in range(len(seq)):
                        for n in range(len(seq)):
                            for l in range(len(seq)):
                                if i < m and m < j and j < k and k < n and n < l:
                                    downstream = Implies(And(p1[i][j],p2[m][n]),Not(p1[k][l]))
                                    upstream = Implies(And(p1[k][l],p2[m][n]),Not(p1[i][j]))
                                    conditions.append(downstream)
                                    conditions.append(upstream)


    return conditions

def solve_SMT(conditions,p1,p2,e1,e2,seq):

    # build SAT/SMT formula as conjunction of all terms
    formula = conditions[0]
    for c in conditions[1:]:
       formula = And(formula,c)

    # find satisfying assignment if one exists
    print('Finding assignment')
    model = get_model(formula,solver_name="z3")

    if model:
        value = get_structure(model,p1,p2,e1,e2,seq)
    else:
        value = None

    return value

def get_structure(model,p1,p2,e1,e2,seq):

    total_energy = 0
    structure = ['.'] * len(seq)
    print('Solution found!')
    chars_set = False
    for i in range(len(seq)-1):
        for j in range(i+1,len(seq)):
            # sum free energies
            energy_a = float(model.get_py_value(e1[i][j]))
            energy_b = float(model.get_py_value(e2[i][j]))
            if energy_a != 0.0:
                print('E[{}][{}] = {}'.format(i,j,energy_a))
            if energy_b != 0.0:
                print('e[{}][{}] = {}'.format(i,j,energy_b))
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
                structure[i] = page1_chars[0]
                structure[j] = page1_chars[1]
                print('X[{}][{}] = {}'.format(i,j,bp_a))
            elif bp_b:
                if not chars_set:
                    page2_chars = ('(',')')
                    page1_chars = ('[',']')
                    chars_set = True
                structure[i] = page2_chars[0]
                structure[j] = page2_chars[1]
                print('Y[{}][{}] = {}'.format(i,j,bp_b))

    structure = ''.join(structure)
    return structure,total_energy

if __name__ == "__main__":

    pkb115_seq = 'CGGUAGCGCGAACCUUAUCGCGCA'
    pkb115_struct = '.(((.[[[[[[)))...]]]]]].'
    predict_structure(pkb115_seq,doubleknots=False)


