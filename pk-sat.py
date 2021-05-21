from pysmt.shortcuts import Symbol, And, Not, Plus, Equals, Real,GE, Implies, get_model, is_sat

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
        return None

def predict_structure(seq):

    conditions = []
    p1 =  []
    p2 = []
    e1 = []
    e2 = []

    # intialize n*n pairing variables for two pages
    for i in range(len(seq)):
        p1_row = []
        p2_row = []
        for j in range(len(seq)):
            x = Symbol('X_{},{}'.format(i,j)) 
            y = Symbol('Y_{},{}'.format(i,j)) 
            p1_row.append(x)
            p2_row.append(y)
        p1.append(p1_row)
        p2.append(p2_row)

    # initialize n*n energy scores for two pages
    for i in range(len(seq)):
        e1_row = []
        e2_row = []
        for j in range(len(seq)):
            e = Symbol('E_{},{}'.format(i,j)) 
            E = Symbol('e_{},{}'.format(i,j)) 
            e1_row.append(e)
            e2_row.append(E)
        e1.append(e1_row)
        e2.append(e2_row)


def build_energy_constraints(seq,conditions,p1,p2,e1,e2):

    # assign scores based on bp stacks
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            nuc1 = seq[i]
            nuc2 = seq[j]
            nuc3 = seq[i+1]
            nuc4 = seq[j-1]
            bp_stack_energy = energy_score(nuc1,nuc2,nuc3,nuc4) 
            e1_stack_score = Implies(And(p1[i][j],p1[i+1][j-1]),Equals(e1[i][j],Real(bp_stack_energy)))
            e2_stack_score = Implies(And(p2[i][j],p2[i+1][j-1]),Equals(e2[i][j],Real(bp_stack_energy)))
            e1_nostack_score = Implies(Not(And(p1[i][j],p1[i+1][j-1])),Equals(e1[i][j],Real(0.0)))
            e2_nostack_score = Implies(And(p2[i][j],p2[i+1][j-1]),Equals(e2[i][j],Real(0.0)))
   
    # check if total energy threshold is exceeded



def build_structural_constraints(seq,conditions,p1,p2):


    # every position can be in at most one bp
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            for k in range(j+1,len(seq)):
                p1_internal = Not(And(p1[i][j],p1[j][k]))
                p2_internal = Not(And(p2[i][j],p2[j][k]))
                p1_first = Not(And(p1[i][j],p2[j][k]))
                p2_first = Not(And(p2[i][j],p1[j][k]))
                new_conditions = [p1_internal,p2_internal,p1_first,p2_first]
                conditions.extend(new_conditions)

    # bp must belong to either page1 or page2
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            unique_page = Not(And(p1[i][j],p2[i][j]))
            conditions.append(unique_page)

    # no internal pseudoknots
    for i in range(len(seq)):
        for k in range(i+1,len(seq)):
            for j in range(k+1,len(seq)):
                for l in range(j+1,len(seq)):
                    p1_pk = Not(And(p1[i][j],p1[k][l]))
                    p2_pk = Not(And(p2[i][j],p2[k][l]))
                    conditions.append(p1_pk)
                    conditions.append(p2_pk)

    # only nested bp (no splits) on page2
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            for k in range(j+1,len(seq)):
                for l in range(k+1,len(seq)):
                    p2_split = Not(And(p2[i][j],p2[k][l]))
                    conditions.append(p2_split)

    # no double-crossing pseudoknots
    for i in range(len(seq)):
        for m in range(i+1,len(seq)):
            for j in range(m+1,len(seq)):
                for k in range(j+1,len(seq)):
                    for n in range(k+1,len(seq)):
                        for l in range(n+1,len(seq)):
                            downstream = Implies(And(p1[i][j],p2[m][n]),Not(p1[k][l]))
                            upstream = Implies(And(p1[k][l],p2[m][n]),Not(p1[i][j]))
                            conditions.append(downstream)
                            conditions.append(upstream)

    return conditions

def solve(conditions):

    # build structural SAT formula
    formula = conditions[0]
    print("# SAT clauses = {}".format(len(conditions)))
    for c in conditions[1:]:
       formula = And(formula,c)

    # find satisfying assignment if one exists
    model = get_model(formula)
    if model:
        print(model)
    else:
        print("No solution found")


if __name__ == "__main__":

    pkb115_seq = 'CGGUAGCGCGAACCUUAUCGCGCA'
    pkb115_struct = ':(((:[[[[[[))):::]]]]]]:'
    build_structural_constraints(pkb115_seq)


	
