import pandas as pd
import sys
from utils import mkLaTexTable

def analyze(trials):

    df = pd.read_csv(trials)
    df['recall_improvement'] = df['recall'] - df['lp_recall']
    df['precision_improvement'] = df['precision'] - df['lp_precision']
    df['f1_improvement'] = df['f1'] - df['lp_f1']
    df['overall_improvement'] = df['recall_improvement'] + df['precision_improvement']    
    df['total_pairs'] = df['pairs1']+df['pairs2'] 
    top_trial_idx = (df['k1'] == 4) & (df['k2'] == 3)
    #top_trial_idx = (df['k1'] == 2) & (df['k2'] == 2)
    top_trial = df[top_trial_idx]
    print(table1_to_latex(top_trial))
    by_k = df.groupby(['k1','k2']).mean()
    
    by_k.to_csv('improvement_by_k.csv',index=False)
    top_trial.to_csv('top_trial.csv',index=False)
    print(table2_to_latex(by_k))
def table1_to_latex(df):
    return df.to_latex(index=False, columns=['name', 'len','total_pairs', 'time', 'iterations','precision','recall','precision_improvement','recall_improvement'])

def table2_to_latex(df):
    return df.to_latex(index=True,columns=['time','iterations','precision','recall','precision_improvement','recall_improvement','overall_improvement'])

if __name__ == "__main__":

    analyze(sys.argv[1])
