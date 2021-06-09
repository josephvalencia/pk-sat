import pandas as pd
import sys

def analyze(trials):

    df = pd.read_csv(trials)
    df['recall_improvement'] = df['recall'] - df['lp_recall']
    df['precision_improvement'] = df['precision'] - df['lp_precision']
    df['f1_improvement'] = df['f1'] - df['lp_f1']
    df['overall_improvement'] = df['recall_improvement'] + df['precision_improvement']    
    
    top_trial_idx = (df['k1'] == 2) & (df['k2'] == 1)
    top_trial = df[top_trial_idx]
    print(top_trial[['k1','k2','time','recall_improvement','precision_improvement']])
    by_k = df.groupby(['k1','k2']).mean()
    print(by_k)
    by_k.to_csv('improvement_by_k.csv',index=False)
    top_trial.to_csv('top_trial.csv',index=False)
if __name__ == "__main__":

    analyze(sys.argv[1])
