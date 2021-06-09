import pandas as pd
import sys

def analyze(trials):

    df = pd.read_csv(trials)
    df['recall_improvement'] = df['recall'] - df['lp_recall']
    df['precision_improvement'] = df['precision'] - df['lp_precision']
    print(df)
    print(df.groupby(['k1','k2']).mean())

if __name__ == "__main__":

    analyze(sys.argv[1])
