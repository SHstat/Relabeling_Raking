import numpy as np
import pandas as pd

def loading_train_test(type, fn, trial = 1, loc = "../data/"):
    t = str(trial)
    train = pd.read_csv(loc+type+"/"+fn+"_train"+t+".csv")
    test = pd.read_csv(loc+type+"/"+fn+"_test"+t+".csv")

    #print(train.shape, test.shape)
    #print("\nBalanced Train Dataset")
    #print(train.iloc[:, -1].value_counts())

    #print("\nImbalanced Test Dataset")
    #print(test.iloc[:, -1].value_counts())
    
    return train, test