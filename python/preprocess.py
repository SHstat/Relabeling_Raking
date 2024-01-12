import pandas as pd

def check_uniqueness(df, change = True):
    if len(set([df[c].nunique() for c in df.columns])) == 1:
        print(f"All column variables have the same number of unique values")
        # this doesn't mean that those unique values are the same for every columns

        tmp = pd.DataFrame([sorted(df[c].unique()) for c in df.columns])
        if all([tmp[c].nunique()==1 for c in tmp]):
            print(f"Even they have the same unique values... {[tmp[c].unique()[0] for c in tmp]}")
            # We will change the values for the special case only, otherwise we will make it by dummies.
            if change & (len(tmp.columns)==2):
                tmp_dict = {tmp[c].unique()[0]: new for new, c in enumerate(tmp)}
                newDf = df.replace(tmp_dict)
                print(f"New rules are applied: {tmp_dict}")
                return newDf
        else:
            print(f"But that doesn't mean that all the unique values are common across all variables.")
    else:
        print(f"There is at least one column which has different length of unique values.")
    
    return df


def check_class(df):
    vc = df.iloc[:, -1].value_counts()    
    if len(vc) > 2:
        print(f"CAUTION: it has {len(vc)} classes")

    classMax = vc[vc.values == vc.max()].index[0]
    classMin = vc[vc.values == vc.min()].index[0]

    print(f"# of {classMax}: {vc[classMax]}")
    print(f"# of {classMin}: {vc[classMin]}")

    # We need to use index when we replace the values
    # This is because one of existing value could be 1 or 0.
    if classMax != 0:
        df.loc[df[df.iloc[:, -1].values == classMax].index, df.columns[-1]] = 0
        print(f"Majority Label changed: {classMax}->{0}")

    if classMin != 1:
        df.loc[df[df.iloc[:, -1].values == classMin].index, df.columns[-1]] = 1
        print(f"Minority Label changed: {classMin}->{1}")

    Label = df.iloc[:, -1].astype(int)
    df = df.iloc[:, :-1].copy()
    df['LABEL'] = Label

    return df    
