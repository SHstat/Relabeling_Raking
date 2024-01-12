import numpy as np
import pandas as pd
import random

from itertools import compress
from sklearn.neighbors import NearestNeighbors


def OREM(data, label, N, q=5):
    minclass = 1
    minSet = list(compress(list(data.index), [l==minclass for l in label]))

    cas, pairs = getCAS(data, label, q=q)
    AS = getAS(cas, pairs, data)
    syn = generate(AS, minSet, N, data).T
    
    newLabel = label + [minclass] * len(syn)
    balanced = pd.concat([data, syn])

    return balanced, newLabel


def get_pairs(array):
    res = {}
    if len(array):
        maj_loc = np.where(array == 0)[0]
        if len(maj_loc):
            new = range(maj_loc[0]+1, len(array))
            res = {n: [m for m in maj_loc if n > m] for n in new}
    return res


def getCAS(data, label, q):
    """
    Get CAS (the candidate assistant seeds)
    
    Parameters
    ----------
    label: Minority Class = 1, Majority Class = 0
    """
    li_df = pd.DataFrame([label, data.index], index = ['Label', 'Index']).T
    idxMin = li_df.loc[li_df.Label.values == 1].index
    label_dict = li_df.T.loc['Label'].to_dict()
    index_dict = li_df.T.loc['Index'].to_dict()

    nbrs = NearestNeighbors(n_neighbors=len(data), metric='euclidean').fit(data)
    distances, indices = nbrs.kneighbors(data)

    ind_df = pd.DataFrame(indices).loc[idxMin]
    ind_df_ori = ind_df.replace(index_dict)
    ind_df_lab = ind_df.replace(label_dict)
    
    successionalMaj = ind_df_lab.iloc[:, 1:].rolling(q, axis=1).sum() == 0
    cas_max = successionalMaj.idxmax(axis=1)-(q)
    
    CAS = {index_dict[i]: ind_df_ori.loc[i, 1:cas_max[i]].values for i in ind_df.index}
    CAS = {k: v for k, v in CAS.items() if v.any()}
    
    pairs = {index_dict[i]: get_pairs(ind_df_lab.loc[i, 1:cas_max[i]].values) for i in ind_df.index}
    
    return CAS, pairs


def getAS(CAS, pairs, data):
    res = {}
    for k, p in pairs.items():
        a = data.loc[k]
        
        if p:
            AS = list(CAS[k][:list(p.keys())[0]])
            for n, os in p.items():
                newIdx = CAS[k][n]
                new = data.loc[newIdx]
                m = (a + new)/2
                
                dist_m_a = np.linalg.norm(m-a)
                assistant = True
                for o in os:
                    oldIdx = CAS[k][o]
                    old = data.loc[oldIdx]

                    dist_m_o = np.linalg.norm(m-old)
                    if dist_m_a >= dist_m_o:
                        assistant = False
                        break
                
                if assistant:
                    AS.append(newIdx)
            res[k] = AS
    
    return res


def generate(AS, minSet, N, data):
    res = []
    N = int(N/len(AS))
    
    for i, assistants in AS.items():
        Xi = data.loc[i]
        for n in range(N):
            s = random.choices(assistants, k=1)[0]
            Xs = data.loc[s]

            gamma = random.random()
            if s not in minSet: gamma /= 2

            res.append(Xi + gamma * (Xs - Xi))
    
    return pd.concat(res, axis=1)