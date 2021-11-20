import numpy as np
import pandas as pd
import scipy.stats as ss
import itertools as it
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.libqsturng import psturng
from pandas import DataFrame, Categorical, Series

def Cnk(a,b):
    iter1 = 1
    iter2 = 1
    iter3 = 1
    for i in range(1,a+1):
        iter1 = iter1*i
    for i in range(1,b+1):
        iter2 = iter2*i
    for i in range(1,a-b+1):
        iter3 = iter3*i
    return int(iter1/(iter2*iter3))
def __convert_to_df(a, val_col=None, group_col=None, val_id=None, group_id=None):
    if not group_col:
        group_col = 'groups'
    if not val_col:
        val_col = 'vals'

    if isinstance(a, DataFrame):
        x = a.copy()
        if not {group_col, val_col}.issubset(a.columns):
            raise ValueError('Specify correct column names using `group_col` and `val_col` args')
        return x, val_col, group_col

    elif isinstance(a, list) or (isinstance(a, np.ndarray) and not a.shape.count(2)):
        grps_len = map(len, a)
        grps = list(it.chain(*[[i+1] * l for i, l in enumerate(grps_len)]))
        vals = list(it.chain(*a))

        return DataFrame({val_col: vals, group_col: grps}), val_col, group_col

    elif isinstance(a, np.ndarray):

        # cols ids not defined
        # trying to infer
        if not(all([val_id, group_id])):

            if np.argmax(a.shape):
                a = a.T

            ax = [np.unique(a[:, 0]).size, np.unique(a[:, 1]).size]

            if np.diff(ax).item():
                __val_col = np.argmax(ax)
                __group_col = np.argmin(ax)
            else:
                raise ValueError('Cannot infer input format.\nPlease specify `val_id` and `group_id` args')

            cols = {__val_col: val_col,
                    __group_col: group_col}
        else:
            cols = {val_id: val_col,
                    group_id: group_col}

        cols_vals = dict(sorted(cols.items())).values()
        return DataFrame(a, columns=cols_vals), val_col, group_col
def posthoc_dunn(a, val_col=None, group_col=None, p_adjust=None, sort=True):
    def compare_dunn(i, j):
        diff = np.abs(x_ranks_avg.loc[i] - x_ranks_avg.loc[j])
        A = n * (n + 1.) / 12.
        B = (1. / x_lens.loc[i] + 1. / x_lens.loc[j])
        z_value = diff / np.sqrt((A - x_ties) * B)
        p_value = 2. * ss.norm.sf(np.abs(z_value))
        
        #print(diff,np.sqrt((A - x_ties) * B),z_value,p_value)
        return p_value,diff,np.sqrt((A - x_ties) * B),z_value,p_value

    x, _val_col, _group_col = __convert_to_df(a, val_col, group_col)

    if not sort:
        x[_group_col] = Categorical(x[_group_col], categories=x[_group_col].unique(), ordered=True)

    x.sort_values(by=[_group_col, _val_col], ascending=True, inplace=True)
    n = len(x.index)
    x_groups_unique = np.unique(x[_group_col])
    x_len = x_groups_unique.size
    x_lens = x.groupby(_group_col)[_val_col].count()

    x['ranks'] = x[_val_col].rank()
    x_ranks_avg = x.groupby(_group_col)['ranks'].mean()

    # ties
    vals = x.groupby('ranks').count()[_val_col].values
    tie_sum = np.sum(vals[vals != 1] ** 3 - vals[vals != 1])
    tie_sum = 0 if not tie_sum else tie_sum
    x_ties = tie_sum / (12. * (n - 1))

    vs = np.zeros((x_len, x_len))
    combs = it.combinations(range(x_len), 2)

    tri_upper = np.triu_indices(vs.shape[0], 1)
    tri_lower = np.tril_indices(vs.shape[0], -1)
    vs[:,:] = 0
    
    Diff = []
    SE = []
    Z = []
    P = []
    for i,j in combs:
        vs[i, j] = compare_dunn(x_groups_unique[i], x_groups_unique[j])[0]
        Diff.append(compare_dunn(x_groups_unique[i], x_groups_unique[j])[1])
        SE.append(compare_dunn(x_groups_unique[i], x_groups_unique[j])[2])
        Z.append(compare_dunn(x_groups_unique[i], x_groups_unique[j])[3])
        P.append(compare_dunn(x_groups_unique[i], x_groups_unique[j])[4])
    if p_adjust:
        vs[tri_upper] = multipletests(vs[tri_upper], method = p_adjust)[1]

    vs[tri_lower] = vs.T[tri_lower]
    np.fill_diagonal(vs, -1)
    return DataFrame(vs, index=x_groups_unique, columns=x_groups_unique), Diff,SE,Z,P