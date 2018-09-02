
## Two-sided
## Version: 20160912

import pandas as pd
import h5py
import numpy as np
import os

NUM_TS = 52
HOMO_PENALTY = 0.33
SIM_THRESHOLD = 0.95
MIN_SEQLEN = 25
LINKER_SEQ = 'CGGATC'

CNT_TAG = 'TOTAL_COUNT'
LSEQ_TAG = 'L_QSEQ'
LLEN_TAG = 'L_QLEN'
LTYP_TAG = 'L_TYPE'
RSEQ_TAG = 'R_QSEQ'
RLEN_TAG = 'R_QLEN'
RTYP_TAG = 'R_TYPE'

tot_prob_cnt = 0

## Clean data
def clean_seq(seq):
    seq = seq.replace('N','')
    loc = seq.find(LINKER_SEQ)
    if loc != -1 and loc >= len(seq) - 40:
        seq = seq[:loc]
    return seq

def calc_dist(seq1, seq2, SIM_THRESHOLD):
    seq1 = seq1[:len(seq2)]
    seq2 = seq2[:len(seq1)]

    len1 = len(seq1)
    len2 = len(seq2)
    last_col = []

    v0 = [0]*(len2+1)  # the previous row
    v1 = [0]*(len2+1)  # the current row
    for k in range(len2+1):
        v0[k] = k

    curr_dist = 0
    for i in range(len1): # i is the true ind for seq1
        v1[0] = i+1
        for j in range(len2): # j is the true ind for seq2
            if i > 0 and seq1[i]==seq1[i-1] and seq1[i]==seq2[j]: # homopolymer ins in seq1
                delMin = v0[j+1] + HOMO_PENALTY
            else:
                delMin = v0[j+1] + 1 # insertion err in seq1
            if j > 0 and seq2[j]==seq2[j-1] and seq1[i]==seq2[j]:
                insMin = v1[j] + HOMO_PENALTY
            else:
                insMin = v1[j] + 1 # insertion err in seq2

            if seq1[i] == seq2[j]: # or seq1[i]=='N' or seq2[j]=='N':
                sameMin = v0[j]
                v1[j+1] = min([sameMin,delMin,insMin])
            else:
                subMin = v0[j] + 1
                v1[j+1] = min([subMin,delMin,insMin])
        last_col.append(v0[-1])
        v0 = list(v1)

        curr_dist = min(min(v0), min(last_col))/float(min(len(seq1), len(seq2)))
        if curr_dist > 1-SIM_THRESHOLD:
            return 1 # Not similar
    return curr_dist # similar

def handle_onesided(filename, side):
    if side == 'R':
        seq_lab = RSEQ_TAG
        seqlen_lab = RLEN_TAG
        match_lab = RTYP_TAG
    else:
        seq_lab = LSEQ_TAG
        seqlen_lab = LLEN_TAG
        match_lab = LTYP_TAG

    ## Read input
    input_data = pd.read_csv(filename, sep='\t')
    sLength = input_data.shape[0]
    input_data['Indl'] = pd.Series(np.random.randn(sLength), index=input_data.index)

    ## Clean data
    input_data.ix[:,seq_lab] = [clean_seq(seq) for seq in input_data[seq_lab]]
    input_data_invalid = input_data[input_data[seqlen_lab] < MIN_SEQLEN]
    input_data_invalid.to_csv('intermediate_input_shorterseqs_'+side+'.txt', index=False, sep='\t')
    input_data = input_data[input_data[seqlen_lab] >= MIN_SEQLEN]
    input_data.to_csv('intermediate_input_valid_'+side+'.txt', index=False, sep='\t')

    USE_UNSUPERVISE = True
    if USE_UNSUPERVISE:
        # So global ind is the real index in the pandas dataframe
        indgs_masts = input_data[(input_data[match_lab].str.contains('Multi'))
                | (input_data[match_lab] == 'Single')].index.tolist()
        num_masts = len(indgs_masts)

        ## Mapping Original inds to Special inds
        indls_masts = range(num_masts)

        input_data.ix[indgs_masts, 'Indl'] = indls_masts

        dict_masts = dict(zip(indls_masts, indgs_masts))

        ##
        # 20160825
        if os.path.isfile('intermediate_dist_matrix_unsup_'+side+'.hdf5'):
            f = h5py.File('intermediate_dist_matrix_unsup_'+side+'.hdf5', 'r')
            dist_matrix_unsup = f['dist_matrix_unsup'][:]
            dist_matrix_unsup.resize(num_masts, num_masts)
            f.close()
        else:
            f = h5py.File('intermediate_dist_matrix_unsup_'+side+'.hdf5', 'w')
            dist_matrix_unsup = np.asarray(
            [calc_dist(seq1, seq2, SIM_THRESHOLD)
            for seq1 in input_data.ix[indgs_masts][seq_lab]
                for seq2 in input_data.ix[indgs_masts][seq_lab]]
            )
            dist_matrix_unsup = dist_matrix_unsup.reshape((num_masts, num_masts))
            f.create_dataset('dist_matrix_unsup', data=dist_matrix_unsup)
            f.close()

        # ranking
        def rank_seqs(seq_rows):
            scores = [0]*len(seq_rows)
            dict_hitTypes = {'Multi':0.67, 'Single':0.33}
            for i in range(len(seq_rows)):
                scores[i] += ('Multi' in seq_rows[i].iloc[0][match_lab]) * 100
                scores[i] += seq_rows[i].iloc[0][CNT_TAG] * 10
                scores[i] += seq_rows[i].iloc[0][seqlen_lab] * 1
            return scores.index(max(scores))

        # grouping, not parallelable
        visited = np.zeros((num_masts,num_masts)) #
        queue = []
        groups = [] # will contain global ind
        for i in indls_masts:
            for j in indls_masts: # j > i!!!
                if j <= i:
                    continue
                if not visited[i][j] and dist_matrix_unsup[i][j] < 1 - SIM_THRESHOLD:
                    group = [] # the ind in the group is local
                    queue.insert(0, [i,j])
                    visited[i][j] = True
                    while len(queue) > 0:
                        curr_node = queue.pop()
                        group.append(curr_node[0])
                        group.append(curr_node[1])
                        # first, search through (i1,j) where i1 < i
                        for i1 in range(curr_node[1]):
                            if (not visited[i1][curr_node[1]]) and \
                                            dist_matrix_unsup[i1][curr_node[1]] < 1 - SIM_THRESHOLD:
                                queue.insert(0, [i1, curr_node[1]])
                                visited[i1][curr_node[1]] = True
                        for j2 in range(curr_node[1]+1, num_masts):
                            if (not visited[curr_node[1]][j2]) and\
                                            dist_matrix_unsup[curr_node[1]][j2] < 1 - SIM_THRESHOLD:
                                queue.insert(0, [curr_node[1], j2])
                                visited[curr_node[1]][j2] = True

                        for j1 in range(curr_node[0]+1, num_masts):
                            if (not visited[curr_node[0]][j1]) and \
                                            dist_matrix_unsup[curr_node[0]][j1] < 1 - SIM_THRESHOLD:
                                queue.insert(0, [curr_node[0], j1])
                                visited[curr_node[0]][j1] = True
                        for i2 in range(curr_node[0]):
                            if (not visited[i2][curr_node[0]]) and\
                                            dist_matrix_unsup[i2][curr_node[0]] < 1 - SIM_THRESHOLD:
                                queue.insert(0, [i2, curr_node[0]])
                                visited[i2][curr_node[0]] = True
                    group = list(set(group))
                    # Above: finished constructing a group of local inds

                    # Next: Ranking
                    curr_best = rank_seqs([input_data.loc[[dict_masts[k]]] for k in group])

                    # Put the best one at the loc 0
                    best_mast_indl = group.pop(curr_best)
                    group.insert(0, best_mast_indl)

                    # Method 1: Add all cnts to the best one
                    best_mast_indg = dict_masts[best_mast_indl]
                    for k in range(1, len(group)):  # local ind
                        curr_mast_indg = dict_masts[group[k]]
                        input_data.ix[best_mast_indg, CNT_TAG] +=\
                            input_data.ix[curr_mast_indg, CNT_TAG]
                        input_data.ix[best_mast_indg, col_names_ts] +=\
                            input_data.ix[curr_mast_indg, col_names_ts]
                        input_data.ix[curr_mast_indg, 'Indl'] = -1

                    groups += [dict_masts[one_local_ind] for one_local_ind in group]
                    groups += ['']

        #input_data.loc[groups].to_csv('intermediate_output_454_masts_grouped_'+side+'.txt', index=False, sep='\t')
        input_data = input_data[(input_data[match_lab].str.contains('NoGoodSpn')
                | input_data[match_lab].str.contains('NoHits')) | (input_data['Indl'] >= 0)] # Get rid of merged stuff
        #input_data.drop('Indl',axis=1,inplace=True)
        input_data.to_csv('intermediate_output_454_masts_merged_'+side+'.txt', index=False,  sep='\t')
        ##
        dict_masts.clear()


    ###### Supervised now ######

    ## Original indices
    indgs_masts = input_data[(input_data[match_lab].str.contains('Multi'))
            | (input_data[match_lab] == 'Single')].index.tolist()
    indgs_slavs = input_data[(input_data[match_lab].str.contains('NoGoodSpn'))
            | (input_data[match_lab].str.contains('NoHits'))].index.tolist()

    num_masts = len(indgs_masts)
    num_slavs = len(indgs_slavs)

    ## Mapping local inds to global inds

    indls_masts = range(num_masts)
    dict_masts = dict(zip(indls_masts, indgs_masts))
    indls_slavs = range(num_slavs)
    dict_slavs = dict(zip(indls_slavs, indgs_slavs))

    ## Mapping global inds to local inds
    input_data.ix[indgs_masts,'Indl'] = indls_masts
    input_data.ix[indgs_slavs,'Indl'] = indls_slavs

    ## Distance matrix
    if os.path.isfile('intermediate_dist_matrix_super_'+side+'.hdf5'):
        f = h5py.File('intermediate_dist_matrix_super_'+side+'.hdf5', 'r')
        dist_matrix_super = f['dist_matrix_super'][:]
        dist_matrix_super.resize(num_masts, num_slavs)
        f.close()
    else:
        f = h5py.File('intermediate_dist_matrix_super_'+side+'.hdf5', 'w')
        dist_matrix_super = np.asarray(
        [calc_dist(seq1, seq2, SIM_THRESHOLD)
        for seq1 in input_data.ix[indgs_masts,seq_lab]
            for seq2 in input_data.ix[indgs_slavs,seq_lab]]
                  )
        dist_matrix_super = dist_matrix_super.reshape((num_masts, num_slavs))
        f.create_dataset('dist_matrix_super', data=dist_matrix_super)
        f.close()

    # The best way to remember this is that the order of for loop inside the list comprehension
    # is based on the order in which they appear in traditional loop approach. Outer most loop comes
    # first, and then the inner loops subsequently.

    ## Merge slavs to masts
    for i in indls_masts:
        Indls_int = [int(one_Indl) for one_Indl in input_data['Indl']]
        slavs_indl = [one_ind for one_ind in
            input_data[input_data.index.isin(indgs_slavs) &
                (dist_matrix_super[i][Indls_int] < 1 - SIM_THRESHOLD)]['Indl']]

        if len(slavs_indl) > 0:
            for j in slavs_indl:
                input_data.ix[dict_masts[i], CNT_TAG] += input_data.ix[dict_slavs[j], CNT_TAG]
                input_data.ix[dict_masts[i], col_names_ts] += input_data.ix[dict_slavs[j], col_names_ts]

    err_slavs = []
    err_slav_masts = []
    ## Get rid of error slavs
    for j in indls_slavs:
        masts_indgs = [dict_masts[i] for i in indls_masts
                       if dist_matrix_super[i][j] < 1 - SIM_THRESHOLD]
        slav_indg = dict_slavs[j]

        if len(masts_indgs) > 1:
            err_slavs.append(slav_indg)
            err_slav_masts.append(masts_indgs)

            for mast_indg in masts_indgs:
                input_data.ix[mast_indg, CNT_TAG] -= input_data.ix[slav_indg, CNT_TAG]
                input_data.ix[mast_indg, col_names_ts] -= input_data.ix[slav_indg, col_names_ts]
            input_data.ix[slav_indg, 'Indl'] = -2
        elif len(masts_indgs) == 1:
            input_data.ix[slav_indg, 'Indl'] = -3 # -input_data.ix[slav_indg, 'Indl']

    ## Write to file
    curr_merged = input_data[input_data.index.isin(indgs_masts)
        | (input_data.index.isin(indgs_slavs) & (input_data['Indl'] >= 0))]
    curr_merged.to_csv('intermediate_output_454_merged_'+side+'.txt', index=False, sep='\t')

    err_slavs_df = pd.DataFrame(columns=input_data.columns.tolist())
    for i in range(len(err_slavs)):
        err_slavs_df = err_slavs_df.append(pd.DataFrame(input_data.loc[[err_slavs[i]]]))
        for one_mast in err_slav_masts[i]:
            err_slavs_df = err_slavs_df.append(pd.DataFrame(input_data.loc[[one_mast]]))

        err_slavs_df = err_slavs_df.append({'CLONE_ID':'-'}, ignore_index=True)
        #print input_data.ix[err_slavs[i], CNT_TAG]
    err_slavs_df.to_csv('intermediate_output_454_problemones_'+side+'.txt', index=False, sep='\t')


DISTTHRESHOLD = 20
def linkable_RL(seq1, seq2):
    dict_turn = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    LINK_THRESHOLD = 0.6
    seq1 = ''.join([dict_turn[n] for n in seq1[:5]]) # code
    seq2 = seq2[:5][::-1] # take the first 5, and reverse
    return calc_dist(seq1, seq2, LINK_THRESHOLD) < 1 - LINK_THRESHOLD

SITE_TAG = 'CHR_SITE1'
ORIE_TAG = 'R_STRAND'
def link_RLseqs(input_data):
    input_data = input_data.sort_values([SITE_TAG], ascending=False)

    for i in range(input_data.shape[0] - 1):
        curr_Rloc = input_data.iloc[i][SITE_TAG]
        next_Rloc = input_data.iloc[i+1][SITE_TAG]

        curr_ORI = input_data.iloc[i][ORIE_TAG]
        next_ORI = input_data.iloc[i+1][ORIE_TAG]

        curr_Rseq = input_data.iloc[i][RSEQ_TAG]
        next_Rseq = input_data.iloc[i+1][RSEQ_TAG]

        curr_Lseq = input_data.iloc[i][LSEQ_TAG]
        next_Lseq = input_data.iloc[i+1][LSEQ_TAG]

        if next_Rloc == -1:
            break
        if (curr_Rloc - next_Rloc > DISTTHRESHOLD)\
                or (curr_ORI != next_ORI):
            continue

        # NOTE: isinstance(curr_Rseq, float) means MISSING data
        if curr_Lseq > 0 and isinstance(curr_Rseq, float)\
                and next_Rseq > 0 and isinstance(next_Lseq, float):
            if linkable_RL(curr_Lseq, next_Rseq):
                curr_index = input_data.index[i]
                input_data.loc[curr_index,RSEQ_TAG] = next_Rseq
                input_data.loc[curr_index,RLEN_TAG] = len(next_Rseq)
                input_data.loc[curr_index,CNT_TAG] += input_data.iloc[i+1][CNT_TAG]
                input_data.loc[curr_index,col_names_ts] += input_data.iloc[i+1][col_names_ts]

                input_data.drop(input_data.index[i+1], inplace=True)

        if curr_Rseq > 0 and isinstance(curr_Lseq, float)\
                and next_Lseq > 0 and isinstance(next_Rseq, float):
            if linkable_RL(curr_Rseq, next_Lseq):
                curr_index = input_data.index[i]
                input_data.loc[curr_index,LSEQ_TAG] = next_Lseq
                input_data.loc[curr_index,LLEN_TAG] = len(next_Lseq)
                input_data.loc[curr_index,CNT_TAG] += input_data.iloc[i+1][CNT_TAG]
                input_data.loc[curr_index,col_names_ts] += input_data.iloc[i+1][col_names_ts]

                input_data.drop(input_data.index[i+1], inplace=True)
    return input_data


input_data = pd.read_csv('pipe15_output.txt', sep='\t')
#input_data = input_data.rename(columns={'Unnamed: 0':'Indl'})

all_colnames = input_data.columns.tolist()
#col_names_ts = []
#for one_colname in all_colnames:
#    if ('CD' in one_colname) or ('PBC' in one_colname):
#        col_names_ts.append(one_colname)
col_names_ts = all_colnames[19:]

#col_names_ts = ['LEFT TOTAL','RIGHT TOTAL',' LRQ55',' LRQ56',' LRQ57',' LRQ58',' RRQ55',' RRQ56',' RRQ57',' RRQ58']

input_data.ix[input_data[RLEN_TAG]>0,RSEQ_TAG] =\
    [clean_seq(seq) for seq in input_data[input_data[RLEN_TAG]>0][RSEQ_TAG]]
input_data.ix[input_data[RLEN_TAG]>0,RLEN_TAG] =\
    [len(seq) for seq in input_data[input_data[RLEN_TAG]>0][RSEQ_TAG]]
input_data.ix[input_data[LLEN_TAG]>0,LSEQ_TAG] =\
    [clean_seq(seq) for seq in input_data[input_data[LLEN_TAG]>0][LSEQ_TAG]]
input_data.ix[input_data[LLEN_TAG]>0,LLEN_TAG] =\
    [len(seq) for seq in input_data[input_data[LLEN_TAG]>0][LSEQ_TAG]]

# TODO: what if both sides < 25?
input_data_invalid = input_data[(input_data[RLEN_TAG]+input_data[LLEN_TAG]) < MIN_SEQLEN]
input_data_invalid.to_csv('intermediate_input_shorterseqs.txt', index=False, sep='\t')
input_data = input_data[(input_data[RLEN_TAG]+input_data[LLEN_TAG]) >= MIN_SEQLEN]
input_data.to_csv('intermediate_input_valid.txt', index=False, sep='\t')
# print 'total valid cnt: ' + str(sum(input_data[CNT_TAG]))

sLength = input_data.shape[0]
input_data['Indl'] = pd.Series(np.random.randn(sLength), index=input_data.index)

DO_GENERATE_FILES = True
if DO_GENERATE_FILES:
    indgs_Tsided = input_data[(input_data[RLEN_TAG] > 0)
            & (input_data[LLEN_TAG] > 0)].index.tolist()
    num_Tsided = len(indgs_Tsided)
    indls_Tsided = range(num_Tsided)
    input_data.ix[indgs_Tsided,'Indl'] = indls_Tsided
    dict_Tsided = dict(zip(indls_Tsided, indgs_Tsided))

    indgs_Rsided = input_data[(input_data[RLEN_TAG] > 0)
            & (input_data[LLEN_TAG] == 0)].index.tolist()
    num_Rsided = len(indgs_Rsided)
    indls_Rsided = range(num_Rsided)
    input_data.ix[indgs_Rsided,'Indl'] = indls_Rsided
    dict_Rsided = dict(zip(indls_Rsided, indgs_Rsided))

    ################
    ## Merge R- to T-
    if os.path.isfile('intermediate_dist_matrix_TRsided.hdf5'):
        f = h5py.File('intermediate_dist_matrix_TRsided.hdf5', 'r')
        dist_matrix_TRsided = f['dist_matrix_TRsided'][:]
        dist_matrix_TRsided.resize(num_Tsided, num_Rsided)
        f.close()
    else:
        f = h5py.File('intermediate_dist_matrix_TRsided.hdf5', 'w')
        dist_matrix_TRsided = np.asarray(
        [calc_dist(seq1, seq2, SIM_THRESHOLD)
        for seq1 in input_data.ix[indgs_Tsided,RSEQ_TAG]
            for seq2 in input_data.ix[indgs_Rsided,RSEQ_TAG]]
        )
        dist_matrix_TRsided = dist_matrix_TRsided.reshape((num_Tsided, num_Rsided))
        f.create_dataset('dist_matrix_TRsided', data=dist_matrix_TRsided)
        f.close()

    ## Merge slavs to masts
    input_data_Rsided = input_data.loc[indgs_Rsided]
    for i in indls_Tsided:
        Indls_R_int = [int(one_Indl) for one_Indl in input_data_Rsided['Indl']]
        slavs_indl = [one_ind for one_ind in
            input_data_Rsided[dist_matrix_TRsided[i][Indls_R_int] < 1 - SIM_THRESHOLD]['Indl']]

        if len(slavs_indl) > 0:
            for j in slavs_indl:
                input_data.ix[dict_Tsided[i], CNT_TAG] += input_data.ix[dict_Rsided[j], CNT_TAG]
                input_data.ix[dict_Tsided[i], col_names_ts] += input_data.ix[dict_Rsided[j], col_names_ts]

    err_Rsided = []
    err_Rsided_Ts = []
    ## Get rid of error slavs
    for j in indls_Rsided:
        masts_indgs = [dict_Tsided[i] for i in indls_Tsided
                       if dist_matrix_TRsided[i][j] < 1 - SIM_THRESHOLD]
        slav_indg = dict_Rsided[j]

        if len(masts_indgs) > 1:
            err_Rsided.append(slav_indg)
            err_Rsided_Ts.append(masts_indgs)

            for mast_indg in masts_indgs:
                input_data.ix[mast_indg, CNT_TAG] -= input_data.ix[slav_indg, CNT_TAG]
                input_data.ix[mast_indg, col_names_ts] -= input_data.ix[slav_indg, col_names_ts]
            input_data.ix[slav_indg, 'Indl'] = -2

        elif len(masts_indgs) == 1:
            input_data.ix[slav_indg, 'Indl'] = -3 #

    ################
    ## Merge L- to T-
    indgs_Lsided = input_data[(input_data[RLEN_TAG] == 0)
            & (input_data[LLEN_TAG] > 0)].index.tolist()
    num_Lsided = len(indgs_Lsided)
    indls_Lsided = range(num_Lsided)
    input_data.ix[indgs_Lsided, 'Indl'] = indls_Lsided
    dict_Lsided = dict(zip(indls_Lsided, indgs_Lsided))

    input_data_Lsided = input_data.loc[indgs_Lsided]

    if os.path.isfile('intermediate_dist_matrix_TLsided.hdf5'):
        f = h5py.File('intermediate_dist_matrix_TLsided.hdf5', 'r')
        dist_matrix_TLsided = f['dist_matrix_TLsided'][:]
        dist_matrix_TLsided.resize(num_Tsided, num_Lsided)
        f.close()
    else:
        f = h5py.File('intermediate_dist_matrix_TLsided.hdf5', 'w')
        dist_matrix_TLsided = np.asarray(
        [calc_dist(seq1, seq2, SIM_THRESHOLD)
        for seq1 in input_data.ix[indgs_Tsided,LSEQ_TAG]
            for seq2 in input_data.ix[indgs_Lsided,LSEQ_TAG]]
        )
        dist_matrix_TLsided = dist_matrix_TLsided.reshape((num_Tsided, num_Lsided))
        f.create_dataset('dist_matrix_TLsided', data=dist_matrix_TLsided)
        f.close()

    ## Merge slavs to masts
    for i in indls_Tsided:
        Indls_L_int = [int(one_Indl) for one_Indl in input_data_Lsided['Indl']]
        slavs_indl = [one_ind for one_ind in
            input_data_Lsided[dist_matrix_TLsided[i][Indls_L_int] < 1 - SIM_THRESHOLD]['Indl']]

        if len(slavs_indl) > 0:
            for j in slavs_indl:
                input_data.ix[dict_Tsided[i], CNT_TAG] += input_data.ix[dict_Lsided[j], CNT_TAG]
                input_data.ix[dict_Tsided[i], col_names_ts] += input_data.ix[dict_Lsided[j], col_names_ts]

    err_Lsided = []
    err_Lsided_Ts = []
    ## Get rid of error slavs
    for j in indls_Lsided:
        masts_indgs = [dict_Tsided[i] for i in indls_Tsided
                       if dist_matrix_TLsided[i][j] < 1 - SIM_THRESHOLD]
        slav_indg = dict_Lsided[j]

        if len(masts_indgs) > 1:
            err_Lsided.append(slav_indg)
            err_Lsided_Ts.append(masts_indgs)

            for mast_indg in masts_indgs:

                input_data.ix[mast_indg, CNT_TAG] -= input_data.ix[slav_indg, CNT_TAG]
                input_data.ix[mast_indg, col_names_ts] -= input_data.ix[slav_indg, col_names_ts]
            input_data.ix[slav_indg, 'Indl'] = -2
        elif len(masts_indgs) == 1:
            input_data.ix[slav_indg, 'Indl'] = -3 # -input_data.ix[slav_indg, 'Indl']

    ################

    ################

    DO_LINK_RLSEQS = True
    if DO_LINK_RLSEQS:
        input_data = link_RLseqs(input_data)
    else:
        input_data = input_data.sort_values([SITE_TAG], ascending=False)
    ################

    ## Write to file, T
    pd_Tsided = input_data[(input_data[RLEN_TAG] > 0) & (input_data[LLEN_TAG] > 0)]
    pd_Tsided.to_csv('intermediate_merged_Tsided.txt', index=False, sep='\t')

    # Handle T-side
    err_slavs_df = pd.DataFrame(columns=input_data.columns.tolist())
    err_slavs = err_Rsided + err_Lsided
    err_slav_masts = err_Rsided_Ts + err_Lsided_Ts
    for i in range(len(err_slavs)):
        err_slavs_df = err_slavs_df.append(pd.DataFrame(input_data.loc[[err_slavs[i]]]))
        for one_mast in err_slav_masts[i]:
            err_slavs_df = err_slavs_df.append(pd.DataFrame(input_data.loc[[one_mast]]))
        err_slavs_df = err_slavs_df.append({'CLONE_ID':'separator'}, ignore_index=True)

        #print input_data.ix[err_slavs[i], CNT_TAG]
    err_slavs_df.to_csv('intermediate_problemones_Tsided.txt', index=False, sep='\t')


    ## Write to file, R
    pd_Rsided = input_data[(input_data[RLEN_TAG] > 0) & (input_data[LLEN_TAG] == 0)\
                            & (input_data['Indl'] >= 0)]
    pd_Rsided.to_csv('intermediate_valid_Rsided.txt', index=False, sep='\t')

    ## Call function for R-side
    if not os.path.isfile('intermediate_output_454_merged_R.txt'):
        handle_onesided('intermediate_valid_Rsided.txt', 'R')

    ## Write to file, L
    pd_Lsided = input_data[(input_data[RLEN_TAG] == 0) & (input_data[LLEN_TAG] > 0)\
                            & (input_data['Indl'] >= 0)]
    pd_Lsided.to_csv('intermediate_valid_Lsided.txt', index=False, sep='\t')

    ## Call function for L-side
    if not os.path.isfile('intermediate_output_454_merged_L.txt'):
        handle_onesided('intermediate_valid_Lsided.txt', 'L')

## Output
problemones = pd.DataFrame(columns=input_data.columns.tolist())
problemones = problemones.append(pd.read_csv('intermediate_problemones_Tsided.txt', sep='\t'))
problemones = problemones.append(pd.read_csv('intermediate_output_454_problemones_R.txt',  sep='\t'))
problemones = problemones.append(pd.read_csv('intermediate_output_454_problemones_L.txt',  sep='\t'))
problemones.drop('Indl',axis=1,inplace=True)
#problemones.to_csv('problemones.txt', index=False, sep='\t')

##TODO: grouped Tmasts? and naming

merged_all = pd.DataFrame(columns=input_data.columns.tolist())
merged_all = merged_all.append(pd.read_csv('intermediate_merged_Tsided.txt', sep='\t'))
merged_all = merged_all.append(pd.read_csv('intermediate_output_454_merged_R.txt', sep='\t'))
merged_all = merged_all.append(pd.read_csv('intermediate_output_454_merged_L.txt', sep='\t'))
merged_all.drop('Indl',axis=1,inplace=True)
merged_all.to_csv('merged_all.txt', index=False, sep='\t')
#print 'total merged cnt: ' + str(sum(merged_all[CNT_TAG]))

# clean up the intermediate files
DO_DELETE_INTERMEDIATE = True
if DO_DELETE_INTERMEDIATE:
    all_files = os.listdir(os.getcwd())
    for one_file in all_files:
        if 'intermediate_' in one_file:
            os.remove(one_file)

