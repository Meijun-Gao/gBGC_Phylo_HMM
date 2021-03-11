from scipy import stats
import numpy as np
import re


def standard_error(vector):
    num = len(vector)
    vector_std = np.std(vector)
    vector_se = vector_std/float(np.sqrt(num))
    return vector_se


def get_LRT_pvalue(pathnew, pathold):
    # pathnew = 'data/simulation/taxa4/h0/new/new-1result.txt'
    # pathold = 'data/simulation/taxa4/h0/simple/old-1result.txt'
    freedom = 8

    lines_old = open(pathold, 'r').readlines()
    lines_new = open(pathnew, 'r').readlines()
    old_LV_begin = float(lines_old[1][lines_old[1].index(':')+1:-1])
    old_LV_end = float(lines_old[2][lines_old[2].index(':') + 1:-1])
    new_LV_begin = float(lines_new[1][lines_new[1].index(':') + 1:-1])
    new_LV_end = float(lines_new[2][lines_new[2].index(':') + 1:-1])
    ch2_value = (old_LV_end - new_LV_end) * (-2)
    pvalue = stats.chi2.sf(ch2_value, freedom)

    return pvalue


def get_hotprob(numTaxa, hotspot_num, loadpath):
    if numTaxa == 4:
        K = 3  # states number in old model
        numSites = 5000
    elif numTaxa == 5:
        K = 15
        numSites = 2000

    hot_pro_mat = np.zeros((numSites, 1))
    back_pro_mat = np.zeros((numSites, 1))
    lines = open(loadpath, 'r').readlines()
    # find the index of posterior probability
    for x in range(200):
        if 'Position' in lines[x]:
            begini = x+1
            break

    for idx in range(numSites):
        temp_pro = list(map(float, re.split('[,]', lines[begini + idx][:-2])))

        pro = np.array(temp_pro)
        pro = np.delete(pro, 0)
        back_pro = np.sum(pro[0:K])
        hot_pro = np.sum(pro[K:2 * K])
        back_pro_mat[idx, 0] = back_pro
        hot_pro_mat[idx, 0] = hot_pro

    return hot_pro_mat, back_pro_mat


numTaxa = 4
hotspot_num = 0
newpath = 'data/simulation/taxa' + str(numTaxa) + '/h0/new/'
simplepath = 'data/simulation/taxa' + str(numTaxa) + '/h0/simple/'

if numTaxa == 4:
    numSites = 5000
    true_flag = np.zeros((1, numSites)).astype(int)
    if hotspot_num == 1:
        true_flag[0, 2000:4000] = 1
    elif hotspot_num == 2:
        true_flag[0, 1000:3000] = 1
        true_flag[0, 4000:4500] = 1
elif numTaxa == 5:
    numSites = 2000
    true_flag = np.zeros((1, numSites)).astype(int)
    if hotspot_num == 1:
        true_flag[0, 500:1300] = 1
    elif hotspot_num == 2:
        true_flag[0, 500:1200] = 1
        true_flag[0, 1500:1800] = 1
Num_trial = 20
hot_pro_all = np.zeros((Num_trial, numSites))
pvalue_list = []
""" average on 20 replicates """
for i in range(Num_trial):
    # if i == 5:
    #     continue
    loadpath_new = newpath + 'new-' + str(i+1) + 'result.txt'
    loadpath_simple = simplepath + 'old-' + str(i+1) + 'result.txt'

    hot_pro_mat, back_pro_mat = get_hotprob(numTaxa, hotspot_num, loadpath_new)
    hot_pro_all[i, :] = hot_pro_mat.T

    pvalue = get_LRT_pvalue(loadpath_new, loadpath_simple)
    pvalue_list.append(pvalue)

print('average of p_value: ' + str(np.mean(pvalue_list)))
print('SE of p_value: ' + str(standard_error(pvalue_list)))

# calculate average classification accuracy
thred1 = 0.95
thred2 = 0.75
acc = []
for t in range(Num_trial):
    hot_flag = []
    for i in range(numSites):
        # hotspot classification
        if hot_pro_all[t][i] >= thred1:
            # judge if this site locates in true hotspot
            hot_flag.append(1)
        else:
            hot_flag.append(0)
    hot_flag = np.array(hot_flag).reshape(1, numSites)
    acc_t = 1 - (np.sum(hot_flag ^ true_flag))/numSites
    acc.append(acc_t)

mean_acc = np.mean(acc)
se_acc = standard_error(acc)
print('average of acc: ' + str(mean_acc))
print('SE of acc: ' + str(se_acc))





