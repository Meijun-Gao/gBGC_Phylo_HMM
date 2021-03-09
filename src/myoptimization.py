import numpy as np
from src.training import *
from src.HMM_algos import Forward
from src.my_tree import *


def objfunc_univ(x, *args):
    """" optimize substitution parameters one by one """
    it = args[0]
    model_name = args[1]
    sh_ac = args[2]
    subs_param_b = args[3]
    subs_param_h = args[4]
    transitions = args[5]
    forest = args[6]
    M = args[7]
    p_len = args[8]

    if it <= 3 + p_len - 1:
        idx = it - 2
        subs_param_b[idx] = x
    else:
        idx = it - 2 - p_len
        subs_param_h[idx] = x

    K = int(len(forest)/2)
    # update substitution model parameters, tree e, f
    new_forest = []
    if sh_ac == 0:   # all tree use the subs_param_b, update all tree
        for ti in range(len(forest)):
            tree_new = tree(forest[ti].D, 1, forest[ti].seqs, subs_param_b, model_name)
            new_forest.append(tree_new)
    else:  # b, h use different value
        if it <= 3 + p_len - 1:
            new_forest[0: K] = [tree(forest[t].D, 1, forest[t].seqs, subs_param_b, model_name) for t in range(K)]
            new_forest[K:2*K] = forest[K:2*K]
        else:
            new_forest[0: K] = forest[0:K]
            new_forest[K: 2*K] = [tree(forest[t].D, 1, forest[t].seqs, subs_param_h, model_name) for t in range(K, 2*K)]

    emissions = np.zeros((2*K, M))  # obtain the new emission matrix, the forest is same
    for k in range(2*K):
        emissions[k, :] = new_forest[k].get_emissions()

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def objfunc_tm(x, *args):
    # optimize the parameters of transition matrix
    if args[0] == 0:
        sb = x
        sh = args[1]
        gamma = args[2]
    elif args[0] == 1:
        sb = args[1]
        sh = x
        gamma = args[2]
    else:
        sb = args[1]
        sh = args[2]
        gamma = x

    K = args[3]
    M = args[4]
    emissions = args[5]
    transitions = get_trans_matrix(sb, sh, gamma, K)  # here K is K, not 2*K

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def tree_scale(x, *args):
    transitions = args[0]
    cu_forest = args[1]
    subs_param_h = args[2]
    K = args[3]   # the half number of the state
    M = args[4]
    model_name = args[5]

    # update necessary tree in the forest
    forest = []
    forest[0:K] = cu_forest[0:K]
    forest[K: 2 * K] = [tree(cu_forest[t].D, x, cu_forest[t].seqs, subs_param_h, model_name) for t in range(K)]

    emissions = np.zeros((2 * K, M))  # obtain the new emission matrix, the forest is same
    for kk in range(2 * K):
        emissions[kk, :] = forest[kk].get_emissions()

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def tree_branch(x, *args):
    # x is one of the tree branch length
    branch_idx = args[0]
    transitions = args[1]
    forest = args[2]
    subs_param_b = args[3]
    subs_param_h = args[4]
    K = args[5]   # the half number of the state
    M = args[6]
    N = args[7]
    all_branch = args[8]
    hscale = args[9]
    model_name = args[10]
    all_branch[branch_idx] = x

    branch_num = 2*N-2
    # generate new D; only generate the new D;
    current_idx = int(np.floor(branch_idx / branch_num))

    temp = all_branch[int(current_idx * branch_num) : int((current_idx + 1) * branch_num)]  # new branch length
    original_D = forest[current_idx].D  # original D
    new_D = gen_new_D(temp, original_D)

    # update forest list
    seqs = forest[0].seqs
    forest[current_idx] = None
    forest[current_idx+K] = None
    forest[current_idx] = tree(new_D, 1, seqs, subs_param_b, model_name)
    forest[current_idx + K] = tree(new_D, hscale, seqs, subs_param_h, model_name)

    emissions = np.zeros((2*K, M))  # obtain the new emission matrix, the forest is same
    for k in range(2*K):
        emissions[k, :] = forest[k].get_emissions()

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def get_trans_matrix(sb, sh, gamma, K):
    new_transitions = np.zeros((2 * K, 2 * K))  # the expectation transition from k to l
    f11 = sb
    f12 = (1 - sb) / float(K - 1)
    f21 = sh
    f22 = (1 - sh) / float(K - 1)
    for i in range(2 * K):  # type:
        for j in range(2 * K):
            if (i < K) and (j < K):
                if i == j:
                    new_transitions[i, j] = (1 - gamma) * f11
                else:
                    new_transitions[i, j] = (1 - gamma) * f12
            if (i < K) and (j >= K):
                new_transitions[i, j] = gamma / float(K)
            if (i >= K) and (j < K):
                new_transitions[i, j] = gamma / float(K)
            if (i >= K) and (j >= K):
                if i == j:
                    new_transitions[i, j] = (1 - gamma) * f21
                else:
                    new_transitions[i, j] = (1 - gamma) * f22

    return new_transitions


def gen_new_D(new_branch, original_D):

    temp = np.triu(original_D, k=1)
    temp.ravel()[np.flatnonzero(temp)] = new_branch
    new_D = temp+temp.T

    return new_D
