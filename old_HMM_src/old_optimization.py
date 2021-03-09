import numpy as np
from src.training import *
from src.HMM_algos import Forward
from src.my_tree import *


def objfunc_univ(x, *args):
    """" optimize substitution parameters one by one """
    idx = args[0]  # the first parameter is scale, we fix it to 1.0
    subs_param = args[1]
    subs_param[idx] = x

    transitions = args[2]
    forest = args[3]
    M = args[4]
    model_name = args[5]
    K = len(transitions[:, 1])
    # update substitution model parameters, tree e, f
    new_forest = []
    for tree_i in forest:
        tree_new = tree(tree_i.D, 1, tree_i.seqs, subs_param, model_name)
        new_forest.append(tree_new)
    forest = new_forest

    emissions = np.zeros((K, M))  # obtain the new emission matrix, the forest is same
    for k in range(K):
        emissions[k, :] = forest[k].get_emissions()

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def objfunc_tm(x, *args):
    # optimize the parameters of transition matrix
    sb = x
    K = args[1]
    M = args[2]
    emissions = args[3]
    transitions = get_trans_matrix(sb, K)  # here K is K, not 2*K

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def tree_branch(x, *args):
    # x is one of the tree branch length    from 14

    branch_idx = args[0] - 11
    transitions = args[1]
    forest = args[2]
    subs_param = args[3]
    K = args[4]   # the half number of the state
    M = args[5]
    all_branch = args[6]
    model_name = args[7]
    all_branch[branch_idx] = x

    branch_num = int(len(all_branch) / K)
    # generate new D; only generate the new D;
    current_modified_tree = int(np.floor(branch_idx / branch_num))

    temp = all_branch[int(current_modified_tree * branch_num):int((current_modified_tree + 1) * branch_num)]  # new branch length
    original_D = forest[current_modified_tree].D  # original D
    new_D = gen_new_D(temp, original_D)

    # update forest list
    seqs = forest[0].seqs
    forest[current_modified_tree] = None
    forest[current_modified_tree] = tree(new_D, 1, seqs, subs_param, model_name)

    # uodate emission
    emissions = np.zeros((K, M))  # obtain the new emission matrix, the forest is same
    for k in range(K):
        emissions[k, :] = forest[k].get_emissions()

    forward, seq_prob = Forward(transitions, emissions, M)
    return -1 * seq_prob


def get_trans_matrix(sb, K):
    new_transitions = np.zeros((K, K))  # the expectation transition from k to l
    f11 = sb
    f12 = (1 - sb) / float(K - 1)
    for i in range(K):  # type:
        for j in range(K):
            if i == j:
                new_transitions[i, j] = f11
            else:
                new_transitions[i, j] = f12

    return new_transitions


def gen_new_D(new_branch, original_D):

    temp = np.triu(original_D, k=1)
    temp.ravel()[np.flatnonzero(temp)] = new_branch
    new_D = temp+temp.T

    return new_D
