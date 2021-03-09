import sys
import numpy as np
from scipy import optimize
from src.my_tree import *
from src.HMM_algos import Forward, Backward, viterbi
from old_HMM_src.old_optimization import *
import os
import random
import copy
import time

trial = 0
like = [[]]
for i in range(100):
    like.append([])


def Baum_Welch(trial, recobj, tree_set, learnk, input_subs_b=None,
               K=None, model_name='HKY', prefix='', writeLog=False):
    if model_name == 'GTR':
        p_len = 10
    else:
        p_len = 5

    if input_subs_b is None:
        input_subs_param_b = np.array([random.uniform(0.001, 1.0) for t in range(p_len)])
    else:
        input_subs_param_b = np.array(input_subs_b)

    try:
        sequences = recobj.alignment
    except:
        sequences = recobj
    like.append([])
    eps = .0005
    out_maxround = 8000
    sb = random.uniform(0.001, 1.0)
    N = len(sequences[:, 1])  # taxa number
    M = len(sequences[1, :])  # alignment length
    subs_param_b = np.array([1.0] * (p_len+1))
    subs_param_b[1:(p_len+1)] = input_subs_param_b

    transitions = get_trans_matrix(sb, K)

    D = np.zeros((K, 2 * N - 1, 2 * N - 1))
    for tree_i in range(len(tree_set)):
        D[tree_i, :, :] = newick2D(tree_set[tree_i][:-1], N)
    test_forest = [tree(D[i, :, :], 1, sequences, subs_param_b, model_name) for i in range(K)]
    forest = test_forest
    # Finished with tree initialization

    old_seq_prob = -np.inf
    emits = np.zeros((K, M))
    for k in range(K):
        emits[k, :] = forest[k].get_emissions()
    emissions = emits
    # Start iterative training of HMM
    iterations = 0
    omit, begin_seq_prob = Forward(transitions, emissions, M)
    if not markov_test(transitions):
        sys.stderr.write('Transitions not summing to 1')

    ''' training process, until convergence diff'''
    while True:
        forward, seq_prob = Forward(transitions, emissions, M)
        # learn the substitution parameters, use the
        diff = seq_prob - old_seq_prob
        like[trial].append(seq_prob)
        if eps > diff >= 0:
            print('convergence outside')
            break
        if abs(diff / seq_prob) < .001:
            if writeLog:
                sys.stderr.write('convergence by fractional increase criterion\n')
            break

        backward, b_seq_prob = Backward(transitions, emissions, M)
        # check for oscillatory behaviour:
        try:
            x = like[trial][-1] - like[trial][-2] - (like[trial][-3] - like[trial][-4])
            if abs(x) < .01 and like[trial][-1] > like[trial][-2]:
                sys.stderr.write('Convergence deduced by oscillatory behaviour\n')
                like[trial].append(seq_prob)
                break
        except:  # it's somewhat OK if something goes wrong here, usually among the first few iters, it can't
            pass

        # Check if forward and backward recursions gave same sequence prob. IMPORTANT!
        if abs(exp(seq_prob) - exp(b_seq_prob)) > .01:
            sys.stderr.write('Forward, backward calculations are inconsistent!  Check calculations\n')

        t_eps = 0.0005
        t_maxround = 8000
        t_round = 0
        """ optimize other parameters one by one and update emission and forest"""
        previous_seq_prob = seq_prob
        current_seq_prob = -np.infty
        round_seq_prob = previous_seq_prob

        current_transitions = transitions
        current_subs_param = subs_param_b
        current_forest = forest
        current_emission = emissions
        current_branches = D2branches(forest)
        current_sb = sb
        # optimize one by one
        inner_like = []
        while True:
            print('begin optimize')
            inner_like.append(current_seq_prob)
            t_diff = current_seq_prob - previous_seq_prob
            if t_eps > t_diff >= 0:
                print('convergence')
                break
            if t_round > t_maxround:
                print('exceed the maxround')
                break
            # check for oscillatory behaviour:
            try:
                x = like[trial][-1] - like[trial][-2] - (like[trial][-3] - like[trial][-4])
                if abs(x) < .01 and like[trial][-1] > like[trial][-2]:
                    print('Convergence deduced by oscillatory behaviour\n')
                    break
            except:  # it's somewhat OK if something goes wrong here, usually among the first few iters, it can't
                pass
            for it in range(learnk):
                if it == 0:   #
                    print('optimize sb')
                    bounds = (0.001, 1.0)
                    solnx = optimize.minimize_scalar(objfunc_tm, bounds=bounds, method='bounded',
                                                     args=(it, K, M, current_emission),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})

                    # if the likelihood increase, update the paramteters, else don't update
                    temp_transitions = get_trans_matrix(solnx.x, K)
                    t_forward, temp_seq_prob = Forward(temp_transitions, current_emission, M)
                    if temp_seq_prob > round_seq_prob:
                        current_sb = solnx.x
                        round_seq_prob = temp_seq_prob
                        current_transitions = temp_transitions
                elif 1 <= it <= p_len:
                    print('optimize subs params')
                    if model_name == 'HKY' and it == 5:
                        bounds = (0.001, 10)  # for beta
                    else:
                        bounds = (0.001, 1.0)
                    solnx = optimize.minimize_scalar(objfunc_univ, bounds=bounds, method='bounded',
                                                     args=(
                                                         it, current_subs_param, transitions,
                                                         current_forest,
                                                         M, model_name),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})

                    # if the likelihood increase, update the paramteters
                    temp_param = current_subs_param
                    temp_param[it] = solnx.x
                    # update substitution model parameters, tree e, f
                    temp_forest = []
                    for tree_i in current_forest:
                        tree_new = tree(tree_i.D, 1, tree_i.seqs, temp_param, model_name)
                        temp_forest.append(tree_new)
                    temp_emissions = np.zeros((K, M))  # obtain the new emission matrix, the forest is same
                    for kk in range(K):
                        temp_emissions[kk, :] = temp_forest[kk].get_emissions()
                    temp_forward, temp_prob = Forward(current_transitions, temp_emissions, M)

                    if temp_prob > round_seq_prob:
                        round_seq_prob = temp_prob
                        current_subs_param = temp_param
                        current_forest = temp_forest
                        current_emission = temp_emissions

                elif it >= 11:  # from 11
                    # branch length in back, denotes as l,  in hotspot can be denoted as h*l;  h is the scale
                    print('begin optimize branch length')
                    bounds = (0.001, 1e3)
                    soln_b = optimize.minimize_scalar(tree_branch, bounds=bounds, method='bounded',
                                                      args=(
                                                          it, current_transitions, current_forest,
                                                          current_subs_param,
                                                          K, M, current_branches, model_name),
                                                      options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})
                    branch_idx = it - 11
                    temp_forest = current_forest
                    temp_branches = current_branches
                    temp_branches[branch_idx] = soln_b.x

                    branch_num = int(len(current_branches) / K)
                    current_idx = int(np.floor(branch_idx / branch_num))

                    temp = temp_branches[current_idx * branch_num: (current_idx + 1) * branch_num]  # new branch length
                    original_D = temp_forest[current_idx].D  # original D
                    new_D = gen_new_D(temp, original_D)

                    # update forest list
                    temp_forest[current_idx] = None
                    temp_forest[current_idx] = tree(new_D, 1, sequences, current_subs_param, model_name)

                    # update emission
                    temp_emissions = np.zeros((K, M))  # obtain the new emission matrix, the forest is same
                    for k in range(K):
                        temp_emissions[k, :] = temp_forest[k].get_emissions()

                    temp_forward, temp_prob = Forward(current_transitions, temp_emissions, M)
                    if temp_prob > round_seq_prob:
                        round_seq_prob = temp_prob
                        current_forest = temp_forest
                        current_emission = temp_emissions
                        current_branches = temp_branches

            previous_seq_prob = current_seq_prob
            current_seq_prob = round_seq_prob
            t_round = t_round + 1

        sb = current_sb
        subs_param_b = current_subs_param
        forest = current_forest  # update the substitution parameters, emission. the topology doesn't change
        emissions = current_emission
        transitions = current_transitions
        old_seq_prob = seq_prob

        iterations += 1

        if not markov_test(transitions):
            sys.stderr.write('Transitions not summing to 1, after iterations: ' + str(iterations) + '\n')
            ### Loop until likelihood changes a small amount, then return maximized matrices
        if iterations >= out_maxround:
            break

    forward, seq_prob = Forward(transitions, emissions, M)
    # update posterior in each site
    backward, b_seq_prob = Backward(transitions, emissions, M)
    posteriors = np.zeros((K, M))
    for k in range(K):
        forest[k].posteriors = []
        for m in range(M):
            post = (exp(forward[k, m] + backward[k, m] - seq_prob))  # the prob. that column came from tree k
            forest[k].posteriors.append(post)
            posteriors[k, m] = post
    final_subs_param = forest[0].model.get_subs_param()  # substitution model parameters
    final_trans_param = sb
    print('finish')
    # final_branch = D2branches(forest)
    return begin_seq_prob, seq_prob, transitions, emissions, forest, posteriors, final_subs_param, final_trans_param


def markov_test(array):
    good = True
    rows, cols = np.shape(array)
    for row in range(rows):
        total = sum(array[row, :])
        if abs(total - 1) > .001:
            good = False
    return good


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


def D2branches(forest):
    # the index in a matrix, count from the first row all elements to next row.
    K = len(forest)
    N = len(forest[0].seqs[:, 0])  # taxa number
    LL = (2 * N - 2)
    branches = np.zeros(K * LL)
    for i in range(K):
        temp1 = np.triu(forest[i].D, k=1)
        branches[LL * i:LL * (i + 1)] = temp1.ravel()[np.flatnonzero(temp1)]

    return branches


def newick2D(tree1, N):
    import ete3
    # N is the number of taxa
    total_node = 2 * N - 1
    D = np.zeros((total_node, total_node))
    begin = 0

    # add the index of inner node
    for i in range(N-1):   # total node - (N)
        idx = tree1.index(')', begin)
        str_list = list(tree1)
        str_list.insert(idx + 1, str(N+1 + i))
        tree1 = ''.join(str_list)
        begin = idx + 1

    tree1 = ete3.Tree(tree1, format=1)
    # if the branch length is 0, let it to 0.001 in order to avoid exception
    for i in range(total_node):
        leaf1 = str(i + 1)
        ancestor_temp = (tree1 & leaf1).get_ancestors()
        if i < total_node - 1:
            parent = ancestor_temp[0].name    # the direct ancestor
        else:
            parent = None    # the root we assume
        childs = (tree1 & leaf1).get_children()
        child_name = []
        for c in range(len(childs)):
            child_name.append(childs[c].name)

        for j in range(i + 1, total_node):   # only compute upper triangle; j from i+1
            leaf2 = str(j + 1)
            if leaf2 == parent or leaf2 in child_name:
                y = tree1.get_distance(leaf1, leaf2)
                if y == 0.0:
                    y = 0.0001
                D[i, j] = y
    D = D.T + D

    return D

