import sys
import numpy as np
from scipy import optimize
from src.my_tree import *
from src.HMM_algos import Forward, Backward, viterbi
from scipy.optimize import LinearConstraint
from src.myoptimization import *
from scipy.optimize import Bounds
from scipy.optimize import minimize
import os
from src.myFun import *
import random

trial = 0
like = [[]]
for i in range(100):
    like.append([])


def Baum_Welch(trial, recobj, tree_set, learnk, input_subs_b=None, input_subs_h=None, sh_ac=0, input_treeheight=None,
               K=None, model_name='GTR', prefix='', writeLog=False):
    # K is the number of one kind of states; For 4 taxa, K=3
    if model_name == 'GTR':
        p_len = 10
        if sh_ac == 0:
            right_i = 12
        else:
            right_i = 22
    else:
        p_len = 5
        if sh_ac == 0:
            right_i = 7
        else:
            right_i = 12

    if input_treeheight is None:
        h_scale = random.uniform(1, 20.0)  # control the substitution rate/ mutation rate
    else:
        h_scale = input_treeheight[1] / input_treeheight[0]
    if input_subs_b is None:
        input_subs_param_b = np.array([random.uniform(0.001, 1.0) for t in range(p_len)])
    else:
        input_subs_param_b = np.array(input_subs_b)

    if input_subs_h is None:
        input_subs_param_h = np.array([random.uniform(0.001, 1.0) for t in range(p_len)])
    else:
        input_subs_param_h = np.array(input_subs_h)

    try:
        sequences = recobj.alignment
    except:
        sequences = recobj
    like.append([])
    eps = .0005
    out_maxround = 8000
    sb = random.uniform(0.6, 1.0)  # sb>sh
    sh = random.uniform(0.001, sb - 0.02)
    gamma = random.uniform(0.001, 1.0)
    trans_param = [sb, sh, gamma]
    N = len(sequences[:, 1])  # taxa number
    M = len(sequences[1, :])  # alignment length
    transitions = get_trans_matrix(sb, sh, gamma, K)
    subs_param_b = np.array([1.0] * (p_len + 1))
    subs_param_h = np.array([1.0] * (p_len + 1))
    subs_param_b[1:p_len + 1] = input_subs_param_b
    subs_param_h[1:p_len + 1] = input_subs_param_h

    option = 0
    forest = None
    if option == 1:
        print('try some different things')
    else:  # get gene trees from arguments
        D = np.zeros((K, 2 * N - 1, 2 * N - 1))  # all nodes
        for tree_i in range(len(tree_set)):
            D[tree_i, :, :] = newick2D(tree_set[tree_i][:-1], N)
        test_forest = [tree(D[i, :, :], 1, sequences, subs_param_b, model_name) for i in range(K)]
        test_forest[K: 2 * K] = [tree(D[i, :, :], h_scale, sequences, subs_param_h, model_name) for i in range(K)]
        forest = test_forest
    if writeLog:
        sys.stderr.write('Finished with tree initialization\n')

    old_seq_prob = -np.inf
    emits = np.zeros((2 * K, M))
    for k in range(2 * K):
        emits[k, :] = forest[k].get_emissions()
    emissions = emits
    # Start iterative training of HMM
    iterations = 0
    omit, begin_seq_prob = Forward(transitions, emissions, M)
    if not markov_test(transitions):
        sys.stderr.write('Transitions not summing to 1')

    ''' training process, until convergence diff'''
    while True:
        print('out begin optimize')
        forward, seq_prob = Forward(transitions, emissions, M)
        diff = seq_prob - old_seq_prob
        like[trial].append(seq_prob)
        if eps > diff >= 0:
            print('convergence outside')
            break
        if abs(diff / seq_prob) < .0001:
            print('convergence by fractional increase criterion\n')
            break

        backward, b_seq_prob = Backward(transitions, emissions, M)
        # check for oscillatory behaviour:
        try:
            x = like[trial][-1] - like[trial][-2] - (like[trial][-3] - like[trial][-4])
            if abs(x) < .01 and like[trial][-1] > like[trial][-2]:
                print('Convergence deduced by oscillatory behaviour\n')
                break
        except:  # it's somewhat OK if something goes wrong here, usually among the first few iters, it can't
            # 'look back' that far
            pass

        # Check if forward and backward recursions gave same sequence prob. IMPORTANT!
        if abs(exp(seq_prob) - exp(b_seq_prob)) > .01:
            sys.stderr.write('Forward, backward calculations are inconsistent!  Check calculations\n')

        """  Calculate transition matrix based on the parameters, gamma, sb, sh  """
        """ @sb: gene stay prob in background area; 
            @sh: gene stay prob in hotspot area;
            @gamma: the probability of transfer to different recombination area """
        t_eps = 0.0005
        t_round = 0
        t_maxround = 8000
        previous_seq_prob = seq_prob
        current_seq_prob = -np.infty
        round_seq_prob = previous_seq_prob

        sb = trans_param[0]
        sh = trans_param[1]
        gamma = trans_param[2]
        current_transitions = transitions
        current_hscale = h_scale
        current_subs_param_b = subs_param_b  # 1*11
        current_subs_param_h = subs_param_h
        current_emission = emissions
        current_forest = forest
        current_branches = D2branches(forest)  # matrix [1, K*(2N-2)]
        inner_like = []
        while True:
            print('begin optimize')
            inner_like.append(current_seq_prob)
            t_diff = current_seq_prob - previous_seq_prob
            if t_eps > t_diff >= 0:
                print('convergence inside')
                break
            if t_round > t_maxround:
                print('exceed the maxround')
                break
            # check for oscillatory behaviour:
            try:
                x = inner_like[-1] - inner_like[-2] - (inner_like[-3] - inner_like[-4])
                if abs(x) < .01 and inner_like[-1] > inner_like[-2]:
                    print('Convergence deduced by oscillatory behaviour\n')
                    break
            except:
                pass

            for it in range(learnk):  # optimize parameters one by one  sb>sh
                if it == 0:  # change: sb;  transition
                    bounds = (0.6, 1.0)
                    solnx = optimize.minimize_scalar(objfunc_tm, bounds=bounds, method='bounded',
                                                     args=(it, sh, gamma, K, M, current_emission),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})

                    # if the likelihood increase, update the paramteters, else don't update
                    temp_transitions = get_trans_matrix(solnx.x, sh, gamma, K)
                    t_forward, temp_seq_prob = Forward(temp_transitions, current_emission, M)
                    if temp_seq_prob > round_seq_prob:
                        sb = solnx.x
                        round_seq_prob = temp_seq_prob
                        current_transitions = temp_transitions
                elif it == 1:  # change: sh;  transition
                    # continue
                    bounds = (0.001, sb - 0.1)
                    solnx = optimize.minimize_scalar(objfunc_tm, bounds=bounds, method='bounded',
                                                     args=(it, sb, gamma, K, M, current_emission),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})
                    # if the likelihood increase, update the paramteters, else don't update
                    temp_transitions = get_trans_matrix(sb, solnx.x, gamma, K)
                    t_forward, temp_seq_prob = Forward(temp_transitions, current_emission, M)
                    if temp_seq_prob > round_seq_prob:
                        sh = solnx.x
                        round_seq_prob = temp_seq_prob
                        current_transitions = temp_transitions
                elif it == 2:  # change: gamma;  transition
                    # continue
                    print('optimize gamma')
                    bounds = (0.001, 1.0)
                    solnx = optimize.minimize_scalar(objfunc_tm, bounds=bounds, method='bounded',
                                                     args=(it, sb, sh, K, M, current_emission),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})
                    # if the likelihood increase, update the paramteters, else don't update
                    temp_transitions = get_trans_matrix(sb, sh, solnx.x, K)
                    t_forward, temp_seq_prob = Forward(temp_transitions, current_emission, M)
                    if temp_seq_prob > round_seq_prob:
                        gamma = solnx.x
                        round_seq_prob = temp_seq_prob
                        current_transitions = temp_transitions

                elif 3 <= it <= right_i:  # 10/20 , it from 3-12/3-22 for GTR or 3-7/3-12   # change: subs parameters; emissions, forest
                    # optimize one by one
                    if model_name == 'HKY' and it == 7:
                        bounds = (0.001, 10)  # for beta
                    elif model_name == 'HKY' and it == 12:
                        bounds = (0.001, 10)  # for beta
                    else:
                        bounds = (0.001, 1.0)
                    print('begin optimize substitution parameters')
                    solnx = optimize.minimize_scalar(objfunc_univ, bounds=bounds, method='bounded',
                                                     args=(
                                                         it, model_name, sh_ac, current_subs_param_b,
                                                         current_subs_param_h,
                                                         current_transitions, current_forest, M, p_len),
                                                     options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})
                    temp_param = None
                    if it <= 3 + p_len - 1:
                        temp_param = current_subs_param_b
                        temp_param[it - 2] = solnx.x
                    else:
                        temp_param = current_subs_param_h
                        temp_param[it - 2 - p_len] = solnx.x

                    # if the likelihood increase, update the paramteters
                    # update substitution model parameters, tree e, f
                    temp_forest = []
                    if sh_ac == 0:  # all tree use the subs_param_b, update all tree
                        for tree_i in current_forest:
                            tree_new = tree(tree_i.D, 1, tree_i.seqs, temp_param, model_name)
                            temp_forest.append(tree_new)
                    else:  # background, hotspot use different substitution model parameters
                        if it <= 3 + p_len - 1:
                            temp_forest[0: K] = [
                                tree(current_forest[t].D, 1, current_forest[t].seqs, temp_param, model_name) for t
                                in range(K)]
                            temp_forest[K:2 * K] = current_forest[K:2 * K]
                        else:
                            temp_forest[0: K] = current_forest[0:K]
                            temp_forest[K: 2 * K] = [
                                tree(current_forest[t].D, 1, current_forest[t].seqs, temp_param, model_name)
                                for t in
                                range(K, 2 * K)]

                    temp_emissions = np.zeros((2 * K, M))  # obtain the new emission matrix, the forest is same
                    for kk in range(2 * K):
                        temp_emissions[kk, :] = temp_forest[kk].get_emissions()
                    temp_forward, temp_prob = Forward(current_transitions, temp_emissions, M)

                    if temp_prob > round_seq_prob:
                        round_seq_prob = temp_prob
                        if sh_ac == 0:
                            current_subs_param_b = temp_param
                            current_subs_param_h = temp_param
                        else:
                            if it <= 3 + p_len - 1:
                                current_subs_param_b = temp_param
                            else:
                                current_subs_param_h = temp_param
                        current_forest = temp_forest
                        current_emission = temp_emissions

                elif it == right_i + 1:  # h_scale;   emissions, forest of last K states
                    # branch length in back, denotes as l,  in hotspot can be denoted as h*l;  h is the scale
                    print('begin optimize h scale')
                    bounds = (1, 100)
                    soln_h = optimize.minimize_scalar(tree_scale, bounds=bounds, method='bounded',
                                                      args=(
                                                          current_transitions, current_forest, current_subs_param_h, K,
                                                          M, model_name),
                                                      options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})

                    temp_forest = []
                    temp_forest[0:K] = current_forest[0:K]
                    temp_forest[K: 2 * K] = [
                        tree(current_forest[t].D, soln_h.x, current_forest[t].seqs, current_subs_param_h, model_name)
                        for t in
                        range(K)]

                    temp_emissions = np.zeros((2 * K, M))  # obtain the new emission matrix, the forest is same
                    for kk in range(2 * K):
                        temp_emissions[kk, :] = temp_forest[kk].get_emissions()

                    temp_forward, temp_prob = Forward(current_transitions, temp_emissions, M)
                    if temp_prob > round_seq_prob:
                        round_seq_prob = temp_prob
                        current_forest = temp_forest
                        current_emission = temp_emissions
                        current_hscale = soln_h.x

                elif it > right_i + 1:  # from 14, 24
                    # branch length in back, denotes as l,  in hotspot can be denoted as h*l;  h is the scale
                    print('begin optimize branch length')
                    bounds = (0.001, 1e3)
                    branch_idx = it - (right_i + 2)
                    soln_b = optimize.minimize_scalar(tree_branch, bounds=bounds, method='bounded',
                                                      args=(
                                                          branch_idx, current_transitions, current_forest,
                                                          current_subs_param_b, current_subs_param_h,
                                                          K, M, N, current_branches, current_hscale, model_name),
                                                      options={'disp': 0, 'maxiter': 5000, 'xatol': 1e-04})

                    temp_forest = current_forest
                    temp_branches = current_branches
                    temp_branches[branch_idx] = soln_b.x

                    branch_num = 2 * N - 2
                    current_idx = int(np.floor(branch_idx / branch_num))

                    temp = temp_branches[current_idx * branch_num: (current_idx + 1) * branch_num]  # new branch length
                    original_D = temp_forest[current_idx].D  # original D
                    new_D = gen_new_D(temp, original_D)

                    # update forest list
                    temp_forest[current_idx] = None
                    temp_forest[current_idx + K] = None
                    temp_forest[current_idx] = tree(new_D, 1, sequences, current_subs_param_b, model_name)
                    temp_forest[current_idx + K] = tree(new_D, current_hscale, sequences,
                                                        current_subs_param_h, model_name)

                    temp_emissions = np.zeros((2 * K, M))
                    for k in range(2 * K):
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

        trans_param = [sb, sh, gamma]
        subs_param_b = current_subs_param_b
        subs_param_h = current_subs_param_h
        h_scale = current_hscale
        forest = current_forest  # update the substitution parameters, emission. the topology doesn't change
        emissions = current_emission
        transitions = current_transitions
        old_seq_prob = seq_prob
        iterations += 1

        if not markov_test(transitions):
            sys.stderr.write('Transitions not summing to 1, after iterations: ' + str(iterations) + '\n')
        if iterations >= out_maxround:
            print('exceed maxround outside')
            break

    # finish optimization
    forward, seq_prob = Forward(transitions, emissions, M)
    backward, b_seq_prob = Backward(transitions, emissions, M)
    # update posterior in each site
    posteriors = np.zeros((2 * K, M))
    for k in range(2 * K):
        forest[k].posteriors = []
        for m in range(M):
            post = (exp(forward[k, m] + backward[k, m] - seq_prob))  # the prob. that column came from tree k
            forest[k].posteriors.append(post)
            posteriors[k, m] = post
    final_subs_param_b = forest[0].model.get_subs_param()  # substitution model parameters
    final_subs_param_h = forest[K + 1].model.get_subs_param()
    final_trans_param = trans_param
    final_hscale = h_scale
    print('finish')
    return begin_seq_prob, seq_prob, transitions, emissions, forest, posteriors, final_subs_param_b, final_subs_param_h, final_trans_param, final_hscale


def markov_test(array):
    good = True
    rows, cols = np.shape(array)
    for row in range(rows):
        total = sum(array[row, :])
        if abs(total - 1) > .001:
            good = False
    return good


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
