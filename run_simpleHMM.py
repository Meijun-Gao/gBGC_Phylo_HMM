import ete3
import sys, getopt
import numpy as np
from src.my_tree import tree2newick
import sys
import re
import os
import time
import multiprocessing as mp
from src.myFun import *
from old_HMM_src.old_optimization import *
from old_HMM_src.old_training import *


def main(argv):
    usage = '''
    General Usage:
    ./run_oldHMM.py alignmentFile.fasta all-genetrees.txt
    Options:
    -h show the basic usage of the algorithm
    -t INT number of independent trials to run
    --modelname HKY or GTR; default HKY
    --prefix PATH, where to put output files
    --SNP 0 use all site;  1 only use SNP
    '''
    """ default setting """
    begin_time = time.time()
    subs_b = None
    input_tree_set = None
    SNP = 0
    path = ''
    t = 10
    model_name = 'HKY'

    if len(argv) <= 2:
        print(usage)
        sys.exit()
    try:
        opts, args = getopt.getopt(argv[3:], "-h-t:", ["modelname=", "prefix=", "SNP="])
    except getopt.GetoptError:
        print('run_oldHMM.py alignmentFile.fasta all-genetrees_topo.txt')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt == "-t":
            t = int(arg)
        elif opt in "--modelname":
            model_name = arg
        elif opt == "--SNP":
            SNP = int(arg)
        elif opt == "--prefix":
            path = arg

    try:
        aligndict = get_fasta(sys.argv[1], log=False)
        tree_set_file = sys.argv[2]
        tree_set = None
        try:
            tree_set = open(tree_set_file, 'r').readlines()
        except:
            print('cannot open the file')
        k = len(tree_set)
        new_tree_set = []
        current_path = sys.argv[1][:-15]
        obs_path = os.getcwd()
        full_path = obs_path + '/' + current_path

        """ generate tree via RAXML"""
        for i in range(len(tree_set)):
            with open(current_path + '/temp_tree', "w+") as f:
                f.write(tree_set[i])
            if model_name == 'GTR':
                os.system('''./RAXML/raxmlHPC-AVX -p 12345 -m GTRCAT -V -f e -s %s -t %s -w %s -n infer%s''' % (
                    sys.argv[1], current_path + '/temp_tree', full_path, str(i)))
            elif model_name == 'HKY':
                os.system('''./RAXML/raxmlHPC-AVX -p 12345 -m GTRCAT --HKY -V -f e -s %s -t %s -w %s -n infer%s''' % (
                    sys.argv[1], current_path + '/temp_tree', full_path, str(i)))
            else:
                print('null model')

            tree_set_i = open(full_path + '/RAxML_result.infer' + str(i), 'r').readlines()
            t1 = ete3.Tree(tree_set_i[0])
            R = t1.get_midpoint_outgroup()  # and set it as tree outgroup
            t1.set_outgroup(R)
            new_tree_set.append(t1.write(format=5) + '\n')

        param_list = open(full_path + '/RAxML_info.infer' + str(0), 'r').readlines()
        begin_idx = 0
        for xt in range(len(param_list)):
            if 'Tree-Length:' in param_list[xt]:
                begin_idx = xt + 1
                break
        print('begin get param')
        temp1 = None
        temp2 = None
        if model_name == 'HKY':
            temp1 = [float(param_list[begin_idx + 1][param_list[begin_idx + 1].index(':') + 2:-1])]
            temp2 = [float(param_list[begin_idx + 7][param_list[begin_idx + 7].index(':') + 2:-1]),
                     float(param_list[begin_idx + 8][param_list[begin_idx + 8].index(':') + 2:-1]),
                     float(param_list[begin_idx + 9][param_list[begin_idx + 9].index(':') + 2:-1]),
                     float(param_list[begin_idx + 10][param_list[begin_idx + 10].index(':') + 2:-1])]
        elif model_name == 'GTR':
            temp1 = [float(param_list[begin_idx][param_list[begin_idx].index(':') + 2:-1]),
                     float(param_list[begin_idx + 1][param_list[begin_idx + 1].index(':') + 2:-1]),
                     float(param_list[begin_idx + 2][param_list[begin_idx + 2].index(':') + 2:-1]),
                     float(param_list[begin_idx + 3][param_list[begin_idx + 3].index(':') + 2:-1]),
                     float(param_list[begin_idx + 4][param_list[begin_idx + 4].index(':') + 2:-1]),
                     float(param_list[begin_idx + 5][param_list[begin_idx + 5].index(':') + 2:-1])]
            temp2 = [float(param_list[begin_idx + 7][param_list[begin_idx + 7].index(':') + 2:-1]),
                     float(param_list[begin_idx + 8][param_list[begin_idx + 8].index(':') + 2:-1]),
                     float(param_list[begin_idx + 9][param_list[begin_idx + 9].index(':') + 2:-1]),
                     float(param_list[begin_idx + 10][param_list[begin_idx + 10].index(':') + 2:-1])]
        else:
            print('null')
        subs_b = temp2 + temp1
        input_tree_set = new_tree_set
        os.system('''rm %s/RAxML*''' % full_path)
        os.system('''rm %s/temp_*''' % full_path)

    except:
        print("Problem: Fasta File could not be imported. Check file location (given as: ", sys.argv[
            1], "), and format.")
        print(usage)
        sys.exit()

    numtaxa = len(aligndict.keys())
    alignArray = [[] for i in range(numtaxa)]
    for i, species in enumerate(aligndict.keys()):
        for letter in aligndict[species]:
            if not letter in ['A', 'C', 'G', 'T', '-', 'N']:
                alignArray[i].append(4)
                sys.stderr.write('unknown letter encountered: ' + letter + '\n')
            else:
                alignArray[i].append(['A', 'C', 'G', 'T', '-', 'N'].index(letter))

    alignment = np.array(alignArray)
    if SNP == 1:
        alignment = get_SNP(alignment)

    begin_seq_problist = []
    seq_problist = []
    postlist = []
    forestlist = []
    paramlist = []  # subs_param
    transitlist = []  # transition matrix
    sb_list = []

    N = len(alignment[:, 1])
    if model_name == 'GTR':
        learnk = 1 + 10 + k * (2 * N - 2)
    else:  # HKY
        learnk = 1 + 5 + k * (2 * N - 2)

    """ use multiple processing """
    num_cores = t
    list_idx = list(range(t))
    pool = mp.Pool(num_cores)
    results = [pool.apply_async(Baum_Welch, args=(i, alignment, input_tree_set, learnk, subs_b,
                                                  k, model_name, path)) for i in list_idx]
    results = [p.get() for p in results]
    # elements in one result
    # begin_seq_prob, seq_prob, transitions, emissions, forest, posteriors, final_subs_param
    for i in range(t):
        begin_seq_prob = results[i][0]
        seq_prob = results[i][1]
        transitions = results[i][2]
        emissions = results[i][3]
        forest = results[i][4]
        posteriors = results[i][5]
        final_subs_param = results[i][6]
        temp_sb = results[i][7]

        forestlist.append(forest)
        begin_seq_problist.append(begin_seq_prob)
        seq_problist.append(seq_prob)
        postlist.append(posteriors.tolist())
        paramlist.append(final_subs_param)  # subs_param
        transitlist.append(transitions.tolist())
        sb_list.append(temp_sb)

    end = time.time()
    total_time = (end - begin_time) / 3600
    max_idx = np.argmax(seq_problist)
    best_begin_prob = begin_seq_problist[max_idx]
    best_prob = seq_problist[max_idx]
    best_posterior = np.array(postlist[max_idx])
    best_forest = forestlist[max_idx]
    best_param = paramlist[max_idx]
    # best_transit = transitlist[max_idx]
    best_sb = sb_list[max_idx]
    # Print out all the trees
    result = 'Optimize time:' + str(total_time) + '\n'
    result += 'Begin Prob:' + str(best_begin_prob) + '\n'
    result += 'End Prob:' + str(best_prob) + '\n'
    # result += 'Transitions:' + str(best_transit) + '\n'
    result += 'Substitution parameters:' + str(best_param) + '\n'
    result += 'Best sb:' + str(best_sb) + '\n'
    result += "Trees:\n"
    for elem in best_forest:
        speciesList = [name.split()[0] for name in aligndict.keys()]
        result += tree2newick(elem, speciesList) + '\n'
    result += '\n'

    # Write out a header for the positions an posteriors
    numtrees = best_posterior.shape[0]
    numsites = best_posterior.shape[1]
    line = 'Position, '
    for number in range(1, numtrees + 1):
        line += 'Posteriors' + str(number) + ', '

    result += line[0:-1] + '\n'

    # Write out all the positions with all the posteriors
    for position in range(0, numsites):
        line = str(position) + ', '
        for treenum in range(0, numtrees):
            line += str(best_posterior.item((treenum, position))) + ', '
        line = line[0:-1]
        line += '\n'
        result += line

    f2 = open(path + 'result.txt', 'w')
    f2.write(result)
    f2.close()


if __name__ == "__main__":

    main(sys.argv)
