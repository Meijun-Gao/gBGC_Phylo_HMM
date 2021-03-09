import sys, getopt
import ete3
import numpy as np
from src.my_tree import *
import sys
import re
import os
import time
from src.myFun import *
import multiprocessing as mp
from src.training import *


def main(argv):
    usage = '''
    General Usage:
    ./run_HMM.py alignmentFile.fasta all-genetrees.txt
    Options:
    -h show the basic usage of the algorithm
    -t INT number of independent trials to run
    --sh_ac sh activate; 1: background and hotspot use different substitution model parameters
    --modelname HKY or GTR; default HKY
    --prefix PATH, where to put output files
    --SNP 0 use all site;  1 only use SNP
    '''
    """ default setting """
    begin_time = time.time()
    subs_h = None
    subs_b = None
    tree_height = None
    SNP = 0
    path = ''
    t = 10
    model_name = 'HKY'
    sh_ac = 1

    if len(argv) <= 2:
        print(usage)
        sys.exit()
    try:
        opts, args = getopt.getopt(argv[3:], "-h-t:", ["modelname=", "sh_ac=", "prefix=", "SNP="])
    except getopt.GetoptError:
        print('run_HMM.py alignmentFile.fasta all-genetrees_topo.txt')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt == "-t":
            t = int(arg)
        elif opt in "--modelname":
            model_name = arg
        elif opt == "--sh_ac":
            sh_ac = int(arg)
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
        os.system('''rm %s/RAxML*''' % full_path)
        os.system('''rm %s/temp_*''' % full_path)
        print('begin raxml')
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
                print('null')
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
        subs_h = subs_b
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
    param_blist = []  # subs_param
    param_hlist = []  # subs_param
    trans_param_list = []  # three parameters for hotspot
    transitlist = []  # transition matrix
    hscale_list = []
    branch_list = []

    N = len(alignment[:, 1])
    learnk = 0
    if model_name == 'GTR':
        if sh_ac == 0:
            learnk = 13 + k * (2 * N - 2) + 1
        else:
            learnk = 23 + k * (2 * N - 2) + 1
    elif model_name == 'HKY':
        if sh_ac == 0:
            learnk = 3 + 5 + k * (2 * N - 2) + 1
        else:
            learnk = 3 + 5 + 5 + k * (2 * N - 2) + 1
    else:
        print('null')

    """ use multiple processing """
    num_cores = t
    list_idx = list(range(t))
    pool = mp.Pool(num_cores)
    results = [pool.apply_async(Baum_Welch, args=(i, alignment, input_tree_set, learnk, subs_b,
                                                  subs_h, sh_ac, tree_height, k, model_name, path)) for i in list_idx]
    results = [p.get() for p in results]

    for i in range(t):
        begin_seq_prob = results[i][0]
        seq_prob = results[i][1]
        transitions = results[i][2]
        emissions = results[i][3]
        forest = results[i][4]
        posteriors = results[i][5]
        final_subs_param_b = results[i][6]
        final_subs_param_h = results[i][7]
        trans_param = results[i][8]  # transition matrix param, s_b, s_h, gamma
        final_hscale = results[i][9]

        forestlist.append(forest)
        begin_seq_problist.append(begin_seq_prob)
        seq_problist.append(seq_prob)
        postlist.append(posteriors.tolist())
        param_blist.append(final_subs_param_b)  # subs_param
        param_hlist.append(final_subs_param_h)  # subs_param
        trans_param_list.append(trans_param)
        transitlist.append(transitions.tolist())
        hscale_list.append(final_hscale)

    end = time.time()
    total_time = (end - begin_time) / 3600
    max_idx = np.argmax(seq_problist)
    best_begin_prob = begin_seq_problist[max_idx]
    best_prob = seq_problist[max_idx]
    best_posterior = np.array(postlist[max_idx])
    best_forest = forestlist[max_idx]
    best_param_b = param_blist[max_idx]
    best_param_h = param_hlist[max_idx]
    best_trans_param = trans_param_list[max_idx]
    best_transit = transitlist[max_idx]
    best_hscale = hscale_list[max_idx]
    # best_branch = branch_list[max_idx]
    # Print out all the info
    result = 'Optimize time:' + str(total_time) + '\n'
    result += 'Begin Prob:' + str(best_begin_prob) + '\n'
    result += 'End Prob:' + str(best_prob) + '\n'
    result += 'Substitution parameters back:' + str(best_param_b) + '\n'
    result += 'Substitution parameters hot:' + str(best_param_h) + '\n'
    result += 'Transition parameters:' + str(best_trans_param) + '\n'
    result += 'hscale parameters:' + str(best_hscale) + '\n'
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
