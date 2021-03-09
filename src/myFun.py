import re
import numpy as np
from math import log, exp


def gap_stripper(alignArray):
    """ Remove gaps...if is deemed necessary ## NOTE, USUALLY THIS ISN'T DONE ANYMORE """
    for i, species in enumerate(alignArray):
        for j, letter in enumerate(species):
            if letter == 4:
                # print "Found gap...deleting column..."
                for k, specys in enumerate(alignArray):
                    alignArray[k][j] = -1
    for i, species in enumerate(alignArray):
        while -1 in species:
            alignArray[i].remove(-1)
    return alignArray


def delete_gap(input_array):
    old_mat = input_array
    for i in range(len(old_mat)):
        temp_mat = old_mat[i, :]
        col = [k for k in range(len(temp_mat)) if temp_mat[k] == 4]
        old_mat[:, col] = -1
    first_row = old_mat[0, :]
    delete_col = [k for k in range(len(first_row)) if first_row[k] == -1]
    new_mat = np.delete(old_mat, delete_col, axis=1)
    return new_mat


def get_SNP(input_array):
    # input_array  N*M
    index_del = []
    seq_len = len(input_array[0])
    for i in range(seq_len):
        temp = set(input_array[:, i])
        if len(temp) == 1:    # all taxa have the same letter in this site, delete it
            index_del.append(i)
    new_mat = np.delete(input_array, index_del, axis=1)
    print(seq_len - len(new_mat[0]))
    return new_mat


def newick2D(tree1, N):
    import ete3
    # N is the number of taxa
    total_node = 2 * N - 1
    D = np.zeros((total_node, total_node))
    begin = 0

    # add the index of inner node
    for i in range(N - 1):  # total node - (N)
        idx = tree1.index(')', begin)
        str_list = list(tree1)
        str_list.insert(idx + 1, str(N + 1 + i))
        tree1 = ''.join(str_list)
        begin = idx + 1

    tree1 = ete3.Tree(tree1, format=1)
    # if the branch length is 0, let it to 0.001 in order to avoid exception
    for i in range(total_node):
        leaf1 = str(i + 1)
        ancestor_temp = (tree1 & leaf1).get_ancestors()
        if i < total_node - 1:
            parent = ancestor_temp[0].name  # the direct ancestor
        else:
            parent = None  # the root we assume
        childs = (tree1 & leaf1).get_children()
        child_name = []
        for c in range(len(childs)):
            child_name.append(childs[c].name)

        for j in range(i + 1, total_node):  # only compute upper triangle; j from i+1
            leaf2 = str(j + 1)
            if leaf2 == parent or leaf2 in child_name:
                y = tree1.get_distance(leaf1, leaf2)
                if y == 0.0:
                    y = 0.0001
                D[i, j] = y
    D = D.T + D

    return D


def gen_new_D(new_branch, original_D):
    temp = np.triu(original_D, k=1)
    temp.ravel()[np.flatnonzero(temp)] = new_branch
    new_D = temp + temp.T

    return new_D


def D2branches(forest):
    # the index in a matrix, count from the first row all elements to next row.
    K = int(len(forest) / 2)
    N = len(forest[0].seqs[:, 0])  # taxa number
    LL = (2 * N - 2)
    branches = np.zeros(K * LL)
    for i in range(K):
        temp1 = np.triu(forest[i].D, k=1)
        branches[LL * i:LL * (i + 1)] = temp1.ravel()[np.flatnonzero(temp1)]

    return branches


def olog(x):
    if x<=0:
        print(x)
        print('calling negative log, in felsenstein module')
    else:
        return log(x)


def o_exp(x):
    try:
        if x<-14:
            return 0
        else:
            return exp(x)
    except OverflowError:
        print("Overflow in o_exp: calling", x)


def mat_log(X):
    if np.any(X == 0):
        print('0 element for prob in my_tree.py module')
    if np.any(X<0):
        print('calling negative log in my_tree.py module')
    else:
        return np.log(X)


def add_logs(x, y):
    "A fast way to add logarithms, the form of x, y is logx, logy"
    """ log(x+y) """
    if not x == 0 and not y == 0:
        return x + olog(1 + o_exp(y - x))
    elif x == 0 and y == 0:
        print('problem with log values')
        return 3500
    elif x == 0:
        return y
    elif y == 0:
        return x


def get_fasta(fasta_input, log=True):
    if isinstance(fasta_input, str):
        fasta = open(fasta_input, 'r')
    else:
        print(type(fasta_input), ': unknown filetype!')
    alignmentStringList=[line.strip() for line in fasta.readlines()]
    fasta.close()
    alignDict={}
    current_taxa=None
    for line in alignmentStringList:
        if len(line)==0:
            continue
        elif line[0]=='>':
            if log:
                print("New taxa:", line[1:])
            alignDict[line[1:]]=''
            current_taxa=line[1:]
        else:
            if log:
                print("Adding to taxa:", current_taxa)
            alignDict[current_taxa]+=line

    return alignDict
