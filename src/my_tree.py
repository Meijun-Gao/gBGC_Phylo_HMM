import numpy as np
import time
from math import log, exp
from copy import deepcopy
import random
from scipy.linalg import expm
from src.my_model import *
from src.myFun import *
import sys
import os

sigma = range(4)


def Prob_t(tree, t):
    """
    define evolution model
    """
    if t<=0:
        print('Calling negative time point')
    return tree.model.prob_mutation(t)


def prior(tree, a):
    # return .25
    return tree.model.prior[a]


def is_even(x):
    return not x % 2


class tree(object):
    def __init__(self, D, h_scale, seq_array, subs_param, modelname):
        # subs_param [prior, transition/tranvers]
        # D is adjacent matrix  denote the relation of all taxa and internal nodes
        self.seqs = seq_array  # sequence
        self.model = make_model(self.seqs, subs_param, modelname)
        seqs = seq_array
        self.M = len(self.seqs[1, :])  # the length of seq
        self.sigma = range(4)  # possible change, make rate-matrix dependent
        self.A = 4
        self.D = D*h_scale
        self.nodes = len(self.D[:, 1])  # 2*N-1  include internal nodes
        self.N = len(seq_array[:, 1])  # the number of taxa
        self.root = self.nodes - 1
        if is_even(len(D[:, 1])):  # to put in extra node to make a rooted tree
            self.adj = alter(self)
        else:
            self.adj = self.D
        self.edges = list_edges(self.adj)
        self.parent, self.children, self.sibling = get_relationships(self)
        if is_even(len(D[:, 1])):
            self.nodes += 1

        self.e, self.f, = Messages(self)

        self.posteriors = None

    def likelihood(self):
        '''
        calculate log-likelihood of a tree/alignment
        '''
        likelihood = 0
        for m in range(self.M):
            temp = self.prob_column(m)
            likelihood += temp

        return likelihood

    def get_emissions(self):
        return [self.prob_column(m) for m in range(self.M)]

    def relatives(self, node):
        children = list(self.children.get(node))
        parent = self.parent.get(node)
        if children and parent:  # has children and parent
            children.append(parent)
            return children
        elif children and not parent:  # has kids, but not parents (root)
            return children
        else:  # has no children i.e. leaf node
            return [self.parent[node]]

    def isLeaf(self, node):
        children = list(self.children.get(node))
        parent = self.parent.get(node)
        if parent and not children:
            return 1
        else:
            return 0

    def prob_column(self, m):
        """
        P(observed column m | tree)
        """
        base0 = np.array(self.model.get_base_freq())
        t = np.exp(self.f[m, self.root, :])  # 1*4
        pro = np.sum(t * base0)
        return np.log(pro)


# Replaced with these two, a bit of a hack but way better
def tree2newick(tree, taxa):
    newick_str = ''
    newick_str += recursive_newick(tree, tree.root, taxa)
    newick_str += ';'
    return newick_str


def recursive_newick(tree, node, alphabet):
    # Get a list of children of the node
    children = list(tree.children[node])

    # Get the distance from this node to the parent
    dist_str = ''
    if (node != tree.root):
        dist_str = str(tree.adj[node, tree.parent[node]])

    # Base case, the node is a leaf...
    if len(children) == 0:
        # Just grab the next available leaf label and distance string
        return alphabet[node] + ':' + dist_str

    # Otherwise we recurse on all the children...
    child_newick_strs = [recursive_newick(tree, x, alphabet) for x in children]

    # ...Take the result an make a coma separated list wraped in parentheses
    result = '(' + ','.join(child_newick_strs) + ')'

    # ...and put the distance on the end (if we have one) before returning
    if (dist_str != ''):
        result += ':' + dist_str
    return result


def list_rels(i, i_array):
    'lists the relatives of a given node in an adjacency matrix'
    if i > len(i_array[:, 1]):
        print('Array size mismatch!')
    neighbors = []
    if type(i_array) == list:
        print('calling list_rels on list')
    for j, entry in enumerate(i_array[i, :]):
        try:
            if j in neighbors or entry == 0 or i == j:
                continue
        except ValueError:
            print(j)
        else:
            neighbors.append(j)
    return neighbors


def list_edges(Tree):
    'returns a list of node-pairs in a tree (as matrix)'
    nodes = len(Tree[:, 1])
    edges = []
    for i in range(nodes):
        if i == nodes - 1:
            break
        for j in range(nodes):
            if Tree[i, j] == 0 or i == j:
                continue
            if not (i, j) in edges:
                if not (j, i) in edges:
                    edges.append((i, j))
    return edges


def alter(self):
    D = self.D
    root = self.root
    if len(D[:, 1]) == root:
        print('changed something')
        root -= 1
    nodes = self.nodes
    eps = .01
    kids = list_rels(root, D)
    new = self.nodes
    adj = np.zeros((self.nodes + 1, self.nodes + 1))
    adj[0:nodes, 0:-1] = D
    adj[root, new] = adj[new, root] = eps
    x, y = kids[1], kids[2]
    adj[new, x] = adj[x, new] = D[root, x]
    adj[new, y] = adj[y, new] = D[root, y]
    adj[root, x] = adj[x, root] = adj[root, y] = adj[y, root] = 0
    return adj


def prob_obs(a, b):  # leaf node observed prob
    if a == 4:  # gap probability
        return 0.999999999
    if a == 5:
        return 0.25
    if a == b:
        return 0.999999999
    else:
        return 0.000000001


def get_relationships(tree):
    """
    returns vectors of children and parent relationships in rooted tree
    """
    children = {}
    parent = {}
    sibling = {}
    relatives = {}
    queue = [tree.root]
    covered_nodes = [tree.root]
    while 1:
        if queue == []:
            break
        # print queue
        i = queue.pop(0)
        covered_nodes.append(i)
        y = list_rels(i, tree.adj)
        # print 'relatives of',i, 'are', y
        x = deepcopy(y)
        for node in y:
            if node in covered_nodes:
                x.remove(node)
                # print 'removed', node, 'as a possible child of ',i
                # if parent[i]!=node:
                #   print 'funny business'
            else:
                parent[node] = int(i)
                # print 'added', int(i),'as parent of', node
        children[i] = tuple(deepcopy(x))
        if len(x) == 2:
            boy, girl = x[0], x[1]
            sibling[boy] = girl
            sibling[girl] = boy
        # print x,'added as children of', i
        queue += x

    return parent, children, sibling


def Messages(tree, log=False):
    # seqs assumed to be global
    root = tree.root
    M = len(tree.seqs[0, :])
    N = len(tree.seqs[:, 0])
    sigma = [0, 1, 2, 3]  #
    A = len(sigma)
    nodes = tree.nodes
    f = np.zeros((M, nodes, A))
    e = np.zeros((M, nodes, A))
    # the number of different observation 4^N, N the species number
    # reduce computational load, do not need compute the same observations
    for m in range(tree.M):
        obs = [int(x) for x in tree.seqs[:, m]]  # one column of the sequence
        obs_list = [str(x) for x in obs]
        temp_e, temp_f = cal_ef(tree, "".join(obs_list))
        f[m, :, :] = temp_f
        e[m, :, :] = temp_e

    return e, f


@lru_cache(maxsize=8192)
def cal_ef(tree, observation, log=False):
    # observation is a column of seqs; when the observation is same, the outcome is same
    obs = [int(x) for x in observation]  # one column in the observe sequence
    M = len(tree.seqs[0, :])
    N = len(tree.seqs[:, 0])
    sigma = [0, 1, 2, 3]
    A = len(sigma)
    nodes = tree.nodes
    f = np.zeros((nodes, A))
    e = np.zeros((nodes, A))
    # the index of nodes counts from leaves to inner nodes, 0,-N-1 are leaf node
    covered_nodes = []
    covered_edges = []
    queue = []
    ###   Declare leaf nodes ###
    for n in range(N):
        for a in sigma:
            f[n, a] = olog(prob_obs(obs[n], a))  # probability of observing something given xn = a
        covered_nodes.append(n)
        j = tree.parent[n]
        queue.append(j)
    while 1:
        if log:
            print('cue:', queue)
        if len(queue) == 0:
            break
        p = queue.pop(0)
        kids = deepcopy(tree.children[p])
        if p in queue:
            # ask if i was added to the queue twice, i.e. both incoming messages have been made
            queue.remove(p)
        else:
            queue.append(p)
            continue
        if log:
            print(p, 'next in upward queue')

        if log:
            print(kids)
        for a in sigma:
            f[p, a] = 0
            for c in kids:
                sum = 0
                prob_mat = Prob_t(tree, tree.adj[p, c])
                for b in sigma:
                    sum = add_logs(sum, olog(prob_mat[a][b]) + f[c, b])
                if log:
                    print(sum)
                e[c, a] = sum
                f[p, a] += sum
                if log:
                    print('added ', c, p)
        if f[p, a] == 0:
            print('an entry of f is zero...')
        new_parent = tree.parent.get(p)
        if new_parent:
            queue.append(new_parent)
        covered_nodes.append(p)
        covered_edges.append((c, p))
    return e, f


def get_siblings(treeobj):
    sibs = []
    for node in range(treeobj.nodes):
        if node == treeobj.root:
            sibs.append(0)
            continue
        sibs.append(treeobj.sibling[node])
    return sibs


def get_parents(treeobj):
    parents = []
    for node in range(treeobj.nodes):
        if node == treeobj.root:
            parents.append(0)
            continue
        parents.append(treeobj.parent[node])
    return parents
