from math import exp, log
import numpy as np
from scipy.linalg import expm
from repoze.lru import lru_cache
import random


def make_zero(N):
    """
    Makes an empty list of lists, size N
    """
    R = [[0 for i in range(N)] for i in range(N)]
    return R


def make_diag(vector, option='array'):
    """
    Creates a list-list with vector as the diagonals.  Or, uses the built in Numeric funtion to do this
    """
    R = make_zero(len(vector))
    for i in range(len(vector)):
        R[i][i] = vector[i]
    if option == 'list':
        return R
    elif option == 'array':
        return R


# calculate the base frequency from data
def cal_base_freq(data):
    nucs = ['A', 'C', 'G', 'T']
    totals = [0, 0, 0, 0]
    total = 0
    # get the prior from data
    if type(data) != dict:
        new_data = {}
        for row in range(len(data[:, 1])):
            new_data[row] = data[row, :]
    else:
        new_data = data
    for seq in new_data.values():
        # figure out data format, be it strings or numbers.  Numbers are a little easier to work with...
        for nucleotide in seq:
            if nucleotide != 4:
                if type(nucleotide) == str:
                    base = nucs.index(nucleotide)
                elif type(nucleotide) == int:
                    base = nucleotide
                else:
                    base = int(nucleotide)
                try:
                    totals[base] += 1
                    total += 1
                except:
                    print("Problem:", base)
    ratios = [totals[i] / float(total) for i in range(4)]

    return ratios


def get_GTR_matrix(input_param):
    """ @t is the Exchangeability parameters, see substitution model in wiki"""
    input_param = np.array(input_param)
    scale = input_param[0]
    # e = softmax(input_param[1:5])
    e = input_param[1:5] / (sum(input_param[1:5]))
    t = input_param[5:11] / input_param[10]
    t = t * scale
    # get the rate matrix
    R = np.array([[-(t[0] * e[1] + t[1] * e[2] + t[2] * e[3]), t[0] * e[1], t[1] * e[2], t[2] * e[3]],
                  [t[0] * e[0], -(t[0] * e[0] + t[3] * e[2] + t[4] * e[3]), t[3] * e[2], t[4] * e[3]],
                  [t[1] * e[0], t[3] * e[1], -(t[1] * e[0] + t[3] * e[1] + t[5] * e[3]), t[5] * e[3]],
                  [t[2] * e[0], t[4] * e[1], t[5] * e[2], -(t[2] * e[0] + t[4] * e[1] + t[5] * e[2])]])
    return R


def get_HKY85_matrix(input_param):
    """
    Builds the HKY85 rate matrix from the nucleotide frequencies and input alpha, beta
    """
    # transitions for beta  AG, CT;   transversions  for alpha =1
    # define rate matrix R
    scale = input_param[0]
    # e = softmax(input_param[1:5])
    e = input_param[1:5] / (sum(input_param[1:5]))
    alpha = 1
    beta = input_param[5]   # b, e in the 6 parameters
    R = np.array([[-(e[1] + beta * e[2] + e[3]),  e[1], beta * e[2], e[3]],
                  [e[0], -(e[0] + e[2] + beta * e[3]), e[2], beta * e[3]],
                  [beta * e[0], e[1], -(beta * e[0] + e[1] + e[3]), e[3]],
                  [e[0], beta * e[1], e[2], -(e[0] + beta * e[1] + e[2])]])
    return R


def exp_diag(mat, t):
    """
    Returns a diagonal matrix with the diagonal elements taken exp()
    """
    m = len(mat)
    array = [[] for i in range(m)]
    for i in range(m):
        try:

            if not i in (0, m - 1):
                array[i] = (mat[i][:i] + (exp(mat[i][i] * t),) + mat[i][i + 1:])
            elif i == 0:
                array[i] = (exp(mat[i][i] * t),) + mat[i][i + 1:]
            elif i == m - 1:
                array[i] = mat[i][:-1] + (exp(mat[i][i] * t),)
        except:
            print(t, array[i][i], array[i][i] * t, 'Overflow here, asking for exp() of too large number')
    return (array)


def array2tuples(array):
    'make an array into a tuple of tuples that is non-mutable.  No longer have to use deepcopy'
    retobj = []
    rows = cols = len(array[:, 1])
    for i in range(rows):
        retobj.append(tuple(array[i, :]))
    return tuple(retobj)


class Model(object):
    """
    define an evolution model based on the rate matrix and nucleotide distribution
    """

    def __init__(self, R, prior, name):
        self.R = R  # rateMatrix
        self.name = name
        self.scale = 1
        self.prior = prior
        self.alphabet = ['A', 'C', 'G', 'T']
        self.alpha = 1
        self.beta = None
        self.trans_param = None  # 6 free parameters

    def getRateMatrix(self):
        return self.R

    def get_base_freq(self):
        return self.prior

    def get_hky_param(self):
        return self.beta

    def get_subs_param(self):
        if self.name == 'HKY':
            subs_param = [0] * 5
            # subs_param[0] = self.scale
            subs_param[0:4] = self.prior
            subs_param[4] = self.beta
        else:
            subs_param = [0] * 10
            # subs_param[0] = self.scale
            subs_param[0:4] = self.prior
            subs_param[4:10] = self.trans_param
        return subs_param

    def get_grt_param(self):
        return self.trans_param

    @lru_cache(maxsize=8192)
    def prob_mutation(self, t):
        afterTime = expm(self.getRateMatrix() * t)  # 4*4
        return afterTime


def add_logs(x, y):
    "A fast way to add logarithms without having to use the exp() function"
    if not x == 0 and not y == 0:
        return x + log(1 + exp(y - x))
    elif x == 0 and y == 0:
        print('problem with log values')
        return None
    elif x == 0:
        return y
    elif y == 0:
        return x


def log_prob(data, model, t):
    """
    Gives the log-probabilty of an alignment given a model of evolution and a time period
    """
    seq1, seq2 = tuple(data.keys())
    prob = model.exp_matrix(t)
    M = range(len(data[seq1]))
    logprob = 0

    # convert letters to numbers...
    for position in M:
        a, b = model.alphabet.index(data[seq1][position]), model.alphabet.index(data[seq2][position])
        logprob += log(prob[a, b])
    return logprob


# def make_model(data, base_freq, alpha, beta, trans_param, name):
def make_model(data, input_param, name):
    """
    A helper function, basically do a few model-making steps in one
    R is the rate matrix; prior is the base frequency
    subs_param, 1*11, scale, base frequency, substitution parameter, input free parameters,
    Creates a GTR model with scaled transition frequencies.
    for base frequency, normalized
    for transition frequency, normalized and scale
    Both equilibrium and pre-scaled transition frequencies need to sum to one.
    """
    if len(input_param) == 0:
        # default setting
        if name == 'HKY':
            input_param = [0.0] * 6
            input_param[0] = 1.0   # scale
            input_param[1:5] = cal_base_freq(data)
            input_param[5] = 2.0   #
        elif name == 'GTR':
            x0 = np.array([random.uniform(0.1, 1.0) for t in range(10)])
            input_param = [1.0] * 11  # input_param[0] = 1.0
            input_param[1:5] = x0[0:4]
            input_param[5:11] = x0[4:10]

    input_param = np.array(input_param)
    # priors = softmax(input_param[1:5])
    priors = input_param[1:5] / (sum(input_param[1:5]))  # if base_freq=[], use the empirical statistics
    if name == 'GTR':
        # GTR model
        R = get_GTR_matrix(input_param)
        model = Model(R, priors, 'GTR')
        model.trans_param = input_param[0] * (input_param[5:11] / input_param[10])
        model.alpha = None
        model.beta = None
        model.scale = input_param[0]
    elif name == 'HKY':
        # HKY model
        R = get_HKY85_matrix(input_param)
        model = Model(R, priors, 'HKY')
        model.trans_param = None
        model.beta = input_param[5]
        model.alpha = 1
    else:
        model = None

    return model


def softmax(x):
    """Compute the softmax in a numerically stable way. x is a vector"""
    x = x - np.max(x)
    exp_x = np.exp(x)
    softmax_x = exp_x / np.sum(exp_x)
    return softmax_x
