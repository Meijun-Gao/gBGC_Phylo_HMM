from numpy import *
from math import log, exp


def add_logs(x, y):
    "A fast way to add logarithms"
    if not x == 0 and not y == 0:
        return x + log(1 + exp(y - x))
    elif x == 0 and y == 0:
        print('problem with log values')
        return None
    elif x == 0:
        return y
    elif y == 0:
        return x


def Forward(transitions, emissions, M):
    # f is a array with K*M
    sequence = range(M)
    K = len(transitions[0, :])
    global begin_transition  # eventually change this??
    begin_transition = [1 / float(K) for k in range(K)]
    f = {}
    for i in range(-1, M):
        for l in range(K):
            f[l, i] = 0
            if i == -1:
                f[l, i] = log(begin_transition[l])
            else:
                for k in range(K):
                    f[l, i] = add_logs(f[l, i], f[k, i - 1] + log(transitions[k, l]))
                f[l, i] += emissions[l, sequence[i]]
    seqprob = 0
    for k in range(K):
        seqprob = add_logs(seqprob, f[k, M - 1] + log(begin_transition[k]))
    # print(seqprob)    # log likelihood score
    return f, seqprob


def Backward(transitions, emissions, M):
    sequence = range(M)
    K = len(transitions[0, :])
    b = {}
    for i in range(M-1,-1,-1):
        for k in range(K):
            if i == M - 1:
                b[k, M - 1] = log(begin_transition[k])
            else:
                b[k, i] = 0
                for l in range(K):
                    b[k, i] = add_logs(b[k, i], log(transitions[k, l]) + emissions[l, sequence[i + 1]] + b[l, i + 1])
    b_seqprob = 0
    for l in range(K):
        b_seqprob = add_logs(b_seqprob, b[l, 0] + log(begin_transition[l]) + emissions[l, sequence[0]])
    return b, b_seqprob


def viterbi(data, init, trans, emit):
    "computes viterbi sequence for a given observed sequence"
    inferred_states = []
    K = len(init)
    for obs in range(len(data[1, :])):
        Probs = {}
        for state in range(K):
            if inferred_states == []:
                Probs[state] = init[state] * exp(emit[state, obs])
            else:
                Probs[state] = trans[inferred_states[-1], state] * exp(emit[state, obs])
        highest_probability = max(Probs.values())
        for state in range(K):
            if Probs[state] == highest_probability:
                next_hidden_state = state
                break
        inferred_states.append(next_hidden_state)
    return inferred_states
