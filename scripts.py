import json
import numpy as np
from random import *

ltok=4
red_aa_dict = {'Z': 'E', 'B': 'D', 'Q': 'E', 'N': 'D', 'U': 'C', 'M': 'L'}
word_dict = {'[PAD]': 0, '[CLS]': 1, '[SEP]': 2, '[GAP]': 3}
len_word_dict = len(word_dict)
wdc = set(list(range(len_word_dict)))

def readfile(jobfile,type):
    inputfile = open(jobfile, "r")
    if type=='t':
        text = inputfile.read()
    if type=='l':
        text = inputfile.readlines()
    inputfile.close()
    return text

def writefile(name, text):
    ff = open(name, 'w')
    ff.write(text)
    ff.close()

def readjson(jobfile):
    with open(jobfile, 'r') as fl:
        data = json.load(fl)
    return data

def writejson(jobfile, data):
    with open(jobfile, 'w') as fl:
        json.dump(data, fl)

def make_ranges(data, k):
    ranges = list(range(0, len(data), k))
    if len(data) - ranges[-1] > k/2:
        ranges = [[ranges[i], ranges[i + 1]] for i in range(0, len(ranges) - 1)] + [[ranges[-1], len(data)]]
    else:
        ranges = [[ranges[i], ranges[i + 1]] for i in range(0, len(ranges) - 2)] + [[ranges[-2], len(data)]]
    chunks = [data[r[0]:r[1]] for rdx, r in enumerate(ranges)]
    return chunks

def rC_red_aa(a):
    m = ''
    for x in a:
        if x.islower():
            m += 'C'
        elif x in red_aa_dict:
            m += red_aa_dict[x]
        else:
            m += x
    return m

def  modify_dict(x):
    x = {**word_dict, **x}
    return x

def generate_inputs(aa, frd, word_dict_inp, word_dict_out, word_dict_outr, max_len, k_prob):
 
    je = set(['X'*m + 'J' + 'X'*(ltok-m-1) for m in range(ltok)])

    def add_sep(x):
        if len(je.intersection(set(x))) > 0:
            x_idx = [ax for ax, a in enumerate(x) if a not in je]
            x = [a if ax in x_idx else '[GAP]' for ax, a in enumerate(x)]
        return x

    def make_fr_i(tokens, thr):
        tokens_fr = []
        tokens_prob = []
        for ix, x in enumerate(tokens):
            if x < len_word_dict:
                t = x
                prob = 1
            else:
                tn = choice(['HHHH', 'CCCC', 'EEEE'])
                t = word_dict_out[tn]
                prob = 0
                y = [a for a in [ix - 4, ix + 4] if (a >= 0) and (a < len(tokens))]
                t_fr1 = set([tokens[a] for a in y])
                t_fr = t_fr1 - wdc
                if str(x) in frd:
                    e = frd[str(x)]
                    e = [e[str(n)] for n in t_fr if str(n) in e]
                    if len(e) > 0:
                        t = {}
                        t_dict = {}
                        for d in e:
                            for tok in d:
                                ltok = word_dict_outr[int(tok)]
                                if ltok in t_dict:
                                    t_dict[ltok] += d[tok]
                                else:
                                    t_dict[ltok] = d[tok]
                        for tok in t_dict:
                            t_ind = word_dict_out[tok]

                            if t_ind not in t:
                                t[t_ind] = t_dict[tok]
                            else:
                                t[t_ind] += t_dict[tok]
                        if len(t) > 0:
                            t_keys = [a for a in t]
                            if len(t_keys) == 1:
                                if random() <= 1.0:
                                    t_tok = t_keys[0]
                                    prob = 1
                                    t = t_tok
                                else:
                                    prob = random()/2
                                    t = randint(len(word_dict), len(word_dict_out)-1)
                            else:
                                prob_sum = sum([t[a] for a in t])
                               
                                if len(t_keys) <= thr:
                                    t_tok = choice(t_keys)
                                    prob = t[t_tok]
                                    t = t_tok
                                else:
                                    t = sorted([[a, t[a]] for a in t], key =lambda x: x[1])
                                    t = choice(t[-thr:])
                                    t, prob = t[0], t[1]
                                if prob_sum != 0:
                                    prob = round(prob / prob_sum, 2)
            #prob = min(0.5, prob)
            tokens_fr.append(t)
            tokens_prob.append(prob)
        return tokens_fr, tokens_prob

    batch = []
    aa_splitted = add_sep([aa[x:x + ltok] for x in range(len(aa) - ltok + 1)])

    aa_splitted_tokens = []
    for x in aa_splitted:
        if x in word_dict_inp:
            aa_splitted_tokens.append(word_dict_inp[x])
        else:
            x_pos = [ax for ax in range(ltok) if x[ax] == 'X']
            xi = 'LLLL'
            if x_pos == []:
                 
                x_vars = [i for i in word_dict_inp if (x[:3] == i[:3]) or (x[1:] == i[1:])]
                if len(x_vars) > 0:
                    xi = choice(x_vars)
                   
            else:
                not_x_pos = [ax for ax in range(ltok) if ax not in x_pos]
                x_vars = [i for i in word_dict_inp if len([ax for ax in not_x_pos if x[ax] == i[ax]]) == len(not_x_pos)]
                if len(x_vars) > 0:
                    xi = choice(x_vars)
            print(x, xi)
            aa_splitted_tokens.append(word_dict_inp[xi])
                

    thr = k_prob
    m_fr, p_fr = make_fr_i(aa_splitted_tokens, thr)
    tokens_a4, tokens_fr4, tokens_fr_prob4 = [
        [[m for mx, m in enumerate(a[r:]) if mx % ltok == 0] for r in range(ltok)] for a in [aa_splitted_tokens, m_fr, p_fr]]
    for r in range(ltok):
        tokens_a, tokens_fr, tokens_fr_prob = tokens_a4[r], tokens_fr4[r], tokens_fr_prob4[r]
        tokens_a, tokens_fr, tokens_fr_prob = tokens_a[:max_len-2], tokens_fr[:max_len-2], tokens_fr_prob[:max_len-2] 
        tokens_a  = [1] + tokens_a + [2]
        tokens_fr = [1] + tokens_fr + [2]
        tokens_fr_prob = [1] + tokens_fr_prob + [1]
        ln_seq = len(tokens_a)
        lt = max_len - ln_seq

        tokens_a.extend([0] * lt)
        tokens_fr.extend([0] * lt)
        tokens_fr_prob.extend([0] * lt)

        batch.append([tokens_a, tokens_fr, tokens_fr_prob])
    batch = [np.array([batch[x][k] for x in range(len(batch))]) for k in range(len(batch[0]))]
    return batch

