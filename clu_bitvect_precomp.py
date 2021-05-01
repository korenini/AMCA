import settings as pset

all_bit_clu = []

idx2vect = {}
vect2idx = {}


for i in range(pset.k):
    bstr_lst = []
    bstr = bin(2**i)[2:]
    bstr_prepended = bstr.zfill(pset.k)
    bstr_rev = bstr_prepended[::-1]
    bstr_lst = [bool(int(x)) for x in bstr_rev]
    all_bit_clu.append(bstr_lst)
    idx2vect[i] = bstr_lst
    vect2idx[tuple(bstr_lst)] = i
