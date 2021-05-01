import numpy as np
import settings as pset
import get_data as gd

# Sextets for Base64
wordlen = 6


BIN2HEX = {'0000':'0', '0001':'1', '0010':'2', '0011':'3',
           '0100':'4', '0101':'5', '0110':'6', '0111':'7',
           '1000':'8', '1001':'9', '1010':'a', '1011':'b',
           '1100':'c', '1101':'d', '1110':'e', '1111':'f' }

HEX2BIN = {'0':'0000', '1':'0001', '2':'0010', '3':'0011',
           '4':'0100', '5':'0101', '6':'0110', '7':'0111',
           '8':'1000', '9':'1001', 'a':'1010', 'b':'1011',
           'c':'1100', 'd':'1101', 'e':'1110', 'f':'1111' }

BIN2BASE64 = {'000000':'A', '000001':'B', '000010':'C', '000011':'D',
              '000100':'E', '000101':'F', '000110':'G', '000111':'H',
              '001000':'I', '001001':'J', '001010':'K', '001011':'L',
              '001100':'M', '001101':'N', '001110':'O', '001111':'P',
              '010000':'Q', '010001':'R', '010010':'S', '010011':'T',
              '010100':'U', '010101':'V', '010110':'W', '010111':'X',
              '011000':'Y', '011001':'Z', '011010':'a', '011011':'b',
              '011100':'c', '011101':'d', '011110':'e', '011111':'f',
              '100000':'g', '100001':'h', '100010':'i', '100011':'j',
              '100100':'k', '100101':'l', '100110':'m', '100111':'n',
              '101000':'o', '101001':'p', '101010':'q', '101011':'r',
              '101100':'s', '101101':'t', '101110':'u', '101111':'v',
              '110000':'w', '110001':'x', '110010':'y', '110011':'z',
              '110100':'0', '110101':'1', '110110':'2', '110111':'3',
              '111000':'4', '111001':'5', '111010':'6', '111011':'7',
              '111100':'8', '111101':'9', '111110':'+', '111111':'/'}


BASE642BIN = {'A':'000000', 'B':'000001', 'C':'000010', 'D':'000011',
              'E':'000100', 'F':'000101', 'G':'000110', 'H':'000111',
              'I':'001000', 'J':'001001', 'K':'001010', 'L':'001011',
              'M':'001100', 'N':'001101', 'O':'001110', 'P':'001111',
              'Q':'010000', 'R':'010001', 'S':'010010', 'T':'010011',
              'U':'010100', 'V':'010101', 'W':'010110', 'X':'010111',
              'Y':'011000', 'Z':'011001', 'a':'011010', 'b':'011011',
              'c':'011100', 'd':'011101', 'e':'011110', 'f':'011111',
              'g':'100000', 'h':'100001', 'i':'100010', 'j':'100011',
              'k':'100100', 'l':'100101', 'm':'100110', 'n':'100111',
              'o':'101000', 'p':'101001', 'q':'101010', 'r':'101011',
              's':'101100', 't':'101101', 'u':'101110', 'v':'101111',
              'w':'110000', 'x':'110001', 'y':'110010', 'z':'110011',
              '0':'110100', '1':'110101', '2':'110110', '3':'110111',
              '4':'111000', '5':'111001', '6':'111010', '7':'111011',
              '8':'111100', '9':'111101', '+':'111110', '/':'111111'}


def clu2bitvect(clu_current):
    bitvect = np.zeros((pset.k, gd.nunits), dtype=bool)
    for i in range(gd.nunits):
        bitvect[:,i][(clu_current[i])] = bool(1)

    return bitvect


def bitvect2clu(bitvect):
    vect = np.array((range(pset.k+1))[1:])[:, np.newaxis]
    bitvect = bitvect * vect
    nbv = np.zeros(gd.nunits, dtype=bool)
    for x in range(pset.k):
        nbv = nbv | bitvect[x]
    nbv = nbv - 1

    return nbv


def bitclu2lstclu(bitclu):
    bitclu_int = np.copy(bitclu.astype(int))
    for row in range(bitclu.shape[0]):
        bitclu_int[row] = bitclu_int[row] * (row+1)
    clu_membership = np.sum(bitclu_int, axis=0)

    return clu_membership


def encode_bin_to_base64_clu(bitvect):
    wlm = len(bitvect[0]) % wordlen
    if wlm == 0:
        nzer = ''
    else:
        nzer = '0' * (wordlen - wlm)
    clu_id = set()

    for i in range(pset.k):
        ith_bitvect_str = nzer + ''.join(['01'[x] for x in bitvect[i].astype(int)])
        clu_id_tmp = ''
        for j in range(0, gd.nunits, wordlen):
            clu_id_tmp += BIN2BASE64[ith_bitvect_str[j:(j+wordlen)]]
        clu_id.add(clu_id_tmp)

    return frozenset(clu_id)


def encode_bin_to_base64_clu_k(bitvect_k):
    wlm = len(bitvect_k) % wordlen
    if wlm == 0:
        nzer = ''
    else:
        nzer = '0' * (wordlen - wlm)
    ith_bitvect_str = nzer + ''.join(['01'[x] for x in bitvect_k.astype(int)])
    clu_id_tmp = ''
    for j in range(0, gd.nunits, wordlen):
        clu_id_tmp += BIN2BASE64[ith_bitvect_str[j:(j+wordlen)]]

    return clu_id_tmp


def decode_base64_to_bin_clu(frzn_set):
    decoded_clu = np.full((pset.k, gd.nunits), False)
    wlm = gd.nunits % wordlen
    for e, encoded_str in enumerate(frzn_set):
        bin_str = ""
        for i in encoded_str:
            bin_str += BASE642BIN[i]
        if wlm != 0:
            bin_str = bin_str[wlm:]
        decoded_clu[e] = [bool(int(x)) for x in bin_str]

    return decoded_clu


def cluster_size_single(clu_bitvect):
    csum = np.sum(clu_bitvect.astype(int), axis=1)

    return csum


def cluster_size_all(clu_bitvect):
    #csum = np.sum(clu_bitvect, axis=1, dtype=np.float64)
    csum = np.sum(clu_bitvect.astype(int), axis=1)

    return csum


def bool2idx(clu):
    for e,x in enumerate(clu):
        if x:
            return e

