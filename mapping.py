import numpy as np
import hashlib
import settings as pset
from compact_data import BIN2HEX, BIN2BASE64, encode_bin_to_base64_clu_k


def mapfunction(clu_k):
    # Map to hash
    if pset.map_function == "md5":
        return hashlib.md5(clu_k.data.tobytes()).hexdigest()
    elif pset.map_function == "sha256":
        return hashlib.sha256(clu_k.data.tobytes()).hexdigest()
    # Map to compact data
    elif pset.map_function == "hex":
        return hashlib.sha256(clu_k.data.tobytes()).hexdigest()
    elif pset.map_function == "base64":
        return encode_bin_to_base64_clu_k(clu_k)


def make_clu_hash(clux):
    hashset = set()
    hashlist = []
    for i in range(pset.k):
        ht = mapfunction(clux[i])
        hashset.add(ht)
        hashlist.append(ht)
    return frozenset(hashset), tuple(hashlist)


def make_clu_hash_tpl(clux):
    hashlist = []
    for i in range(pset.k):
        ht = mapfunction(clux[i])
        hashlist.append(ht)
    return tuple(hashlist)


def make_clu_hash_set(clux):
    hashset = set()
    for i in range(pset.k):
        ht = mapfunction(clux[i])
        hashset.add(ht)
    return frozenset(hashset)


def make_pf_fs2tpl_hashes(clu):
    chashset2chashtpl = {}
    for c in clu:
        chash, chashtpl = make_clu_hash(c)
        chashset2chashtpl[chash] = chashtpl
    return chashset2chashtpl


def make_metadata(pclu):
    # Move history {hash: aclu with bool(1) indicating move attempt}
    rmtah = {}
    # Interchange history
    ritah = {}
    ritahset = {}
    # chashset:chashtpl (preserve order when needed) (needed to translate previous unsuccesful tranfsortmation attempts from peruted clustering format)
    chashset2chashtpl = {}
    #
    for i, clx in enumerate(pclu):
        chash, chashtpl = make_clu_hash(clx)
        # See if there is a duplicated clustering
        try:
            duplicated_clu = rmtah[chash]
            print("Duplicate found i=%s" %i)
        except KeyError:
            rmtah[chash] = np.copy(clx)
            ritah[chash] = np.copy(clx)
            ritahset[chash] = np.full(clx.shape, set(), dtype=object)
            chashset2chashtpl[chash] = chashtpl

    return rmtah, ritah, ritahset, chashset2chashtpl
