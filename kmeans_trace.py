import numpy as np
import random
from copy import deepcopy
import settings as pset
import get_data as gd
import compact_data as cdat
import ecs_smpl as ecs
from clu_bitvect_precomp import all_bit_clu, idx2vect

import mapping


def randominit():
    lclu = np.empty([pset.k, gd.nunits], dtype=bool)
    for i in range(gd.nunits):
        lclu[:, i] = all_bit_clu[random.randint(0, pset.k-1)]
    return(lclu)


def lkmbv(crit, crit_optimal_clx, is_random):

    flag_continue = True

    itercounter = 0

    lclu = np.copy(crit_optimal_clx)
    units = list(range(gd.nunits))
    clulst = list(range(pset.k))

    trace_clu = []
    trace_wss = []

    cur_crit_wss = ecs.ecswss_bvect_scrit(lclu, crit)

    while flag_continue:

        itercounter += 1
        if itercounter >= pset.iter_kmeans_max:
            print("Transformation limit reached. Try increasing 'iter_kmeans_max'.")
            break

        csd_val = cdat.cluster_size_all(lclu)
        csd_tf = csd_val <= pset.mspcs
        flag_clustersize = any(csd_tf)

        flag_continue = False

        if pset.doshuffle:
            random.shuffle(units)

        for mvun in units:
            clu_from = np.copy(lclu[:, mvun])

            if pset.doshuffle:
                random.shuffle(clulst)
            else:
                pass

            for clk in clulst:

                clu_to = idx2vect[clk]
                if all(clu_from == clu_to):
                    continue
                elif flag_clustersize:
                    if any(clu_from & csd_tf):
                        continue
                else:
                    pass

                new_wss, flag_improvement = ecs.ecswss_bvect_scrit_imprv(lclu, cur_crit_wss, crit, mvun, clu_to)

                if flag_improvement:
                    flag_continue = True
                    lclu[:, mvun] = np.copy(clu_to)
                    cur_crit_wss = deepcopy(new_wss)

                    print("Improvement in single criterion MOVE. WSS: %s" % str(cur_crit_wss))

                    if not is_random:
                        trace_clu.append(np.copy(lclu))
                        trace_wss.append(ecs.allecswss_bvect(lclu))
                    break
                else:
                    pass
            if flag_continue:
                break

    if is_random:
        return lclu, ecs.allecswss_bvect(lclu)
    else:
        return trace_clu, trace_wss


