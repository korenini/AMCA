import numpy as np
import random
from copy import deepcopy
import settings as pset
import get_data as gd
import compact_data as cdat
from mapping import make_clu_hash, make_clu_hash_tpl, make_clu_hash_set
import paretof as pf
import ecs_smpl as ecs
from clu_bitvect_precomp import all_bit_clu, idx2vect, vect2idx
from do_logging import logme_name


def improve_initial_workhorse(aclux, do_interchange, do_trace):

    aclu = np.copy(aclux)

    flag_continue = True
    flag_improved = False

    units = list(range(gd.nunits))
    clulst = list(range(pset.k))

    trace_clu = []
    trace_wss = []

    cur_wss = ecs.allecswss_bvect(aclu)
    while_loopcounter = 0

    while flag_continue:

        csd_val = cdat.cluster_size_all(aclu)
        csd_tf = csd_val <= pset.mspcs

        flag_clustersize = any(csd_tf)

        flag_continue = False
        flag_break = False

        if pset.doshuffle:
            random.shuffle(units)

        while_loopcounter += 1
        print("IMPROVING AN INITIAL CLUSTERING. Starting iteration: %s" % while_loopcounter)

        for mvun in units:

            if flag_break:
                break
            clu_from = np.copy(aclu[:, mvun])
            clu_from_idx = vect2idx[tuple(clu_from)]

            if pset.doshuffle:
                random.shuffle(clulst)

            for cl in clulst:
                cltidx = cl

                if flag_break:
                    break

                clu_to = np.copy(idx2vect[cltidx])

                flag_move_allowed = bool(1)

                if clu_from_idx == cltidx:
                    continue
                elif flag_clustersize:

                    if any(clu_from & csd_tf):
                        if do_interchange:
                            flag_move_allowed = bool(0)
                        else:
                            print("Cluster is too small, transformation not allowed.")
                            continue
                else:
                    pass
                move_wss, flag_improvement_mv = ecs.allecswss_barithm(aclu, cur_wss, mvun, clu_to)
                #
                if flag_improvement_mv and flag_move_allowed:
                    flag_continue = True
                    flag_improved = True
                    aclu[:, mvun] = np.copy(clu_to)

                    csd_val = cdat.cluster_size_all(aclu)

                    cur_wss = move_wss
                    print("Improvement in MOVE. WSS: %s" % str(cur_wss))

                    if do_trace:
                        trace_clu.append(np.copy(aclu))
                        trace_wss.append(ecs.allecswss_bvect(aclu))

                    flag_break = True
                    break
                if do_interchange:
                    ac = aclu[clu_to][0]
                    int_idx = [x for x in range(gd.nunits) if ac[x] and x is not mvun]

                    if pset.doshuffle:
                        random.shuffle(int_idx)
                    for intchgun in int_idx:
                        #
                        interchange_wss, flag_improvement_intchg = ecs.allecswss_barithm_intchg(aclu, cur_wss, mvun,
                                                                                    intchgun, clu_from, clu_to)
                        if flag_improvement_intchg:
                            flag_continue = True
                            flag_improved = True

                            aclu[:, mvun] = np.copy(clu_to)
                            aclu[:, intchgun] = np.copy(clu_from)

                            csd_val = cdat.cluster_size_all(aclu)
                            cur_wss = interchange_wss

                            print("Improvement in INTERCHANGE. WSS: %s" % str(cur_wss))

                            if do_trace:
                                trace_clu.append(aclu)
                                trace_wss.append(cur_wss_sum)

                            flag_break = True
                            break
                else:
                    pass

    if do_trace:
        return(trace_clu, trace_wss, flag_improved)
    else:
        return(aclu, cur_wss, flag_improved)


def improve_initial_mv(pareto_clu, pareto_wss, chashset2chashtpl):

    pareto_clu_queue = []
    pareto_wss_queue = []

    for e, i in enumerate(pareto_clu):
        pareto_clu_queue.append(np.copy(i))
        pareto_wss_queue.append(deepcopy(pareto_wss[e]))

    do_interchange = False
    do_trace = False

    lcounter = 0

    for e, initial in enumerate(pareto_clu_queue):

        lcounter += 1
        print("Starting from a new clustering no.: %s" % lcounter)

        new_clu, new_wss, flag_improved_mv = improve_initial_workhorse(initial, do_interchange, do_trace)

        chash, chashtpl = make_clu_hash(new_clu)

        if flag_improved_mv and chash not in chashset2chashtpl.keys():

            chash, chashtpl = make_clu_hash(new_clu)
            chashset2chashtpl[chash] = chashtpl

            dominated, flag_pareto = pf.eval_paretof(pareto_wss, new_wss)

            for dom_idx in dominated:
                dhash = make_clu_hash_set(pareto_clu[dom_idx])
                del chashset2chashtpl[dhash]
                del pareto_clu[dom_idx]
                del pareto_wss[dom_idx]

            pareto_clu.append(np.copy(new_clu))
            pareto_wss.append(deepcopy(tuple(new_wss)))

            if flag_pareto:
                print("IMPROVEMENT made in MOVE")
            else:
                print("Improvement made in MOVE but no impact on Pareto front.")
        else:
            print("No improvement in MOVE transformation.")

    #print("Move-only improvement of Pareto front of clusterings obtained via tracing.")
    print("Number of all Pareto clusterings (candidates): %s" % len(pareto_clu))

    return pareto_clu, pareto_wss


def improve_initial_intchg(pareto_clu, pareto_wss, chashset2chashtpl):

    pareto_clu_queue = []
    pareto_wss_queue = []

    for e, i in enumerate(pareto_clu):
        pareto_clu_queue.append(np.copy(i))
        pareto_wss_queue.append(deepcopy(pareto_wss[e]))

    do_interchange = True
    do_trace = False

    lcounter = 0

    for e, initial in enumerate(pareto_clu_queue):

        lcounter += 1
        print("Investigating neighborhood of clustering no. %s" % e)
        #
        print("Starting from a new clustering.")

        new_clu, new_wss, flag_improved_intchg = improve_initial_workhorse(initial, do_interchange, do_trace)

        chash, chashtpl = make_clu_hash(new_clu)

        if flag_improved_intchg and chash not in chashset2chashtpl.keys():

            chash, chashtpl = make_clu_hash(new_clu)
            chashset2chashtpl[chash] = chashtpl

            dominated, flag_pareto = pf.eval_paretof(pareto_wss, new_wss)

            for dom_idx in dominated:
                dhash = make_clu_hash_set(pareto_clu[dom_idx])
                del chashset2chashtpl[dhash]

                del pareto_clu[dom_idx]
                del pareto_wss[dom_idx]

            pareto_clu.append(np.copy(new_clu))
            pareto_wss.append(deepcopy(tuple(new_wss)))
            #
            if flag_pareto:
                print("IMPROVEMENT made in MOVE/INTERCHANGE.")
            else:
                print("Improvement made in MOVE/INTERCHANGE but no impact on Pareto front.")
        else:
            print("No improvement in MOVE/INTERCHANGE transformation.")

    print("Number of all Pareto clusterings (candidates): %s" % len(pareto_clu))

    return pareto_clu, pareto_wss
