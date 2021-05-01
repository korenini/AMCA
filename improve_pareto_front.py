import numpy as np
import random
from copy import copy
from copy import deepcopy
import settings as pset
import get_data as gd
import compact_data as cdat
import mapping
import paretof as pf
import ecs_smpl as ecs
from clu_bitvect_precomp import all_bit_clu, idx2vect, vect2idx
import improve_pareto_front_toolbox as pftbox

pset.do_interchange = False

move_endpoints = []
while_loopcounter = 0


def improve_pareto(aclui, do_interchange, do_trace):

    global dominated_enc
    global paretofront_clu, paretofront_wss
    global chash, chashtpl, chashset2chashtpl
    global rmtah, ritah, ritahset

    aclu = np.copy(aclui)

    units = list(range(gd.nunits))
    clulst = list(range(pset.k))

    flag_continue = True
    flag_improvement = False

    num_clu_new = 0
    num_clu_dominated = 0

    trace_clu = []
    trace_wss = []

    cur_wss = ecs.allecswss_bvect(aclu)

    rmtah_chash = rmtah[chash]
    try:
        ritah_chash = ritah[chash]
        ritahset_chash = ritahset[chash]
    except KeyError:
        ritah[chash] = np.copy(aclu)
        ritahset[chash] = np.full(aclu.shape, set(), dtype=object)
        ritah_chash = ritah[chash]
        ritahset_chash = ritahset[chash]

    while_loopcounter = 0

    while flag_continue:

        csd_val = cdat.cluster_size_all(aclu)
        csd_tf = csd_val <= (pset.mspcs)
        flag_clustersize = any(csd_tf)

        flag_continue = False
        flag_break = False

        if pftbox.check_clu_mv_transformations(rmtah_chash, aclu.size):
            if do_interchange:
                if pftbox.check_clu_mvi_transformations(ritah_chash, aclu.size):
                    print("All move-interchange transformations for clustering already tested.")
                    flag_continue = False
                    continue
            else:
                print("All move transformations for clustering already tested.")
                flag_continue = False
                continue

        if pset.doshuffle:
            random.shuffle(units)

        while_loopcounter += 1
        print("Iteration: %s" % while_loopcounter)

        for mvun in units:

            if flag_break:
                break

            if not do_interchange:
                if pftbox.check_unit_mv_transformations(rmtah_chash, mvun):
                    print("All move transformations tested for unit: %s." % mvun)
                    continue
            else:
                if pftbox.check_unit_mvi_transformations(ritah_chash, mvun):
                    print("All move-interchange transformations tested for unit: %s." % mvun)
                    continue
            csd_val = cdat.cluster_size_all(aclu)
            csd_tf = csd_val <= (pset.mspcs)
            flag_clustersize = any(csd_tf)
            #
            # MOVE
            #
            clu_from = np.copy(aclu[:, mvun])
            clu_from_idx = vect2idx[tuple(clu_from)]
            #
            if pset.doshuffle:
                random.shuffle(clulst)

            for clx in clulst:
                if flag_break:
                    break

                clu_to = np.copy(idx2vect[clx])
                cltidx = clx

                if all(clu_from == clu_to):
                    continue

                flag_move_allowed = bool(1)
                if flag_clustersize:
                    if any(clu_from & csd_tf):
                        if do_interchange:
                            flag_move_allowed = bool(0)
                        else:
                            print("Cluster too small. Skipping.")
                            continue
                # Move scenario
                if not do_interchange:
                    if rmtah_chash[cltidx, mvun]:
                        print("Move transformation already tested - unit: %s, cluster: %s" % (mvun, cltidx))
                        continue
                    else:
                        # print("Testing transformation - unit: %s to cluster: %s" % (mvun, cltidx))
                        dominated, move_wss, flag_improvement_mv = ecs.allecswss_bpareto_mv(aclu, cur_wss, mvun, clu_to,
                                                                                            paretofront_wss)
                        rmtah_chash[cltidx, mvun] = True
                    if not flag_move_allowed:
                        print("Cluster too small, transformation is not allowed")
                        continue
                # Move-interchange scenario
                else:
                    if rmtah_chash[cltidx, mvun] and ritah_chash[cltidx, mvun]:
                        print("All interchange transformations already tested - unit: %s, cluster: %s" % (mvun, cltidx))
                        continue
                    else:
                        # print("Testing transformation. unit: %s to cluster: %s" % (mvun, cltidx))
                        dominated, move_wss, flag_improvement_mv = ecs.allecswss_bpareto_mv(aclu, cur_wss, mvun, clu_to,
                                                                                            paretofront_wss)
                        rmtah_chash[cltidx, mvun] = True
                        rmtah[chash][cltidx, mvun] = True

                        if flag_improvement_mv and not flag_move_allowed:
                            print("MOVE IMPROVEMENT but cluster to small for transformation.")
                            flag_improvement_mv = False
                            dominated = []
                # (MOVE) IMPROVEMENT
                if flag_improvement_mv:
                    flag_continue = True
                    flag_improvement = True

                    num_clu_new += 1
                    num_clu_dominated += len(dominated)

                    rmtah[chash] = np.copy(rmtah_chash)
                    ritah[chash] = np.copy(ritah_chash)
                    ritahset[chash] = np.copy(ritahset_chash)

                    aclu[:, mvun] = np.copy(clu_to)

                    new_hash, new_hash_tpl = mapping.make_clu_hash(aclu)
                    if new_hash in rmtah.keys():
                        aclu[:, mvun] = np.copy(clu_from)
                        flag_continue = False
                        flag_improvement = False
                        continue
                    else:
                        chash, chashtpl = new_hash, new_hash_tpl

                    csd_val = cdat.cluster_size_all(aclu)

                    csd_tf = csd_val <= pset.mspcs
                    flag_clustersize = any(csd_tf)

                    cur_wss = deepcopy(move_wss)

                    chashset2chashtpl[chash] = chashtpl

                    inactive_mask_tran, inactive_mask_val = pftbox.mask_inactive(aclu, clu_from, clu_to)

                    rmtah_chash = np.where(inactive_mask_tran, rmtah_chash, False)
                    ritah_chash = np.where(inactive_mask_tran, ritah_chash, False)
                    ritahset_chash = np.where(inactive_mask_val, ritahset_chash, set())

                    for dom_idx in dominated:
                        denc = cdat.encode_bin_to_base64_clu(paretofront_clu[dom_idx])
                        dominated_enc.add(denc)
                        dhash, dhash_list = mapping.make_clu_hash(paretofront_clu[dom_idx])
                        if chash == dhash:
                            continue

                        # Delete from Pareto front
                        del paretofront_clu[dom_idx]
                        del paretofront_wss[dom_idx]

                    # Append aclu, cur_wss to Pareto front
                    paretofront_clu.append(np.copy(aclu))
                    paretofront_wss.append(copy(tuple(cur_wss)))
                    print("Improvement in MOVE. WSS: %s" % str(cur_wss))

                    # TRACING
                    if do_trace:
                        trace_clu.append(np.copy(aclu))
                        trace_wss.append(ecs.allecswss_bvect(aclu))

                    #print("Move improvement. unit: %s, cluster: %s" % (mvun, cltidx+1))
                    flag_break = True
                    break
                elif do_interchange:
                    mvback_candidates_idx, movedbackset = pftbox.get_mvbackunits(ritahset_chash, cltidx, mvun,
                                                                                 aclu[clu_to][0])

                    if pset.doshuffle:
                        random.shuffle(mvback_candidates_idx)

                    for intchgun in mvback_candidates_idx:

                        if flag_break:
                            break

                        dominated, interchange_wss, interchange_improvement = ecs.allecswss_bpareto_itchg(aclu,
                                                                              cur_wss, mvun, intchgun, clu_from,
                                                                              clu_to, paretofront_wss)
                        movedbackset.add(intchgun)
                        if all([bool(1) for mci in mvback_candidates_idx if mci in movedbackset]):
                            ritah_chash[cltidx, mvun] = True
                            ritahset_chash[cltidx, mvun] = set()
                            ritah[chash] = np.copy(ritah_chash)
                            ritahset[chash] = np.copy(ritahset_chash)

                        if interchange_improvement:
                            flag_continue = True
                            flag_improvement = True

                            num_clu_new += 1
                            num_clu_dominated += len(dominated)

                            rmtah[chash] = np.copy(rmtah_chash)
                            ritah[chash] = np.copy(ritah_chash)
                            ritahset[chash] = np.copy(ritahset_chash)

                            aclu[:, mvun] = np.copy(clu_to)
                            aclu[:, intchgun] = np.copy(clu_from)

                            check_hash = mapping.make_clu_hash_set(aclu)

                            try:
                                blind_check = rmtah[check_hash]
                                aclu[:, mvun] = np.copy(clu_from)
                                aclu[:, intchgun] = np.copy(clu_to)
                                flag_continue = False
                                flag_improvement = False
                                continue
                            except KeyError:
                                pass

                            cur_wss = interchange_wss
                            print("Improvement in INTERCHANGE. WSS: %s" % str(cur_wss))

                            # Remove dominated
                            for dom_idx in dominated:
                                denc = cdat.encode_bin_to_base64_clu(paretofront_clu[dom_idx])
                                dominated_enc.add(denc)

                                dhash, dhash_list = mapping.make_clu_hash(paretofront_clu[dom_idx])
                                if chash == dhash:
                                    continue
                                # Delete from pareto
                                del paretofront_clu[dom_idx]
                                del paretofront_wss[dom_idx]

                            # Append aclu, cur_wss to Pareto front
                            paretofront_clu.append(np.copy(aclu))
                            paretofront_wss.append(copy(tuple(cur_wss)))

                            chash, chashtpl = mapping.make_clu_hash(aclu)

                            chashset2chashtpl[chash] = chashtpl

                            inactive_mask_tran, inactive_mask_val = pftbox.mask_inactive(aclu, clu_from, clu_to)
                            rmtah_chash = np.where(inactive_mask_tran, rmtah_chash, False)

                            ritah_chash = np.where(inactive_mask_tran, ritah_chash, False)
                            ritahset_chash = np.where(inactive_mask_val, ritahset_chash, set())

                            rmtah[chash] = np.copy(rmtah_chash)
                            ritah[chash] = np.copy(ritah_chash)
                            ritahset[chash] = np.copy(ritahset_chash)

                            if do_trace:
                                trace_clu.append(aclu)
                                trace_wss.append(cur_wss_sum)
                            flag_break = True
                            break

        rmtah[chash] = np.copy(rmtah_chash)
        ritah[chash] = np.copy(ritah_chash)
        ritahset[chash] = np.copy(ritahset_chash)

    return aclu, cur_wss, flag_improvement, num_clu_new, num_clu_dominated


def translate_alignment(chashtpl, chashtpl_new):
    available = [bool(1) for x in range(pset.k)]
    match_ft = [0 for x in range(pset.k)]

    for e1 in range(pset.k):
        for e2 in range(pset.k):
            if available[e2]:
                if chashtpl[e1] == chashtpl_new[e2]:
                    match_ft[e1] = e2
                    available[e2] == bool(0)

    return match_ft


def improve_pareto_mv(clu, wss, metapack):

    global dominated_enc
    global paretofront_clu, paretofront_wss
    global chash, chashtpl, chashset2chashtpl
    global rmtah, ritah, ritahset

    rmtah, ritah, ritahset, chashset2chashtpl, dominated_enc = metapack

    do_interchange = False
    do_trace = False

    clu_intchg_queue = []
    wss_intchg_queue = []

    for e, i in enumerate(clu):
        clu_intchg_queue.append(np.copy(i))
        wss_intchg_queue.append(deepcopy(wss[e]))

    paretofront_clu = []
    paretofront_wss = []

    for e, i in enumerate(clu):
        paretofront_clu.append(np.copy(i))
        paretofront_wss.append(deepcopy(wss[e]))

    lcounter = 0

    flag_newclusterings = True

    while flag_newclusterings:

        num_clu_new_all = 0
        num_clu_dominated_all = 0

        for e, initial in enumerate(clu_intchg_queue):

            lcounter += 1

            if pf.eval_clust_nondominance(paretofront_wss, wss_intchg_queue[e]):
                print("Starting from a new clustering: %s." % lcounter)

                chash, chashtpl = mapping.make_clu_hash(initial)

                rmtah_chash = rmtah[chash]

                if np.sum(rmtah_chash.astype(int)) == initial.size:
                    print("Prevented from entering the queue. All move transitions already tested.")
                    continue

                try:
                    chashtpl_known = chashset2chashtpl[chash]
                    if chashtpl == chashtpl_known:
                        pass
                    else:
                        talign = translate_alignment(chashtpl, chashtpl_known)
                        rmtah[chash] = np.copy(rmtah[chash][talign])
                        chashset2chashtpl[chash] = chashtpl
                except KeyError:
                    print("KeyError in processing initial @ fast move - this should not ever happen!!!")
                    rmtah[chash] = np.copy(initial)
                    chashset2chashtpl[chash] = chashtpl
                #
                new_clu_mv_tip, new_wss_k_mv_tip, moved, num_clu_new, num_clu_dominated = improve_pareto(initial,
                                                                                                         do_interchange,
                                                                                                         do_trace)
                if moved:
                    num_clu_new_all += num_clu_new
                    num_clu_dominated_all += num_clu_dominated
            else:
                print("Initial clustering is DOMINATED")

        print("Number of new clustertings: %s" % num_clu_new_all)
        print("Number of dominated clustertings: %s" % num_clu_dominated_all)

        clu_intchg_queue = []
        wss_intchg_queue = []

        for enum, ithc in enumerate(paretofront_clu):
            clu_intchg_queue.append(np.copy(ithc))
            wss_intchg_queue.append(np.copy(paretofront_wss[enum]))

        # Loop continuation logic
        if num_clu_new_all == 0:
            flag_newclusterings = False

    metapack = [rmtah, ritah, ritahset, chashset2chashtpl, dominated_enc]

    return paretofront_clu, paretofront_wss, metapack


def improve_pareto_intchg(clu, wss, metapack):

    global dominated_enc
    global paretofront_clu, paretofront_wss
    global chash, chashtpl, chashset2chashtpl
    global rmtah, ritah, ritahset

    rmtah, ritah, ritahset, chashset2chashtpl, dominated_enc = metapack

    do_interchange = True
    do_trace = False

    clu_intchg_queue = []
    wss_intchg_queue = []

    for e, i in enumerate(clu):
        clu_intchg_queue.append(np.copy(i))
        wss_intchg_queue.append(deepcopy(wss[e]))

    paretofront_clu = []
    paretofront_wss = []

    for e, i in enumerate(clu):
        paretofront_clu.append(np.copy(i))
        paretofront_wss.append(deepcopy(wss[e]))

    flag_newclusterings = True

    while flag_newclusterings:

        lcounter = 0

        num_clu_new_all = 0
        num_clu_dominated_all = 0

        for e, initial in enumerate(clu_intchg_queue):

            initial_enc = cdat.encode_bin_to_base64_clu(initial)
            if initial_enc in dominated_enc:
                print("The queued clustering no.: %s is now a local pareto. Skipping." % e)
                continue

            chash, chashtpl = mapping.make_clu_hash(initial)

            rmtah_chash = rmtah[chash]

            try:
                ritah_chash = ritah[chash]
                if pftbox.check_clu_mvi_transformations(ritah_chash, initial.size):
                    print("All move-interchange transformations already checked.")
                    continue
            except KeyError:
                ritah[chash] = np.copy(initial)
                ritahset[chash] = np.full(initial.shape, set(), dtype=object)

            try:
                chashtpl_known = chashset2chashtpl[chash]
                if chashtpl == chashtpl_known:
                    pass
                else:
                    talign = translate_alignment(chashtpl, chashtpl_known)
                    rmtah[chash] = np.copy(rmtah[chash][talign])
                    chashset2chashtpl[chash] = chashtpl
                    ritah[chash] = np.copy(ritah[chash][talign])
                    ritahset[chash] = np.copy(ritahset[chash][talign])
            except KeyError:
                print("KeyError in processin initial @ fast move - this should not ever happen!!!")
                rmtah[chash] = np.copy(rmtah_chash)
                chashset2chashtpl[chash] = chashtpl
                ritah[chash] = np.copy(initial)
                ritahset[chash] = np.full(aclu.shape, set(), dtype=object)

            lcounter += 1
            print("FAST MOVE-INTERCHANGE iteration: %s" % lcounter)

            if pf.eval_clust_nondominance(paretofront_wss, wss_intchg_queue[e]):
                print("Starting from a new clustering.")

                new_clu_mv_tip, new_wss_k_mv_tip, moved, num_clu_new, num_clu_dominated = improve_pareto(initial,
                                                                                          do_interchange, do_trace)
                if moved:
                    new_pareto_encoded = cdat.encode_bin_to_base64_clu(new_clu_mv_tip)

                    num_clu_new_all += num_clu_new
                    num_clu_dominated_all += num_clu_dominated
                else:
                    pass
            else:
                print("Initial clustering is DOMINATED")
        #
        #
        #
        print("Number of new clustertings: %s" % num_clu_new_all)
        print("Number of dominated clustertings: %s" % num_clu_dominated_all)

        clu_intchg_queue = []
        wss_intchg_queue = []

        for e, i in enumerate(paretofront_clu):
            clu_intchg_queue.append(np.copy(i))
            wss_intchg_queue.append(np.copy(paretofront_wss[e]))

        # Loop continuation logic
        if num_clu_new_all == 0:
            flag_newclusterings = False

    metapack = [rmtah, ritah, ritahset, chashset2chashtpl, dominated_enc]

    return paretofront_clu, paretofront_wss, metapack
