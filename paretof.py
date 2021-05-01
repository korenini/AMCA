import itertools
import settings as pset


def paretof(clu, wss):
    accepted = [True for i in range(len(clu))]

    for e, p in enumerate(wss):
        if accepted[e]:
            for f, r in enumerate(wss):
                if e == f:
                    continue
                else:
                    pass
                if accepted[f]:
                    flag_pareto = False
                    for pcri, rcri in zip(p, r):
                        if (rcri+pset.epsilon) < pcri:
                            flag_pareto = True
                            break
                        else:
                            pass
                    if not flag_pareto:
                        accepted[f] = False
                else:
                    pass
        else:
            pass

    return list(itertools.compress(clu, accepted)), list(itertools.compress(wss, accepted))


def paretof_add(clu, wss, clu_candidates, wss_candidates):
    cutoff = len(wss_candidates) * -1
    clu.extend(clu_candidates)
    wss.extend(wss_candidates)
    accepted = [True for i in range(len(clu))]
    for e, p in enumerate(wss):
        if accepted[e]:
            for f, r in enumerate(wss):
                if e == f:
                    continue
                else:
                    pass
                if accepted[f]:
                    flag_pareto = False
                    for pcri, rcri in zip(p, r):
                        if (rcri+pset.epsilon) < pcri:
                            flag_pareto = True
                            break
                        else:
                            pass
                    if not flag_pareto:
                        accepted[f] = False
                else:
                    pass
        else:
            pass
    num_new = sum(accepted[cutoff:]*1)

    return list(itertools.compress(clu, accepted)), list(itertools.compress(wss, accepted)), num_new


def paretof_trace_add(clu, wss, clu_candidates, wss_candidates, coclu_ft):
    cutoff = len(wss_candidates) * -1
    clu.extend(clu_candidates)
    wss.extend(wss_candidates)
    accepted = [True for i in range(len(clu))]
    #
    for e, p in enumerate(wss):
        if accepted[e]:
            for f, r in enumerate(wss):
                if e == f:
                    continue
                else:
                    pass
                if accepted[f]:
                    flag_pareto = False
                    for pcri, rcri in zip(p, r):
                        if (rcri+pset.epsilon) < pcri:
                            flag_pareto = True
                            break
                        else:
                            pass
                    if not flag_pareto:
                        accepted[f] = False
                else:
                    pass
        else:
            pass
    num_new = sum(accepted[cutoff:]*1)

    return list(itertools.compress(clu, accepted)), list(itertools.compress(wss, accepted)), list(itertools.compress(coclu_ft, accepted)), num_new


def eval_clust_nondominance(pareto_wss, cur_wss):
    flag_nondominated = True
    for p in pareto_wss:
        #
        if all([(pcri+pset.epsilon) < ccri for pcri, ccri in zip(p, cur_wss)]):
            flag_nondominated = False
            break
    return flag_nondominated


def eval_tip(paretofront_wss, move_wss):
    flag_pareto = False
    flag_borderline = False
    dominated = []

    for e, p in enumerate(paretofront_wss):
        if flag_pareto == False:
            if all([ccri-pset.epsilon > pcri for pcri, ccri in zip(p, move_wss)]):
                break
        if all([ccri+pset.epsilon <= pcri for pcri, ccri in zip(p, move_wss)]):
            if any([ccri+pset.epsilon < pcri for pcri, ccri in zip(p, move_wss)]):
                if not flag_pareto:
                    flag_pareto = True
                dominated.append(e)
            else:
                pass
    if dominated:
        dominated.sort(reverse=True)

    return dominated, flag_pareto, flag_borderline


def eval_paretof(paretofront_wss, move_wss):
    flag_pareto = True
    dominated = []
    for e, p in enumerate(paretofront_wss):
        if any([ccri + pset.epsilon < pcri for pcri, ccri in zip(p, move_wss)]):
            if all([ccri+pset.epsilon <= pcri for pcri, ccri in zip(p, move_wss)]):
                if any([ccri+pset.epsilon < pcri for pcri, ccri in zip(p, move_wss)]):
                    dominated.append(e)
            else:
                pass
        else:
            flag_pareto = False
            break
    if dominated:
        dominated.sort(reverse=True)

    return dominated, flag_pareto


def eval_paretof_quick(paretofront_wss, move_wss_sum):
    flag_pareto = False
    dominated = []

    for e, p in enumerate(paretofront_wss):
        #
        if flag_pareto == False:
            if all([ccri-pset.epsilon > pcri for pcri, ccri in zip(p, move_wss_sum)]):
                flag_nondominated = False
                break
        if all([ccri+pset.epsilon <= pcri for pcri, ccri in zip(p, move_wss_sum)]):
            if any([ccri+pset.epsilon < pcri for pcri, ccri in zip(p, move_wss_sum)]):
                if not flag_pareto:
                    flag_pareto = True
                dominated.append(e)
            else:
                pass
    if dominated:
        dominated.sort(reverse=True)

    return dominated, flag_pareto


def remove_dominated(pareto_clu, pareto_wss, cur_clu, cur_wss):
    accepted = [True for i in range(len(pareto_wss))]
    improved = False
    for e, p in enumerate(pareto_wss):
        if all([(cwcri+pset.epsilon) <= pcri for pcri, cwcri  in zip(p, cur_wss)]):
            if any([(cwcri+pset.epsilon) < pcri for pcri, cwcri  in zip(p, cur_wss)]):
                accepted[e] = bool(0)
                if not improved:
                    improved = True
    if improved:
        pareto_clu_new = list(itertools.compress(pareto_clu, accepted))
        pareto_wss_new = list(itertools.compress(pareto_wss, accepted))
        pareto_clu_new.append(cur_clu)
        pareto_wss_new.append(cur_wss)
    else:
        pareto_clu_new = pareto_clu
        pareto_wss_new = pareto_wss
    return pareto_clu_new, pareto_wss_new, improved



