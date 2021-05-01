import numpy as np
import get_data as gd
import settings as pset
import paretof as pf


def allecswss_bvect(aclu):
    cwss = []
    for i in gd.mcrits:
        sum_dist = 0
        for clm in range(aclu.shape[0]):
            members = aclu[clm]
            lenmembers = members.astype(int).sum()
            sum_dist += (gd.npdist[i][np.ix_(members, members)].sum()) / float(lenmembers * 2.0)
        cwss.append(sum_dist)
    return tuple(cwss)


def allecswss_bvect_k(aclu):
    cwss = np.zeros((gd.ncrits, aclu.shape[0]))
    for i in gd.mcrits:
        for clm in range(aclu.shape[0]):
            members = aclu[clm]
            lenmembers = members.astype(int).sum()
            cwss[i, clm] = (gd.npdist[i][np.ix_(members, members)].sum()) / float(lenmembers * 2.0)
    return cwss


def ecswss_bvect_active_k(aclu, crit, active_lst):
    adist = np.zeros(len(active_lst))
    for en, clm in enumerate(active_lst):
        members = aclu[clm]
        lenmembers = members.astype(int).sum()
        adist[en] = (gd.npdist[crit][np.ix_(members, members)].sum()) / float(lenmembers * 2.0)
    return adist


def allecswss_barithm(aclu, old_wss, mvun, clu_to):
    yclu = np.copy(aclu)
    yclu[:, mvun] = np.copy(clu_to)
    new_wss = allecswss_bvect(yclu)
    if np.all([(a + pset.epsilon) <= o for a, o in zip(new_wss, old_wss)]):
        if np.any([(a + pset.epsilon) <= o for a, o in zip(new_wss, old_wss)]):
            return new_wss, True
    return new_wss, False


def allecswss_barithm_intchg(aclu, old_wss, mvun, intchgun, clu_from, clu_to):
    yclu = np.copy(aclu)
    yclu[:, mvun] = np.copy(clu_to)
    yclu[:, intchgun] = np.copy(clu_from)
    new_wss = allecswss_bvect(yclu)
    if np.all([(a + pset.epsilon) <= o for a, o in zip(new_wss, old_wss)]):
        if np.any([(a + pset.epsilon) <= o for a, o in zip(new_wss, old_wss)]):
            return new_wss, True
    return new_wss, False


def allecswss_bpareto_mv(aclu, cur_wss, mvun, clu_to, paretofront_wss):
    yclu = np.copy(aclu)
    yclu[:, mvun] = np.copy(clu_to)
    new_wss = allecswss_bvect(yclu)

    if any([x+pset.epsilon < y for x, y in zip(new_wss, cur_wss)]):
        if pset.search_mode == "thorough":
            dominated, flag_pareto = pf.eval_paretof(paretofront_wss, new_wss)
        else:
            dominated, flag_pareto = pf.eval_paretof_quick(paretofront_wss, new_wss)

        return dominated, new_wss, flag_pareto
    else:
        return [], new_wss, False


def allecswss_bpareto_itchg(aclu, cur_wss, mvun, intchgun, clu_from, clu_to, paretofront_wss):
    yclu = np.copy(aclu)
    yclu[:, mvun] = np.copy(clu_to)
    yclu[:, intchgun] = np.copy(clu_from)
    new_wss = allecswss_bvect(yclu)

    if np.any([x+pset.epsilon < y for x,y in zip(new_wss, cur_wss)]):
        if pset.search_mode == "thorough":
            dominated, flag_pareto = pf.eval_paretof(paretofront_wss, new_wss)
        else:
            dominated, flag_pareto = pf.eval_paretof_quick(paretofront_wss, new_wss)

        return dominated, new_wss, flag_pareto
    else:
        return [], new_wss, False


def ecswss_bvect_scrit(aclu, crit):
    sum_dist = 0
    for clm in range(aclu.shape[0]):
        members = aclu[clm]
        lenmembers = members.astype(int).sum()
        sum_dist += (gd.npdist[crit][np.ix_(members, members)].sum()) / float(lenmembers * 2.0)
    return sum_dist


def ecswss_bvect_scrit_imprv(aclu, cur_crit_wss, crit, mvun, clu_to):
    yclu = np.copy(aclu)
    yclu[:, mvun] = np.copy(clu_to)
    sum_dist = 0
    for clm in range(aclu.shape[0]):
        members = yclu[clm]
        lenmembers = members.astype(int).sum()
        sum_dist += (gd.npdist[crit][np.ix_(members, members)].sum()) / float(lenmembers * 2.0)
    if sum_dist + pset.epsilon < cur_crit_wss:
        return sum_dist, True
    else:
        return sum_dist, False
