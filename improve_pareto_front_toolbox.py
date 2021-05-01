import numpy as np
import settings as pset
import get_data as gd


def check_clu_mv_transformations(rmtah_chash, aclu_size):
    if np.sum(rmtah_chash.astype(int)) == aclu_size:
        return True
    else:
        return False


def check_clu_mvi_transformations(ritah_chash, aclu_size):
    if np.sum(ritah_chash.astype(int)) == aclu_size:
        return True
    else:
        return False


def check_unit_mv_transformations(rmtah_chash, mvun):
    if np.sum(rmtah_chash[:, mvun].astype(int)) == pset.k:
        return True
    else:
        return False


def check_unit_mvi_transformations(ritah_chash, mvun):
    if np.sum(ritah_chash[:, mvun].astype(int)) == pset.k:
        return True
    else:
        return False


def mask_inactive_mv(aclu, clu_from, clu_to, mvun):
    if pset.k <= 3:
        inactive_mask_tran = np.copy(aclu)
        inactive_mask_val = np.full(aclu.shape, False)
        return inactive_mask_tran, inactive_mask_val

    inactive_mask_tran = np.full(aclu.shape, True)
    inactive_mask_val = np.full(aclu.shape, True)
    inactive_mask_tran[clu_from] = False
    inactive_mask_tran[clu_to] = False
    inactive_mask_val[clu_from] = False
    inactive_mask_val[clu_to] = False
    active_clusters = np.logical_or(clu_from, clu_to)
    active_units_idx = [x for x in range(gd.nunits) if np.any(aclu[active_clusters,x])]
    for aun in active_units_idx:
        inactive_mask_tran[:, aun] = np.copy(aclu[:, aun])
        inactive_mask_val[:, aun] = False

    return inactive_mask_tran, inactive_mask_val


def mask_inactive_intchg(aclu, clu_from, clu_to, mvun, intchgun):
    if pset.k <= 3:
        inactive_mask = np.copy(aclu)
        return inactive_mask

    inactive_mask = np.full(aclu.shape, True)
    inactive_mask[clu_from] = False
    inactive_mask[clu_to] = False
    active_clusters = np.logical_or(clu_from, clu_to)
    active_units_idx = [x for x in range(gd.nunits) if np.any(aclu[active_clusters,x])]
    for aun in active_units_idx:
        inactive_mask[:,aun] = np.copy(aclu[:,aun])
    inactive_mask[clu_from, mvun] = True
    inactive_mask[clu_to, intchgun] = True

    return inactive_mask


def mask_inactive(aclu, clu_from, clu_to):
    if pset.k <= 3:
        inactive_mask_tran = np.copy(aclu)
        inactive_mask_val = np.full(aclu.shape, False)
        return inactive_mask_tran, inactive_mask_val

    inactive_mask_tran = np.full(aclu.shape, True)
    inactive_mask_val = np.full(aclu.shape, True)
    inactive_mask_tran[clu_from] = False
    inactive_mask_tran[clu_to] = False
    inactive_mask_val[clu_from] = False
    inactive_mask_val[clu_to] = False
    active_clusters = np.logical_or(clu_from, clu_to)
    active_units_idx = [x for x in range(gd.nunits) if np.any(aclu[active_clusters, x])]
    for aun in active_units_idx:
        inactive_mask_tran[:,aun] = np.copy(aclu[:,aun])
        inactive_mask_val[:,aun] = False

    return inactive_mask_tran, inactive_mask_val


def get_mvbackunits(ritahset_chash, cltidx, mvun, ac):
    cell_content = ritahset_chash[cltidx, mvun]

    if not cell_content:
        mvback_candidates_idx = [x for x in range(gd.nunits) if ac[x] and x is not mvun]
        return mvback_candidates_idx, set()
    elif isinstance(cell_content, set): # Cell contains a set which contains units already moved back to clu_from
        movedbackset = cell_content
        mvback_candidates_idx = [x for x in range(gd.nunits) if ac[x] and x not in movedbackset and x is not mvun]
        return mvback_candidates_idx, cell_content
    else:
        mvback_candidates_idx = [x for x in range(gd.nunits) if ac[x] and x is not mvun]
        return mvback_candidates_idx, set()
