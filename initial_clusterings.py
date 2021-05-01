import numpy as np
import itertools
import random
import get_data as gd
import kmeans_trace as kmtrace
import paretof as pf
import settings as pset
from do_logging import logme, logme_name


def find_initial():

    clu = []
    wss = []

    coclu_ft = []

    per_crits = list(itertools.permutations(range(gd.ncrits), 2))

    crit_optimal_clu = [0 for l in range(gd.ncrits)]
    crit_optimal_wss = [0 for l in range(gd.ncrits)]

    print("")
    print("Searching for criterion optimal clusterings.")
    print("")

    for crit in range(gd.ncrits):
        crit_optimal_clu[crit], crit_optimal_wss[crit] = kmtrace.lkmbv(crit, kmtrace.randominit(), True)
        clu.append(np.copy(crit_optimal_clu[crit]))
        wss.append(crit_optimal_wss[crit])
        coclu_ft.append((crit, crit))
        print("Inspecting criterion %s. WSS: %s" % (crit+1,  crit_optimal_wss[crit]))

    print("")
    print("")
    print("Searching for best initial clusterings using tracing.")
    print("")

    counter = 0

    # Keeps track of the number of improvements for each criterion achieved in the loop
    len_new_container = [1 for x in range(len(per_crits))]

    while any(len_new_container) and counter < pset.iter_kmeans_max:

        counter += 1
        for en, ft in enumerate(per_crits):
            f, t = ft
            print("")
            print("From criterion %s optimal to criterion %s optimal." % (f+1, t+1))
            trace_clu, trace_wss = kmtrace.lkmbv(f, crit_optimal_clu[t], False)
            trace_coclu_ft = [(f, t)] * len(trace_wss)
            coclu_ft = coclu_ft + trace_coclu_ft
            clu, wss, coclu_ft, len_new = pf.paretof_trace_add(clu, wss, trace_clu, trace_wss, coclu_ft)
            len_new_container[en] = len_new


    # Report on the number of new clusterings in each 'per_crits' path
    print("F: from clustering optimized to criterion X")
    print("T: to clustering optimized to criterion Y")
    print("Number of iterations: %s" % counter)

    # From this module
    for e, i in enumerate(per_crits):
        f, t = i
        print("F:%s, T:%s, new clusterings: %s" % (f+1, t+1, len_new_container[e]))

    repetition_counter = 0
    clu_trcq = []

    for e, i in enumerate(clu):
        clu_trcq.append(np.copy(i))

    if gd.ncrits > 2:
        while repetition_counter < pset.iter_cross_max:
            randclu = random.randint(0, (len(coclu_ft)-1))
            f, t = coclu_ft[randclu]
            if f == t:
                continue
            repetition_counter += 1
            print("")
            print("Starting repetition %s" % (repetition_counter))
            #
            to_crits = [x for x in range(gd.ncrits) if x not in coclu_ft[randclu]]
            for to_crit in to_crits:
                trace_clu, trace_wss = kmtrace.lkmbv(to_crit, clu_trcq[randclu], False)
                clu, wss, len_new = pf.paretof_add(clu, wss, trace_clu, trace_wss)
                print("Number of new Pareto clusterings: %s" % len_new)

    print("TRACING: Number of all Pareto clusterings (candidates): %s" % len(clu))

    return clu, wss
