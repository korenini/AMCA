#!/usr/bin/env python3

""" AMCCA, is an alternative multicriteria clustering algorithm.
    Copyright (C) 2020  Bojan Korenini
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>."""


import sys
from copy import deepcopy
import settings as pset

import mapping
import get_data as gd
import initial_clusterings as initclu
import improve_initial as improveinit
import improve_pareto_front as ipf

import check_deps
import save_results
import plotme


__author__ = "Bojan Korenini"
__credits__ = ["", ]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bojan Korenini"
__email__ = "bojan@korenini.net"
__status__ = "Development"


def main():

    # Check consistency
    if pset.k < 2 or pset.k > gd.nunits:
        print("Number of clusters can not be %s." % pset.k)
        sys.exit(1)

    if pset.spcs < 1 or pset.mspcs < 1:
        print("Empty clusters are not allowed. Please correct spcs and mspcs values.")
        sys.exit(1)

    if pset.ad_doctype == "asciidoc" and not check_deps.adc():
        print("Can not make an asciidoc output without asciidoc installed.")
        sys.exit(1)

    if pset.do_plot and not check_deps.check_mpl():
        print("Can not plot without Matplotlib installed.")

    # Search for good initial clusterings
    clu, wss = initclu.find_initial()
    clu_init, wss_init = deepcopy(clu), deepcopy(wss)

    # Improve initial clusterings
    chashset2chashtpl = mapping.make_pf_fs2tpl_hashes(clu)
    clu, wss = improveinit.improve_initial_mv(clu, wss, chashset2chashtpl)
    clu, wss = improveinit.improve_initial_intchg(clu, wss, chashset2chashtpl)

    dominated_enc = set()
    rmtah, ritah, ritahset, chashset2chashtpl = mapping.make_metadata(clu)
    metapack = [rmtah, ritah, ritahset, chashset2chashtpl, dominated_enc]

    clu, wss, metapack = ipf.improve_pareto_mv(clu, wss, metapack)
    clu, wss, metapack = ipf.improve_pareto_intchg(clu, wss, metapack)

    # Save results
    save_results.format_and_dump(clu, wss)

    # Plot results
    if pset.do_plot:
        if gd.ncrits == 2:
            plotme.plotme2d(wss_init, wss)
        elif gd.ncrits == 3:
            plotme.plotme3d(wss_init, wss)
        else:
            print("Can only plot data up to three dimensions.")


if __name__ == "__main__":
    
    # Check Python version. Must be Python 3!
    if (sys.version_info > (3, 0)):
        pass
    else:
        print("Please use Python version 3 to run this program.")
        sys.exit(1)

    main()