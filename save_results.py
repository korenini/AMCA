import numpy as np
import pickle
import datetime
import subprocess
from os import linesep
import compact_data as cdat
import check_deps as adc


from settings import output_name, output_format, ad_doctype, output_precision
from settings import k, data_files



def format_and_dump(paretofront_clu, paretofront_wss):
    ## Name files containing str results
    output_clu_str = output_name + "_clu.txt"
    output_wss_str = output_name + "_wss.txt"
    output_data = output_name + ".pickle"

    ## Initialize files with headers
    if output_format == "asciidoc":
        header_wss = "= The value of the criterion function for each criterion" + linesep
        header_wss += "Clustering no.,\t"
        header_wss += ",\t".join(["Crit."+str(crt+1) for crt in range(len(data_files))]) + linesep
        header_wss += linesep + "----" + linesep
        footer_wss = linesep + "----" + linesep
        header_clu = "= Cluster membership for each Pareto clustering" + linesep
        header_clu += linesep + "----" + linesep
        footer_clu = linesep + "----" + linesep
    else:
        header_wss = str("=" * 80) + linesep
        header_wss += "The value of the criterion function for each criterion" + linesep
        header_wss += str("=" * 80) + linesep
        header_wss += "Clustering no.,\t"
        header_wss += ",\t".join(["Crit."+str(crt+1) for crt in range(len(data_files))]) + linesep
        header_wss += str("-" * 80) + linesep
        header_clu = str("=" * 80) + linesep
        header_clu += "Cluster membership for each Pareto clustering"
        header_clu += linesep
        header_clu += str("=" * 80) + linesep

    # Write header to file
    with open(output_clu_str, 'w') as clu_str:
        clu_str.write(header_clu)

    with open(output_wss_str, 'w') as wss_str:
        wss_str.write(header_wss)

    ## Sort results

    sortarray = np.empty((len(data_files), len(paretofront_clu)))

    for e, w in enumerate(paretofront_wss):
        sortarray[:,e] = w[::-1]

    sortkey = np.lexsort(sortarray)

    ## Write data to file

    # Clusterings
    for e,sk in enumerate(sortkey):
        clu_no = "Clustering no. %s %s" % (e+1, linesep)
        clu_lst0 = cdat.bitvect2clu(paretofront_clu[sk]).tolist()
        clu_lst1 = [clumem+1 for clumem in clu_lst0]
        clu_str = "("
        for e,c in enumerate(clu_lst1):
            clu_lencounter = len(clu_str.split(linesep)[-1])
            if clu_lencounter % 80 < 70:
                clu_str += str(c)
                if e != len(clu_lst1) -1:
                    clu_str += ", "
                else:
                    clu_str += ")" + linesep
            else:
                clu_str += linesep + " "
                clu_str += str(c)
                if e != len(clu_lst1) -1:
                    clu_str += ", "
                else:
                    clu_str += ")" + linesep
        clu_str += linesep
        clu_str_out = clu_no + clu_str
        with open(output_clu_str, 'a') as clu_str:
            clu_str.write(clu_str_out)

    # Criterion function
    precision_str = "%." + str(output_precision) + "f"

    for e,sk in enumerate(sortkey):
        wss_rounded = np.around(paretofront_wss[sk], output_precision)
        wss_formated = [precision_str % cf for cf in wss_rounded]
        outstring_wss = "%s)\t" % str(e+1)
        outstring_wss += ",\t".join([w for w in wss_formated])
        outstring_wss += linesep
        # append to output file
        with open(output_wss_str, 'a') as wss_str:
            wss_str.write(outstring_wss)

    # Asciidoc report footers & execution
    if output_format == "asciidoc" and adc.adc():
        with open(output_clu_str, 'a') as clu_str:
            clu_str.write(footer_clu)
        with open(output_wss_str, 'a') as wss_str:
            wss_str.write(footer_wss)
        # Execute asciidoc files
        subprocess.run(["asciidoc", "-b", ad_doctype, output_clu_str])
        subprocess.run(["asciidoc", "-b", ad_doctype, output_wss_str])


    # Dump Pareto clusterings data for later retrieval & analysis
    pf_clu_sorted_idx = []
    pf_wss_sorted_idx = []

    for sk in sortkey:
        clu_idx = cdat.bitvect2clu(paretofront_clu[sk]).tolist()
        wss_idx = paretofront_wss[sk]
        pf_clu_sorted_idx.append(clu_idx)
        pf_wss_sorted_idx.append(wss_idx)


    data_dict = {}
    data_dict["created"] = datetime.datetime.now().strftime("%Y-%m-%d %T")
    data_dict["data_files"] = data_files
    data_dict["ncrits"] = len(data_files)
    data_dict["ncases"] = len(data_files[0])
    data_dict["k"] = k
    data_dict["clu"] = pf_clu_sorted_idx
    data_dict["wss"] = pf_wss_sorted_idx

    with open(output_data, 'wb') as out_clu_dat:
        out_clu_dat.write(pickle.dumps(data_dict))


# INSTRUCTIONS: Restore the pickled output file in a new Python session
"""
import pickle
from settings import output_name

output_data = output_name + ".pickle"

pf_clu_restored = []
pf_wss_restored = []

try:
    with open(output_data, 'rb') as data_dumped:
        data_dict = pickle.load(data_dumped)
    print("Data loaded.")
except FileNotFoundError:
    print("File %s not found" % output_data)

paretofront_clu = data_dict['clu']
paretofront_wss = data_dict['wss']
"""




