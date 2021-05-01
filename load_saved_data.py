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

