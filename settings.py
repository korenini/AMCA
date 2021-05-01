
## Data
#data_files = ["data/crit1.csv", "data/crit2.csv", "data/crit3.csv"]
data_files = ["data/crit1.csv", "data/crit2.csv"]

# Delimiter used in data files, e.g. ',' in case of csv.
delimiter = ','


# Analysis
# Number of clusters
k = 3

# Smallest permitted cluster size (spcs) used in single criterion clustering
spcs = 3

# Smallest permitted cluster size (mspcs) used in mulicriteria criterion clustering
mspcs = 3

# Hash/map function ['md5', 'sha256', 'base64']
map_function = 'base64'

# Maximum number of k-means iterations when searching for initial clusterings
iter_kmeans_max = 20

# Repetitions of cross tracing (random selected to crit optimal)
iter_cross_max = 20

# Random order of units in move & interchange transformation attempts; [ True | False ]
doshuffle = True

# Search mode [ quick | thorough ]
search_mode = "quick"


# Enable / disable logging [ True | False ]
logging = False
#logging = True



## Output
# Save data and metadata to an output file for later processing [ True | False ]
do_export = True

# Output filename (without extension)
 output_name = "pareto"

# Precision used in output (doesn't affect precision in calculations)
output_precision = 4

# Report output format [asciidoc | plain ascii].
# Asciidoc, http://asciidoc.org, is preferred.
output_format = "asciidoc"

# Asciidoc doctypes; default is html 
ad_doctype = "html"

## Plot
# Plot initial and Pareto clusterings [ True | False ]
do_plot = True


# Allowed number of iterations in tracing
trace_iter_max = 50


# Logging
# loglevel 0: there is no logging at all.
# loglevel 1: basic k-means and multicriterial logging
# loglevel 2: extended logging
# loglevel 3: elaborate logging which includes all local pareto clusterings
loglevel = 0

logfile = "logfile.txt"


epsilon = 1e-11


# Debugging
debug = False

# Log clustering procedure
log = True

