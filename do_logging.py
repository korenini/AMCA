import copy
from settings import logfile


def logclean():
    with open(logfile, "w") as f:
        f.write("")


def logme(mystring):
    with open(logfile, "a") as f:
        f.write("%s \n" % str(mystring))


def logme_name(mystring, filename):
    with open(filename, "a") as f:
        f.write("%s \n" % str(mystring))


def log_strings(mystring):
    with open(logfile, "a") as f:
        f.write("%s \n" % mystring)


def log_crits(crits):
    p1, p2 = crits
    with open(logfile, "a") as f:
        f.write("%.2f, %.2f " % (p1, p2))
            
            
def log_list(cl_members):
    tmp_list = []
    for i in cl_members:
        tmp_list.append(sorted(copy.copy(i)))
    with open(logfile, "a") as f:
        f.write("%s \n" % tmp_list)  
