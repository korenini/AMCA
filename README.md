# AMCCA – alternative multicriteria clustering algorithm


Traditionally, clustering problems have been addressed by optimizing a criterion function over a set of variables that described a facet of a phenomenon in some meaningful way. Most often, algorithms belonging to the *k*-means clustering family have been used in the analysis. This way of cluster analysis can be referred to as a single–criterion approach. There are also many clustering problems that cannot be adequately addressed using the mentioned traditional approach. The research problem may consist of investigated units on which several sets of variables have been measured, each set describing a separate facet of the phenomenon. The measured variables can not be joined together into a single dataset.

For example, countries may be partitioned into clusters based on their ability to face the Covid-19 problem. This ability may be assessed from different angles, such as: the current state of the pandemic, measured as the number of infected people, the number of dead people due to infection (*criterion 1*), the current state of the healthcare system, measured as the number of doctors, hospital beds and ventilators available for every thousand inhabitants (*criterion 2*), the economic situation in the country, measured as GDP and credit ratings (*criterion 3*)...  
It is not to be expected that clusterings according to every single criterion are identical. Therefore clustering of units based on only a single criterion insufficiently describes the phenomenon. On the other hand, in multicriteria clustering, all criteria are optimized simultaneously. Most often, this does not result in a single clustering but in a set of clusterings which represent possible trade-offs to each single-criterion optimal clustering. From this set of clusterings, a clustering may be selected, which represents the best trade-off from the perspective of the research problem at hand.

Various, but not many, algorithms aimed at solving the multicriteria clustering problems exist. Most of them are designed to solve a specific problem bound to a research domain and not fit to be used elsewhere.

Direct multicriteria clustering algorithm [PDF](http://vlado.fmf.uni-lj.si/vlado/papers/multclu.pdf) represents a more general approach to solving multicriteria clustering problem. AMCCA is an extension of direct multicriteria clustering algorithm and aims to obtain the same results but (much) faster.  
Please keep in mind that, in general, the clustering problem is NP-hard, and multicriteria clustering problem is even harder. It is recommended that anyone wishing to familiarize themselves with this software uses small datasets. A small dataset (n=45) consisting of 6 variables, 3 criteria, can be found in the data folder.

## Code execution
The current implementation of AMCDA takes several datasets as input. Each dataset must be written to a separate file and specified in ``settings.py``. All other settings must be specified in ``settings.py``, which also contain sensible default values.
Python3 and NumPy are needed to run the software, optionally Matplotlib and Asciidoc, to get a nicer output.
After the input files and the desired number of clusters are specified in settings, the software can be executed.

```python
python3 pareto.py
```

## Sample outputs

### Cluster membership, *k*=3.
![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/clu_members.png "Cluster membership for each clustering.")

### The values of the criterion function.
![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/criterion_funct.png "The values of the criterion function.")

### Pareto clusterings, 2 criteria.
Initial clusterings are represented as triangles, and final clusterings are circles.

![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/paretofront-2crit.png "Pareto clusterings, 2 criteria.")

### Pareto clusterings, 3 criteria.
Initial clusterings are represented as triangles, and final clusterings are circles.

![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/paretofront-3crit-angle1.png "Pareto clusterings, 2 criteria.")
![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/paretofront-3crit-angle2.png "Pareto clusterings, 2 criteria.")
![alt text](https://github.com/korenini/AMCCA/blob/main/docs/_images/paretofront-3crit-angle3.png "Pareto clusterings, 2 criteria.")


