# WGCNA_shiny

![](https://github.com/chrisclarkson/pics/blob/master/Poster.png)

# Introduction
This is work from my masters thesis. I had to make a website that took biological data as input and infer a network. I chose WGCNA because the output was interpretable. (For any information of WGCNA please see: https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/)

# Purpose
The code herein turns the main features of the WGCNA workflow into an interactive shiny page allowing for co-expression of cluster identification, network visualisation and correlation of eigengeme values to phenotypes.

You can visualise the network on the shiny page or have a more dynamic network visualisation output can be saved.

# Sample data
I have provided sample data in 'microbial_abundance.txt' which contains compositional microbiotic co-abundance data. I should mention that this data is not suited to WGCNA and the output can not be interpretted as biologically meaningful. However I am supplying this data as it will give a fully annotated example of the output that this code can produce.


# Installation

git clone https://github.com/chrisclarkson/WGCNA_shiny/
cd WGCNA_shiny/
R
runApp('./app.R')

# Methodology

# Note: Necessary file format
The initially uploaded file (containing gene-expression data or compositional microbiotic co-abundance data) must contain the genes/genera that are to be correlated and clustered in the columns and the samples must be specified in the rows. For labeling the nodes and to maximize the effectiveness of the displayed network the file must contain a header that specifies the column names, otherwise the nodes will be labeled independently.

In the website sidebar, there is a navigation-panel that allows the user to specify if the files to be uploaded are text files or CSV files. Users can also specify how the files are delimited.


# Step 1: Filtering

![](https://github.com/chrisclarkson/pics/blob/master/filtering.png)

The first step in the pipeline is to filter out samples and nodes that are not informative. The first round of HC is done to remove samples too far displaced from all others and would perturb the inference of clusters if included. The user will choose the threshold via an interactive slider entitled "Minimum cluster size".


# Step 2: Soft threshold selection
The user will choose the threshold via an interactive slider entitled "Choose softthreshold for scale-free topology".

The command ‘pickSoftThreshold’ helps the user to choose a soft-threshold for the network to reach a scale-free topology. In accordance with the equation: 

sij=cor|(xi,xj)|

weights the adjacency-matrix and hence penalizes weak correlations. The value of beta is essential for the network to reach a scale-free topology.
The recommended threshold is output to the page. Furthermore, in the event that the user has an intuition that beta value should be different than the recommended power- the R2 fit-value to scale-free topology is plotted for each power.

# Step 3: An adjacency-matrix is made via Pearson correlation.


After correlating the expression-values via Pearson-correlation, to assign weights to the interactions. This is implemented via the command: ‘adjacency’.


# Step 4: Cluster inference


The second round of HC allows for the detection of modules within the dataset. This is done with the command ‘TOMsimilarity’, which documents the number of shared neighbors among multiple nodes. Hence, instead of directly utilizing the pair-wise adjacencies from Pearson correlation, multi-node comparisons are used and stronger connections within the pair-wise correlation values are identified. Subsequently HC is performed to infer Topolgical Overlap and hence modules. The calculated similarity values are subtracted from 1 to convert them to dissimilarity measures. This is done to eliminate noise caused by spurious or erroneous adjacencies.
 

![](https://github.com/chrisclarkson/pics/blob/master/Workflow.png)

Part A- A filtration cutoff level is selected and marked via a red line. Any samples that fall above that are line are considered outliers and are omitted from the dataset.


Part B- A heatmap represents the adjacency matrix produced by HC of the filtered dataset where the columns are treated with Pearson correlation to produce an adjacency matrix.


Part C- A soft-threshold is selected via a plot that documents the R2 fit to scale-free topology. Above the soft threshold value for reaching scale-free topology can be selected in accordance with the above graph where the R2 fit-value is plotted against each individual threshold-value.


Part D- HC of the TOM dissimilarity values identifies modules and colours are assigned to each simulated module.


# Step 5: Network-display


Using the package ‘igraph’ the network is plotted and the cluster members are displayed. The module colors of WGCNA-assigned module clusters are displayed in the network.


# Node-size and Edge-thickness

These two attributes were governed both in igraph and networkD3 by ‘enrichment factors’ whereby a user can choose a value to multiply the logarithmic degree per node and edge-strength to enrich the node-size and edge-thickness consecutively. The user will choose the enrichment level via an interactive slider entitled "Choose edge-enrichment".


# Layout Algorithms

Both in igraph and networkD3, the node-placement/ layout algorithms are ‘force-directed’ this is to say that there is an automatic repulsion between each node and yet the edges that exist between them cause an attractive force. The stronger the edges between a pair of nodes, the more closely that pair is placed together and this promotes the clustering of closely interacting nodes.
However the selection is defaulted to ‘layout.auto’- a command that selects the optimal of three algorithm choices based on the number of nodes.


# Step 7: Interactive graph 


Using the package ‘networkD3’ an interactive network is made from the igraph data and stored in a HTML file, saved using the "Save interactive graph" button and the page will be saved to the working directory. 

The node labels are on permanent display and the user can zoom in and investigate particular nodes. Again node-sizes are in accordance with the number of edges that a node has. The opacity of the nodes is set such that both the node labels and the nodes themselves are simultaneously visible.


# Step 8: Assign relatedness to modules

To relate the assigned modules to one another, a final round of clustering is applied to the Eigengenes that summarise the expression data of each cluster and this gives a measure of how related the modules are amongst each other. This can be compared to both the module-identifying dendogram and the network where module colours that are positioned closely will be shown to have close relation in this dendogram.


![](https://github.com/chrisclarkson/pics/blob/master/graph_eigengene.png)

Part A: displays a portion of the HTML output of the networkD3 package (as mentioned above in “Step 7”, a legend is necessary to translate the independent networD3 coloring-scheme to that of WGCNA). 



Part B: shows the hierarchical clustering of the corresponding Module-Eigengene values to show which are the most closely related. The red arrows illustrate how (as mentioned in “Step 8 “) the node-placement of networkD3 concurs with the clustering of the Eigengene-values in showing that the brown and green modules are the most closely interacting.



# Step 9: Assessment of Network Edges


A table was displayed below the igraph output and it documents the nodes of the highest degree of betweeness and the module that the concerned node resides in.


To give the user an illustration of the distribution of the strength of the connections between nodes and allows them to selectively eliminate nodes that fall below a particular connectedness and observe the changes in the graph.




# Step 10: Incorporate external information
#sample data supplied in 'info_final.txt'
To upload an external-information file the style of row names must be common to that of the initially puloaded file. The rows can have been subjected to rearrangement between each file. 

After a network has been inferred, it is possible to then upload another file that characterizes the row-samples from the initially uploaded file. Such characterizations are in the form of numeric metrics on phenotypes that were experimentally measured in each of the taken samples when the data were being generated.


Initially a single traits correlation can be measured with respect to each module-Eigengene to find significant relationships on an individual level.

![](https://github.com/chrisclarkson/pics/blob/master/relating_phenotype1.png)

Subsequently if the uploaded sample-characterizing has more than one column i.e. a measure for more than one phenotype on each of the sample rows, then multiple columns can be selected and a Heatmap will be rendered documenting the correlation and accompanying P-value of each trait and each module Eigengene.

![](https://github.com/chrisclarkson/pics/blob/master/relating_phenotype2.png)


# Running-time
The time taken for this pipeline to output results is stratified according to the following file-sizes: <1 MB- 2minutes, <4 MB- 3minutes, <8 MB- 5minutes and <12 MB- 8 minutes.

