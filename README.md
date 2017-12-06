# WGCNA_shiny



<p>
    <img src="pics/Poster.png" />
</p>

Turning the main features of WGCNA into an interactive shiny page allowing for co-expression of cluster identification, network visualisation and correlation of eigengeme values to phenotypes.


You can visualise the network on the shiny page or have a more dynamic network visualisation output can be saved.

# Methodology


Necessary file format
The initially uploaded file (containing gene-expression data or compositional microbiotic co-abundance data) must contain the genes/genera that are to be correlated and clustered in the columns and the samples must be specified in the rows. For labeling the nodes and to maximize the effectiveness of the displayed network the file must contain a header that specifies the column names, otherwise the nodes will be labeled independently.
![](pics/Poster.png)

To upload an external-information file the style of row names must be common between each file. The rows can have been subjected to rearrangement between each file. 


In the website sidebar, there is a navigation-panel that allows the user to specify if the files to be uploaded are text files or CSV files. Users can also specify how the files are delimited.


# Step 1: Filtering


The first step in the pipeline is to filter out samples and nodes that are not informative. The first round of HC is done to remove samples too far displaced from all others and would perturb the inference of clusters if included. The user will choose the threshold via an interactive slider.


# Step 2: Soft threshold selection


The command ‘pickSoftThreshold’ helps the user to choose a soft-threshold for the network to reach a scale-free topology. In accordance with the equation: 


sij=cor|(xi,xj)|

weights the adjacency-matrix and hence penalizes weak correlations. The value of beta is essential for the network to reach a scale-free topology.
The recommended threshold is output to the page. Furthermore, in the event that the user has an intuition that beta value should be different than the recommended power- the R2 fit-value to scale-free topology is plotted for each power as shown in Figure 5 below in “Results- Scale-Free topology dependent on Filtration“.

# Step 3: An adjacency-matrix is made via Pearson correlation.


After correlating the expression-values via Pearson-correlation, to assign weights to the interactions. This is implemented via the command: ‘adjacency’.


# Step 4: Cluster inference


The second round of HC allows for the detection of modules within the dataset. This is done with the command ‘TOMsimilarity’, which documents the number of shared neighbors among multiple nodes. Hence, instead of directly utilizing the pair-wise adjacencies from Pearson correlation, multi-node comparisons are used and stronger connections within the pair-wise correlation values are identified. Subsequently HC is performed to infer Topolgical Overlap and hence modules. The calculated similarity values are subtracted from 1 to convert them to dissimilarity measures. This is done to eliminate noise caused by spurious or erroneous adjacencies (Ai Li, S Horvath).
 
# Figure 2: 

Part A- A filtration cutoff level is selected and marked via a red line. Any samples that fall above that are line are considered outliers and are omitted from the dataset.


Part B- A heatmap represents the adjacency matrix produced by HC of the filtered dataset where the columns are treated with Pearson correlation to produce an adjacency matrix.


Part C- A soft-threshold is selected via a plot that documents the R2 fit to scale-free topology. Above the soft threshold value for reaching scale-free topology can be selected in accordance with the above graph where the R2 fit-value is plotted against each individual threshold-value.


Part D- HC of the TOM dissimilarity values identifies modules and colours are assigned to each simulated module.


# Step 5: Network-display


Using the package ‘igraph’ the network is plotted and the cluster members are displayed. The module colors of WGCNA-assigned module clusters are displayed in the network.


# Node-size and Edge-thickness


These two attributes were governed both in igraph and networkD3 by ‘enrichment factors’ whereby a user can choose a value to multiply the logarithmic degree per node and edge-strength to enrich the node-size and edge-thickness consecutively.


# Layout Algorithms


Both in igraph and networkD3, the node-placement/ layout algorithms are ‘force-directed’ this is to say that there is an automatic repulsion between each node and yet the edges that exist between them cause an attractive force. The stronger the edges between a pair of nodes, the more closely that pair is placed together and this promotes the clustering of closely interacting nodes.
However the selection is defaulted to ‘layout.auto’- a command that selects the optimal of three algorithm choices based on the number of nodes.


# Step 7: Interactive graph email


Using the package ‘networkD3’ an interactive network is made from the igraph data and stored in a HTML file, which can then be emailed to the user. The node labels are on permanent display and the user can zoom in and investigate particular nodes. Again node-sizes are in accordance with the number of edges that a node has. The opacity of the nodes is set such that both the node labels and the nodes themselves are simultaneously visible (Yang and Liu, 2013).


# Step 8: Assign relatedness to modules


To relate the assigned modules to one another, a final round of clustering is applied to the Eigengenes that summarise the expression data of each cluster and this gives a measure of how related the modules are amongst each other. This can be compared to both the module-identifying dendogram and the network where module colours that are positioned closely will be shown to have close relation in this dendogram.


# Figure 3: Interactive networkD3 Package Output Display and Subsequent Relation of the Modules


Part A: displays a portion of the HTML output of the networkD3 package (as mentioned above in “Step 7”, a legend is necessary to translate the independent networD3 coloring-scheme to that of WGCNA). 



Part B: shows the hierarchical clustering of the corresponding Module-Eigengene values to show which are the most closely related. The red arrows illustrate how (as mentioned in “Step 8 “) the node-placement of networkD3 concurs with the clustering of the Eigengene-values in showing that the brown and green modules are the most closely interacting.



# Step 9: Assessment of Network Edges


A table was displayed below the igraph output and it documents the nodes of the highest degree of betweeness and the module that the concerned node resides in.


To give the user an illustration of the distribution of the strength of the connections between nodes and allows them to selectively eliminate nodes that fall below a particular connectedness and observe the changes in the graph.


# Figure 4: Network-edges assessment


# Step 10: Incorporate external information


After a network has been inferred, it is possible to then upload another file that characterizes the row-samples from the initially uploaded file. Such characterizations are in the form of numeric metrics on phenotypes that were experimentally measured in each of the taken samples when the data were being generated.


Initially a single traits correlation can be measured with respect to each module-Eigengene to find significant relationships on an individual level. These relationships are measure as below in Figures 16 and 17 in “Results”.


Subsequently if the uploaded sample-characterizing has more than one column i.e. a measure for more than one phenotype on each of the sample rows, then multiple columns can be selected and a Heatmap will be rendered documenting the correlation and accompanying P-value of each trait and each module Eigengene.


# Running-time


The time taken for this pipeline to output results is stratified according to the following file-sizes: <1 MB- 2minutes, <4 MB- 3minutes, <8 MB- 5minutes and <12 MB- 8 minutes.

# Figure 9: Network and subsequent Barplot of inferred correlation between Module-Eigengenes with a cutoff-threshold of 700 (the lowest possible threshold).  

Figure 9: A barchart indicates the significance of each cluster-module with regard to the given trait. Unlike in the case where a cutoff-threshold of 10,000 (which resulted in a much larger and more sparse dataset) was used, the colors the discovered modules show much higher gene-significance values with each of the given traits, in this case the time-elapsed difference in mgCO2 levels in the sampled populations, measurements taken after 7 and 14 days.
