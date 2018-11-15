# Graph Theory and Networks for Drug Discovery 

![iris](https://raw.githubusercontent.com/Vanabins28/db_network/images/P3.png)


Why Graphs?
----------
Graph data structures provide for a way to store and represent complex, heterogenous data sources. Unlike other data structures, graphs focuses on relationships between objects. This is different from relational data models (SQL) where relationships between two different objects have to be inferred by the use of foreign keys to perform complex joins. For applications that require quick lookup of long-distance interactions, graph data structures provide an improvement over relational data models. 


Graph data structures are naturally well-suited to be applied to represent and store the complex relationships between drugs/compounds their protein targets and their associated diseases. 

![iris](https://raw.githubusercontent.com/Vanabins28/db_network/images/Drug_example.png)


Applications in Drug Discovery
----------
Aside from a convenient way to represent and store complex data, multiple algorithms have been developed to study graphs. These methods include graph traversal methods (BFS, DFS, Dijkstra's algorithm etc...) for searching specific relationships between objects in the graphs, community detection (Girvan-Newman Algorithm), clustering, label-propagation, and anomaly detection. 
These methods are already being used in drug discovery. For example anomaly detection of drug-target networks can be used to identify drug molecules that may display polypharmacological activity. Community detection and clustering can be used to identify groups of protein targets that are strongly associated with a disease. Label propagation can be used as a semi-supervised learning approach to predict drug and side-effect interactions.
Since graph data structures already store relationships between drugs/compounds, protein targets, and diseases, graphs can be used as a graph database that can be used to quickly retrieve select protein targets and the compounds that bind to them. Graph-databases like Neo4j have already been applied to store data related to drug discovery.
Finally graphs can be used to easily and clerly visualize the complex relationships of the drug/target network.

![iris](https://raw.githubusercontent.com/Vanabins28/db_network/images/P2.png)


This repo describes my work in applying graph and network based methods to drug discovery. Using data from public sources like the ChemBL, MedDRA, UniProt, and OMIM, we can store the complex relationships between these different data sources as a graph database and apply the methods mentioned above to help drug discovery.


Requires
--------
* `python 2.7`
* `networkx`
* `numpy`
* `rdkit`
* `scipy`
* `sklearn`

Usage
---------
jupyter notebooks are provided that show the basic usage of the db_network model. These can be used to generate datasets for training, for visualization, and for anomaly detection
