import cPickle, gzip, glob, os, scipy
import numpy as np
def pickle(fname, obj):
    cPickle.dump(obj=obj, file=gzip.open(fname, "wb", compresslevel=3), protocol=2)
def unpickle(fname):
    return cPickle.load(gzip.open(fname, "rb"))
    
import networkx as nx

from networkx.algorithms.traversal.depth_first_search import dfs_edges
from networkx.algorithms.isolate import isolates


class db_network(object):
    def __init__(self, graph):
        self.nodes = graph.nodes().keys()
        self.graph = graph

    def top_n_degree(self,n):
        degrees = [(node,val) for (node, val) in self.graph.degree()]
        degrees = sorted(degrees, key=lambda x: x[1], reverse=True)
        return degrees[0:n]
    
    def get_neighbors(self,node):
        neighbor_dict={}
        for ng in self.graph.neighbors(node):
            neighbor_dict[ng]=self.graph.get_edge_data(node,ng)
        return neighbor_dict
    
    def get_subsample(self,node_list,from_node_list_only): ## from_node_list_only is a option to prevent nodes not in the nodelist from being included
        sub_network = nx.MultiGraph()
        node_set = set(node_list)
        if not from_node_list_only:
            for node in node_list:
                sub_network.add_node(node)
                neighbor_dict=self.get_neighbors(node)
                for ng in neighbor_dict.keys():
                    if not sub_network.has_node(ng):
                        sub_network.add_node(ng)          
                    for ngg in neighbor_dict[ng].items():
                        sub_network.add_edge(node,ng,**ngg[1])
            return db_network(sub_network)
        else:
            for node in node_list:
                sub_network.add_node(node)
                neighbor_dict=self.get_neighbors(node)
                for ng in neighbor_dict.keys():
                    if ng in node_set:
                        sub_network.add_node(ng)          
                        for ngg in neighbor_dict[ng].items():
                            sub_network.add_edge(node,ng,**ngg[1])
            return db_network(sub_network)
                

    def sel_by_edge_criteria(self, edge_filter_dict):  ### returns a db_network object with the filtered properties 
        sub_network = nx.MultiGraph()
        filter_list = edge_filter_dict.keys()
        for node in self.nodes:
            sub_network.add_node(node)
        for edge, feat in self.graph.edges.items():
            bool_arr = []
            for fil in filter_list:
                bool_arr.append(edge_filter_dict[fil](feat[fil]))
            if np.all(bool_arr):
                sub_network.add_edge(edge[0],edge[1],**feat)
        sub_network.remove_nodes_from(list(isolates(sub_network)))
        return db_network(sub_network)
        

    def get_n_neighbors(self,node,depth):
        node_set = set()
        for edge in dfs_edges(self.graph,source=node, depth_limit=depth):
            node_set.add(edge[0])
            node_set.add(edge[1])
        return list(node_set)   

class db_node_dict(object):

    def __init__(self, prop_dict):
        self.prop_dict = prop_dict
        self.nodes = prop_dict.keys()
        if self.nodes:
            self.properties = prop_dict[self.nodes[0]].keys()
        else:
            self.properties = None
        
    def add_new_data(self, new_dict):
        for node in new_dict:
            if node in self.nodes:
                for prop in new_dict[node]:
                    self.prop_dict[node][prop]=new_dict[node][prop]
            else:
                self.nodes.append(node)
                self.prop_dict[node]={}
                for prop in new_dict[node]:
                    self.prop_dict[node][prop]=new_dict[node][prop]         
    
    def filter_by_prop(self, prop_filter_dict):
        select_nodes={}
        for node in self.nodes:
            bool_arr = []
            for prop in prop_filter_dict.keys():
                bool_arr.append(prop_filter_dict[prop](self.prop_dict[node][prop]))
            if np.all(bool_arr):
                select_nodes[node]=self.prop_dict[node]
        return db_node_dict(select_nodes)
    
    def update_info(self):
        if self.prop_dict:
            self.nodes = self.prop_dict.keys()
            self.properties = self.prop_dict[self.nodes[0]].keys()

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem

class rd_mol(object):

    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        self.mol_H = None
        self.prop_2d_collector = dict()
        self.prop_3d_collector = dict()
        self.fingerprint_collector = dict()
        self.rdkit_base_info = {'2D_Props': [t[0] for t in Descriptors.descList],
                               '3D_Props': ['Asphericity','Eccentricity','InertialShapeFactor',
                                            'NPR1','NPR2','PMI1','PMI2','PMI3',
                                            'RadiusOfGyration','SpherocityIndex']}
        
    def get_2d_properties(self,prop_list):
        for prop in prop_list:
            self.prop_2d_collector[prop]=eval('Descriptors.'+prop+'(self.mol)')
            
    def get_fingerprint(self) : #### Current supported fingerprint is MACCS
        self.fingerprint_collector['MACCS']=DataStructs.BitVectToText(MACCSkeys.GenMACCSKeys(self.mol))
    
    def gen_3d_conf(self,num_confs):
        self.mol_H = Chem.AddHs(self.mol)
        AllChem.EmbedMultipleConfs(self.mol_H, numConfs=num_confs)
        
    def get_3d_properties(self,prop_list):
        if self.mol_H:    
            for prop in prop_list:
                self.prop_3d_collector[prop]=eval('Descriptors3D.'+prop+'(self.mol_H)')
        else:
            self.gen_3d_conf(10)
            for prop in prop_list:
                self.prop_3d_collector[prop]=eval('Descriptors3D.'+prop+'(self.mol_H)')
     
from collections import defaultdict

def get_class_set(db_network,filter_criteria,sel_nodes_list,null_nodes_list,is_multilabel):
    node_set = defaultdict(set)
    filter_list = filter_criteria.keys()
    for node in sel_nodes_list:
        temp_nodes = db_network.get_neighbors(node)
        for ddd in temp_nodes:
            bool_arr=[]
            for dddd in temp_nodes[ddd]:
                for fil in filter_list:
                    bool_arr.append(filter_criteria[fil](temp_nodes[ddd][dddd][fil]))
            if np.any(bool_arr):
                node_set[node].add(ddd)
    for node in null_nodes_list:    
        temp_nodes = db_network.get_neighbors(node)
        for ddd in temp_nodes:
            node_set['null'].add(ddd)
    sel_drug_list = set()
    for node in sel_nodes_list:
        sel_drug_list.update(node_set[node])
    null_drug_list = list(node_set['null'] - set(sel_drug_list).intersection(node_set['null']))
    sel_drug_list = list(sel_drug_list)
    invert_set = defaultdict(list)
    for node in sel_nodes_list:
        for mol in node_set[node]:
            invert_set[mol].append(node)
    for mol in null_drug_list:
        invert_set[mol].append('null')        
    if is_multilabel:   ## if not multibabel, then multiclass
        X_mol,y_labels = [],[]
        for x,y in invert_set.items():
            X_mol.append(x)
            y_labels.append(y)
    else:
        X_mol,y_labels = [],[]
        for x,y in invert_set.items():
            if len(y)==1:
                X_mol.append(x)
                y_labels.append(y)
    return X_mol,y_labels
        
from sklearn.pipeline import Pipeline
from sklearn import base


class BitVectorizer(base.BaseEstimator, base.TransformerMixin):
    
    def __init__(self,rd_mol_dict):
        self.fp_type = 'MACCS'
        self.finger_dict = rd_mol_dict
#        return None
  # We will need these in transform()
    
    def fit(self, X, y=None):
        # This transformer doesn't need to learn anything about the data,
        # so it can just return self without any further processing
        return self
    
    def transform(self, X):
        X_col = []
        for xi in X:
            X_col.append(map(int,self.finger_dict[xi][self.fp_type]))
        return X_col
        # Return an array with the same number of rows as X and one
        # column for each in self.col_names
        
        
class FeatureCollector(base.BaseEstimator, base.TransformerMixin):
    
    def __init__(self,feature_names,feature_dicts):
        self.feature_names = feature_names
        self.feature_dicts = feature_dicts
#        return None
  # We will need these in transform()
    
    def fit(self, X, y=None):
        # This transformer doesn't need to learn anything about the data,
        # so it can just return self without any further processing
        return self
    
    def transform(self, X):
        X_col = []
        for xi in X:
            feat_list = []
            for feat in self.feature_names:
                feat_list.append(self.feature_dicts[xi][feat])
            X_col.append(feat_list)
        return X_col
        
        
if __name__=="__main__":
    print True

 