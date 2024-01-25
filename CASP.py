import pandas as pd
import dowhy
from dowhy import CausalModel

import numpy as np
import pandas as pd
import graphviz
import networkx as nx

from cdt.causality.graph import LiNGAM

import networkx as nx

def make_graph(adjacency_matrix, labels=None):
    idx = np.abs(adjacency_matrix) > 0.01
    dirs = np.where(idx)
    d = graphviz.Digraph(engine='dot')
    names = labels if labels else [f'x{i}' for i in range(len(adjacency_matrix))]
    for name in names:
        d.node(name)
    for to, from_, coef in zip(dirs[0], dirs[1], adjacency_matrix[idx]):
        d.edge(names[from_], names[to], label=str(coef))
    return d

def str_to_dot(string):
    '''
    Converts input string from graphviz library to valid DOT graph format.
    '''
    graph = string.strip().replace('\n', ';').replace('\t','')
    graph = graph[:9] + graph[10:-2] + graph[-1] # Removing unnecessary characters from string
    return graph


class CASP:
    def __init__(self, adata):
        self.adata = adata 
        self.gene_dict = {item.split(':')[0].split(' - ')[0]: item for item in self.adata.var_names}
        self.idx_dict = {item.split(':')[0].split(' - ')[0]: i for i, item in enumerate(self.adata.var_names)}
        self.df = pd.DataFrame(adata.X, columns=[item.split(':')[0].split(' - ')[0] for item in self.adata.var_names])
        self.graph = None
        self.casul_model = None
        self.identified_estimand = None
        
    def causal_discovery(self, show=True):
        graphs = {}
        labels = [f'{col}' for i, col in enumerate(self.df.columns)]
        functions = {
            'LiNGAM' : LiNGAM
        }

        for method, lib in functions.items():
            obj = lib()
            output = obj.predict(self.df)
            adj_matrix = nx.to_numpy_array(output)
            adj_matrix = np.asarray(adj_matrix)
            
            graph_dot = make_graph(adj_matrix, labels)
            graphs[method] = graph_dot

        # Visualize graphs
        for method, graph in graphs.items():
            print("Method : %s"%(method))
            if show:
                display(graph)
        self.graph = output
        self.graph_io = graph
        
    def get_ancestors_and_weights(self, node):
        ancestors = nx.ancestors(self.graph, node)
        ancestor_weights = {}
        for ancestor in ancestors:
            if self.graph.has_edge(ancestor, node):
                weight = self.graph[ancestor][node]['weight']
                ancestor_weights[ancestor] = weight
        return ancestor_weights
    
    def get_all_neighbors(self, node):
        predecessors = list(self.graph.predecessors(node))
        successors = list(self.graph.successors(node))
        all_neighbors = predecessors + successors
        return all_neighbors

    def get_leaf_nodes(self):
        leaf_nodes = [node for node in self.graph.nodes() if self.graph.out_degree(node) == 0]
        leaf_index = [self.idx_dict[node] for node in leaf_nodes]
        return leaf_nodes, leaf_index


    def refute_relation(self, treatment, outcome, method='random'):
        
        causal_model = CausalModel(self.df, treatment=treatment, outcome=outcome, graph=str_to_dot(self.graph_io.source))

        self.causal_model = causal_model
        
        identified_estimand  = causal_model.identify_effect()
        
        self.identified_estimand = identified_estimand
        
        estimate = causal_model.estimate_effect(identified_estimand,
            method_name="backdoor.linear_regression",
            test_significance=True,
            #confidence_intervals=True
            target_units="ate"
        )

        print(estimate)
        
        if method == 'random':     
            refutation = causal_model.refute_estimate(identified_estimand, estimate, method_name="random_common_cause",random_seed=42)

        if method == 'subset':
            refutation = causal_model.refute_estimate(identified_estimand, estimate,
                method_name="data_subset_refuter", subset_fraction=0.9, random_seed=42)

        if method == 'placebo':
            refutation = causal_model.refute_estimate(identified_estimand, estimate, method_name="placebo_treatment_refuter",
                            placebo_type="permute", num_simulations=20, random_seed=42)
        
        print(refutation)
        
        return refutation
    
    def get_sub_graph(self, node):
        neighbors = self.get_all_neighbors(self.graph, node)
        print(f"All neighbors of node {node}: {neighbors}")

        neighbors.append(node)

        subgraph = output.subgraph(neighbors)

        new_graph = nx.DiGraph(subgraph)

        '''name_dict = {'CD44':'CD44', 'p53':'TP53', 'T-bet':'TBX21', 
                'CD4':'CD4', 'CD21':' CR2', 
                'CD5':'CD5', 'IDO-1':'IDO-1', 
                'CD11b':'ITGAM', 'CD56':'NCAM1', 
                'PD-1':'PDCD1', 'HOECHST1': 'HOECHST1',
                'GFAP': 'GFAP', 'MMP12':'MMP12'}'''
                
        '''name_dict = {'p53':'TP53', 'Vimentin':'VIM', 'Beta catenin':'CTNNB1', 
                'HLA-DR':'HLA-DRA', 'CD45':'PTPRC', 
                'H3K9ac':'H3K9ac', 'Pan-Keratin':'Pan-Keratin', 
                'H3K27me3':'H3K27me3', 'HLA_Class_1':'HLA_Class_1', 
                'dsDNA':'dsDNA', 'SMA': 'SMN1',
                'Ki67': 'MKI67', 'CD138':'SDC1', 'CD4':'CD4', 'CD11b': 'ITGAM'}'''

        labels = [self.name_dict[term] for term in list(new_graph.nodes())]

        adj_matrix = nx.to_numpy_array(new_graph)
        adj_matrix = np.asarray(adj_matrix)

        graph_dot = make_graph(adj_matrix, labels)

        display(graph_dot)
        
    def linear_refute(self):
        res_unobserved_range = self.causal_model.refute_estimate(self.identified_estimand, self.estimate, method_name="add_unobserved_common_cause",
                                           confounders_effect_on_treatment="binary_flip", confounders_effect_on_outcome="linear",
                                           effect_strength_on_treatment=[0.001, 0.005, 0.01, 0.02],
                                           effect_strength_on_outcome=[0.001, 0.005, 0.01,0.02])
        print(res_unobserved_range)
    