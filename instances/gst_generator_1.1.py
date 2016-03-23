# -*- coding: utf-8 -*-
"""
--------------------------------------
Group Steiner Tree Problem generator
-------------------------------------
DESCRIPTION:


"""

# using networkX interface for graph
import networkx as nx
import random
import sys,getopt
from collections import defaultdict
from itertools import *

class GST_Generator:
    """GST instance generator"""
          
        
    def __init__(self,**kwargs):
        """generate GST problem instance with V vertices, E edges, and k number of Steiner groups"""
        args = set([str(x) for x in kwargs.keys()])
        self.G = nx.Graph()
        self._groups_of_node = defaultdict(set)
        self._nodes_of_group = defaultdict(set)
                            
        self._cmd_call_help = '\rpython gst_generator_1.01.py --file=<stp_file>\n \rpython gst_generator_1.01.py --nodes=<V> --edges==<E> --groups=<k> --bnd=<gamma> --graph=<connected|small-world>'
        if 'stp_file' in args:
              self.stp2graph(kwargs['stp_file'])
              self.fig_name = kwargs['stp_file'][-3]
              self.draw()
              return
        
        _needed_args = set(['nodes','edges','groups', 'bound', 'graph'])
        
       
        if _needed_args.issubset(args):
            self.V          = kwargs['nodes']
            self.E          = kwargs['edges']
            self.k          = kwargs['groups']
            self.bnd        = kwargs['bound']
            self.graph_type = kwargs['graph']
        else:
            print('Some parameters are missing')
            print(self._cmd_call_help)
            raise ValueError
            sys.exit(1)
            
      
    
           
        # check graph type
        if self.graph_type=='connected':
            # generate random connected graph
            self.generate_rand_graph()
            
        elif self.graph_type=='small-world': 
            print("Generate small-world graph with %d vertices." % self.V)
            # watts strogatz small world graph with V vertices, E edges, and avg. degree 3
            self.G = nx.watts_strogatz_graph(self.V,3,0.2)
            self.E = self.G.number_of_edges()
		elif self.graph_type=="duin2"
			self.generate_rand_graph()
                        
            
        self.fig_name = "GST_%s_V=%d_E=%d_k=%d_bnd=%d.png" % (self.graph_type,self.V,self.E,self.k, self.bnd)
        self.file_name = "GST_%s_V=%d_E=%d_k=%d_bnd=%d.stp" % (self.graph_type,self.V,self.E,self.k, self.bnd)
        self._info_groups = ""
    
        if self.graph_type=="duin2"
			self.assign_edge_weights(1,100)
			freeGroupSize=random.uniform(0,0.5)
			
			
			
		elif
			self.assign_edge_weights(1,50)
			self.assign_node_groups()
        
        self.save2stp()
        
                
    def __repr__(self):
        """string representation of GSP instance"""
        return (nx.info(self.G)+self._info_groups)
       
            
    def generate_rand_graph(self):
        """method for generating weighted undirect graph with V vertices and E edges."""
        # self.G = nx.dense_gnm_random_graph(self.V,self.E)
        
        self.G.add_edge(0,1)
        for v in range(2,self.V):
            # pick random vertex
            u = random.choice(self.G.nodes())
            self.G.add_edge(u,v)
        
        print "cvorovi: ", self.G.nodes()
        
        if self.E-self.G.number_of_edges()>0:
            cG = nx.complement(self.G)        
            self.G.add_edges_from(random.sample(cG.edges(), self.E-self.G.number_of_edges()))
            
        
        if nx.is_connected(self.G) : print("Generated graph is connected.")
       
        
        
    
    def assign_edge_weights(self,a,b):
        """assign weights to edges from set {a,a+1,...,b}, a<b."""
        weights = {e:random.randint(a,b) for e in self.G.edges()}
        nx.set_edge_attributes(self.G,"weight",weights)
        
        
    
    def assign_node_groups(self):
        """      
            sample vertices in k groups
            
        """
        self.nodes = self.G.nodes()
        self.groups = range(0,self.k)
        print self.nodes
        self._nodes_of_group = {g : set() for g in self.groups}
        self._groups_of_node = {v : set() for v in self.nodes}
        
        # for every group place atmost gamma vertices in it        
        for g in self.groups:
            nodes_sample = random.sample(self.nodes,random.randint(1,self.bnd))
    
            # update groups of v
            for v in nodes_sample: self._groups_of_node[v].add(g)
            self._nodes_of_group[g] = nodes_sample
                  
            
        self.bnd = max(map(len,self._nodes_of_group.values()))
        print self.bnd
        print self._nodes_of_group
            
                   
            
        
        
    def get_group_by_node(self,v):
        """identify group of node"""
        return self._nodes_of_group[v]

    def get_nodes_by_group(self,g):
        """"return list of nodes in some group"""
        return self._group_of_node[g]
        
        
    def draw(self):
        """Draw graph representation."""
        
             
        self._info_groups = (
                        "\n -------------------------\n" + 
                         repr(self._nodes_of_group) + "\n" +
                        "Bound on the size of groups: " + str(self.bnd)
                        )   
                        
        print(nx.info(self.G)+self._info_groups)
      
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8,8))
        pos=nx.random_layout(self.G) # positions for all nodes
        
        # nodes within same group are of the same color
        for g in self._nodes_of_group.keys():
            
            nx.draw_networkx_nodes(self.G,pos,
                       nodelist = self._nodes_of_group[g],
                       node_color = [random.random()]*len(self._nodes_of_group[g]),
                       vmin = 0,
                       vmax = 1,
                       node_size=700,
                   alpha=0.5)
                   
        # edges
        nx.draw_networkx_edges(self.G,pos,width=1.0,alpha=0.5)
        
           # labels
        group_labels = {v:str(list(self._groups_of_node[v])).strip('[]') for v in self._groups_of_node.keys()}
        
        loc_group_labels = zip(pos.values(),group_labels.values())
            
      
        for _pos,_l in loc_group_labels:
            
            plt.text(_pos[0],_pos[1]+0.05,s=_l)        
        
        nx.draw_networkx_labels(self.G,pos,font_size=14,font_family='sans-serif')
        edge_labels = {(e[0],e[1]):self.G[e[0]][e[1]]['weight'] for e in self.G.edges()}
        nx.draw_networkx_edge_labels(self.G,pos,edge_labels)
        nx.draw
        
        # save graph image        
        plt.axis('off')
        plt.savefig(self.fig_name) # save as png
        plt.show() # display       
          
       
        
        
    def stp2graph(self,stp_file):
        """ read data from stp file    """
        
        if not stp_file.endswith('.stp'): 
            print('File is missing.')
            sys.exit(1)
            
        self.G = nx.Graph()
          
        
        stp_handler = open(stp_file,'r')
        g = 0
        
        for line in stp_handler:
                            
            if line.startswith('Nodes'):
                self.V = line.partition(' ')[2]
            if line.startswith('Edges'):
                self.E = line.split()[1]
            
            
            if line.startswith('E '):
                _args = line.split()
                u,v,wt = int(_args[1]),int(_args[2]), float(_args[3]) 
                self.G.add_edge(u,v,weight=wt)
                
            if line.startswith('SECTION Terminals'):
                self.nodes = self.G.nodes()
                              
                
            if line.startswith('T '):
                _args = line.split()
              
                for x in _args[1:]:
                    self._nodes_of_group[g].add(int(x))              
                                            
                 # assign group to node
                for v in self._nodes_of_group[g]:
                    self._groups_of_node[v].add(g)
                 
                g += 1
        # close file handler
        stp_handler.close() 
        
        self.bnd = max(len(_nodes) for _nodes in self._nodes_of_group.values() )
        
        print "podaci:"
        print self.bnd
        
        print self._nodes_of_group
        
                            
                
            
            
    
    
    def save2stp(self):
        """
            output GST to .stp to file format            
        """
        stp_handler = open(self.file_name,'w')
        
        stp_handler.write("   0 string 33D32945  STP Steiner Tree Problem File\n")
        stp_handler.write("SECTION Comment\n")
        stp_handler.write("Name 'GSP instanca'\n")
        stp_handler.write("Remark 'First example for GST problem.'\n")
        stp_handler.write("END\n\n")
        
        # problem description
        stp_handler.write("SECTION Graph\n")
        stp_handler.write("Nodes %d\n" % self.G.order())
        stp_handler.write("Edges %d\n" % self.G.number_of_edges())
    
        # output edges
        for e in self.G.edges():
            u,v = e[0],e[1]
            w = self.G[u][v]['weight']            
            stp_handler.write("E %d %d %f\n" % (u,v,w))
        stp_handler.write("END\n\n")
        
        stp_handler.write("SECTION Terminals\n")
        stp_handler.write("Terminals %d\n" % self.k)
        
        for g in self._nodes_of_group.values():
            stp_handler.write("T ")
            for v in g : stp_handler.write("%d " % v)
            stp_handler.write("\n")
        stp_handler.write("END\n\n")
        stp_handler.write("EOF")
        
        # end of file
        stp_handler.close()
        
        
            
           

    
    
def main(argv):
    """ parsing command-line arguments """
    try:
        opts,args = getopt.getopt(argv,"ht:V:E:g:b:f:",
                                  ["help","type=","vertices=","edges=","groups=","bound=","file="])
    except getopt.GetoptError:
        print "python gst_generator.py -V <int> -E <int> -k <int> -b <None|int>"
        sys.exit(1)
        
    file_stp = None    
    for opt,arg in opts:
        if opt == '-h':
            print "python gst_generator.py -t <connected>   -V <int> -E <int> -k <int> -b <None|int> -f <stp_file>"
            print "python gst_generator.py -t <small-world> -V <int>          -k <int> -b <None|int> -f <stp_file>"
            sys.exit(1)
        elif opt in ("-T","--type"):
            graph_type = str(arg)
            print graph_type
        elif opt in ("-V","--vertices"):
            V = int(arg)
        elif opt in ("-E","--edges"):
            E = int(arg)
        elif opt in ("-g","--groups"):
            k = int(arg)
        elif opt in ("-b","--bound"):
            gamma = int(arg)
        elif opt in ("-f","--file"):
            file_stp = str(arg)
            

    # define GST instance
    if file_stp:
        GSP_instance = GST_Generator(stp_file=file_stp)
        
    else:
        print "\t"+graph_type
        GSP_instance = GST_Generator(nodes=V,edges=E,groups=k,bound=gamma,graph=graph_type)    
        GSP_instance.draw()
        
    
if __name__ == '__main__':
    main(sys.argv[1:])
    
            
    
    
    
    