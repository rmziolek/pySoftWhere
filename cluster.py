import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances

class largest_cluster:
    """Returns the resids of molecules in the largest cluster in the system (starting from 0!). THIS ACTUALLY RETURNS ANYTHING BIGGER THAN 1... FIX  """ 

    
    def __init__(self, u, cluster_atoms, cluster_cutoff,frame):

        self.u = u
        self.cluster_atoms = cluster_atoms
        self.cluster_cutoff = cluster_cutoff
        
        self.frame = frame
       
        self._cluster_gt()

        
    def _cluster_gt(self):
        
        
        """
        do the clustering
        """
#        for ts in self.u.trajectory[self.frame:self.frame+1:]: 
        ts=self.u.trajectory[self.frame]
        cluster_sizes=[]

        box_ = self.u.dimensions

        dist_array=mda.analysis.distances.contact_matrix(self.cluster_atoms.atoms.positions,cutoff=self.cluster_cutoff,box=box_)
        pairs=[]
        source_id=0
        for source in dist_array:

            for target in range(len(dist_array)):
                if source[target]==True and source_id!=target:
                    pairs.append([source_id,target])

            source_id+=1

        G=nx.Graph()
        for pair in pairs:
            G.add_edge(int(pair[0]),int(pair[1]))


        clusters=[h.nodes for h in nx.connected_component_subgraphs(G)]

        cluster_sizes=[len(h) for h in nx.connected_component_subgraphs(G)]

        no_ones=len(self.cluster_atoms)-sum(cluster_sizes)

        [cluster_sizes.append(1) for i in range(no_ones)]
        
        self.clusters=clusters
        
    def __call__(self):
        
        """cluster.
        """
        return self.clusters


###################


class cluster_distribution:
    """Returns the number of molecules in each cluster""" 

    
    def __init__(self, u, cluster_atoms, cluster_cutoff,frame):

        self.u = u
        self.cluster_atoms = cluster_atoms
        self.cluster_cutoff = cluster_cutoff
        
        self.frame = frame
       
        self._cluster_gt()

        
    def _cluster_gt(self):
        
        
        """
        do the clustering
        """
#        for ts in self.u.trajectory[self.frame:self.frame+1:]: 
        ts=self.u.trajectory[self.frame]
        cluster_sizes=[]

        box_ = self.u.dimensions

        dist_array=mda.analysis.distances.contact_matrix(self.cluster_atoms.atoms.positions,cutoff=self.cluster_cutoff,box=box_)
        pairs=[]
        source_id=0
        for source in dist_array:

            for target in range(len(dist_array)):
                if source[target]==True and source_id!=target:
                    pairs.append([source_id,target])

            source_id+=1

        G=nx.Graph()
        for pair in pairs:
            G.add_edge(int(pair[0]),int(pair[1]))


        clusters=[h.nodes for h in nx.connected_component_subgraphs(G)]

        cluster_sizes=[len(h) for h in nx.connected_component_subgraphs(G)]

        no_ones=len(self.cluster_atoms)-sum(cluster_sizes)

        [cluster_sizes.append(1) for i in range(no_ones)]
        
        self.clusters=clusters
        
        self.cluster_sizes = cluster_sizes
         
    def __call__(self):
        
        """cluster.
        """
        return self.cluster_sizes





