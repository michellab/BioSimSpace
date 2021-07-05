# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP: Graph generation
=====

Alchemical free energy calculations hold increasing promise as an aid to drug 
discovery efforts. However, applications of these techniques in discovery 
projects have been relatively few, partly because of the difficulty of planning 
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an 
automated algorithm to plan efficient relative free energy calculations between 
potential ligands within a substantial of compounds.

"""

# *****************************************************************************
# Lomap2: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2016  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
#
# This part of the code has been originally made by Jonathan Redmann, 
# and Christopher Summa at Summa Lab, Dept. of Computer Science, 
# University of New Orleans and it has just been adapded to the new Lomap code
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see http://www.gnu.org/licenses/
# *****************************************************************************


# ****************
# MODULE IMPORTS
# ****************

import networkx as nx
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import copy
from operator import itemgetter
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import os.path
import logging
import tempfile
import shutil
import traceback

__all__ = ['GraphGen']


# *************************
# Graph Class
# *************************


class GraphGen(object):
    """
    This class is used to set and generate the graph used to plan
    binding free energy calculation
    """

    def __init__(self, dbase):

        """
        Inizialization function
    
        Parameters
        ----------

        dbase : dbase object
            the molecule container
       
        """

        self.dbase = dbase

        self.maxPathLength = dbase.options.max

        self.maxDistFromActive = dbase.options.max_dist_from_actives

        self.similarityScoresLimit = dbase.options.cutoff

        self.requireCycleCovering = not dbase.options.allow_tree

        if dbase.options.radial:
            self.lead_index = self.pick_lead()
        else:
            self.lead_index = None

        # A set of nodes that will be used to save nodes that are not a cycle cover for a given subgraph
        self.nonCycleNodesSet = set()

        # A set of edges that will be used to save edges that are acyclic for given subgraph
        self.nonCycleEdgesSet = set()

        # A count of the number of nodes that are not within self.maxDistFromActive edges
        # of an active
        self.distanceToActiveFailures = 0

        # Draw Parameters

        # THIS PART MUST BE CHANGED

        # Max number of displayed chemical compound images as graph nodes
        self.max_images = 2000

        # Max number of displayed nodes in the graph
        self.max_nodes = 100

        # The maximum threshold distance in angstroms unit used to select if a molecule is depicted
        self.max_mol_size = 50.0

        self.edge_labels = True

        # The following Section has been strongly copied/adapted from the original implementation

        # Generate a list related to the disconnected graphs present in the initial graph 
        if dbase.options.fast and dbase.options.radial:
            # only enable the fast map option if use the radial option
            self.initialSubgraphList = self.generate_initial_subgraph_list(fast_map=True)
        else:
            self.initialSubgraphList = self.generate_initial_subgraph_list()

        # A list of elements made of [edge, weights] for each subgraph
        self.subgraphScoresLists = self.generate_subgraph_scores_lists(self.initialSubgraphList)

        # Eliminates from each subgraph those edges whose weights are less than the hard limit
        self.remove_edges_below_hard_limit()

        # Make a new master list of subgraphs now that there may be more disconnected components
        self.workingSubgraphsList = self.generate_working_subgraphs_list()

        # Make a new sorted list of [edge, weights] for each subgraph now that there may be new subgraphs
        self.workingSubgraphScoresLists = self.generate_subgraph_scores_lists(self.workingSubgraphsList)

        # Remove edges, whose removal does not violate constraints, from the subgraphs,
        # starting with lowest similarity score first

        if dbase.options.fast and dbase.options.radial:
            # if we use the fast and radial option, just need to add the surrounding edges from the initial graph
            self.resultGraph = self.add_surrounding_edges()
            # after adding the surround edges, some subgraphs may merge into a larger graph and so need to update the
            # current subgraphs
            # self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)
            # merge all Subgraphs together for layout
            # self.resultGraph = self.merge_all_subgraphs()
        else:
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>ISSUE ORDER PROBLEM<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            self.minimize_edges()
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>ISSUE ORDER PROBLEM<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            # Collect together disjoint subgraphs of like charge into subgraphs
            self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)

            # Combine separate subgraphs into a single resulting graph
            self.resultGraph = self.merge_all_subgraphs()

            # Make a copy of the resulting graph for later processing in connectResultingComponents()
            self.copyResultGraph = self.resultGraph.copy()

            # Holds list of edges that were added in the connect components phase
            self.edgesAddedInFirstTreePass = []

            # Add edges to the resultingGraph to connect its components
            self.connect_subgraphs()

        return

    def pick_lead(self):
        if (self.dbase.nums() * (self.dbase.nums() - 1) / 2) != self.dbase.strict_mtx.size:
            raise ValueError("There are errors in the similarity score matrices")
        if not self.dbase.options.hub == "None":
            # hub radial option. Use the provided reference compound as a hub
            hub_index = None
            for i in range(0, self.dbase.nums()):
                if os.path.basename(self.dbase[i].getName()) == self.dbase.options.hub:
                    hub_index = i
            if hub_index is None:
                logging.info(
                    "Warning: the specified center ligand %s is not in the ligand database, will not use the radial option." % self.dbase.options.hub)
            return hub_index
        else:
            # complete radial option. Pick the compound with the highest total similarity to all other compounds to use as a hub
            all_sum_i = []
            for i in range(0, self.dbase.nums()):
                sum_i = 0
                for j in range(0, self.dbase.nums()):
                    sum_i += self.dbase.strict_mtx[i, j]
                all_sum_i.append(sum_i)
            max_value = max(all_sum_i)
            max_index = [i for i, x in enumerate(all_sum_i) if x == max_value]
            max_index_final = max_index[0]
            return max_index_final

    def generate_initial_subgraph_list(self, fast_map=False):

        """
        This function generates a starting graph connecting with edges all the 
        compounds with a positive strict similarity score  
        
        Returns
        -------
        
        initialSubgraphList : list of NetworkX graph
            the list of connected component graphs     
        
        """
        compound_graph = nx.Graph()

        if (self.dbase.nums() * (self.dbase.nums() - 1) / 2) != self.dbase.strict_mtx.size:
            raise ValueError("There are errors in the similarity score matrices")

        if not fast_map:
            # if not fast map option, connect all possible nodes to generate the initial graph
            for i in range(0, self.dbase.nums()):
                if i == 0:
                    compound_graph.add_node(i, ID=self.dbase[i].getID(),
                            fname_comp=os.path.basename(self.dbase[i].getName()),
                            active=self.dbase[i].isActive())

                for j in range(i + 1, self.dbase.nums()):

                    if i == 0:
                        compound_graph.add_node(j, ID=self.dbase[j].getID(),
                            fname_comp=os.path.basename(self.dbase[j].getName()),
                            active=self.dbase[j].isActive())

                    wgt = self.dbase.strict_mtx[i, j]

                    if wgt > 0.0:
                        compound_graph.add_edge(i, j, similarity=wgt, strict_flag=True)
        else:
            # if fast map option, then add all possible radial edges as the initial graph
            for i in range(0, self.dbase.nums()):
                # add the node for i
                compound_graph.add_node(i, ID=self.dbase[i].getID(),
                                        fname_comp=os.path.basename(self.dbase[i].getName()))
                if i != self.lead_index:
                    wgt = self.dbase.strict_mtx[i, self.lead_index]
                    if wgt > 0:
                        compound_graph.add_edge(i, self.lead_index, similarity=wgt, strict_flag=True)

        initialSubgraphGen = [compound_graph.subgraph(c).copy() for c in nx.connected_components(compound_graph)]
        initialSubgraphList = [x for x in initialSubgraphGen]

        return initialSubgraphList

    def generate_subgraph_scores_lists(self, subgraphList):

        """
        This function generate a list of lists where each inner list is the 
        weights of each edge in a given subgraph in the subgraphList, 
        sorted from lowest to highest 

        
        Returns
        -------
        
        subgraphScoresLists : list of lists
            each list contains a tuple with the graph node indexes and their 
            similatiry as weigth
        
        """

        subgraphScoresLists = []

        for subgraph in subgraphList:
            weightsDictionary = nx.get_edge_attributes(subgraph, 'similarity')

            subgraphWeightsList = [(edge[0], edge[1], weightsDictionary[edge]) for edge in weightsDictionary.keys()]

            subgraphWeightsList.sort(key=lambda entry: entry[2])

            subgraphScoresLists.append(subgraphWeightsList)

        return subgraphScoresLists

    def remove_edges_below_hard_limit(self):
        """
        
        This function removes edges below the set hard limit from each subGraph 
        and from each weightsList
        
        """

        totalEdges = 0

        for subgraph in self.initialSubgraphList:

            weightsList = self.subgraphScoresLists[self.initialSubgraphList.index(subgraph)]

            index = 0

            for edge in weightsList:

                if edge[2] < self.similarityScoresLimit:
                    subgraph.remove_edge(edge[0], edge[1])

                    index = weightsList.index(edge)

            del weightsList[:index + 1]

            totalEdges = totalEdges + subgraph.number_of_edges()

    def generate_working_subgraphs_list(self):
        """
        After the deletition of the edges that have a weigth less than the 
        selected threshould the subgraph maybe disconnected and a new master 
        list of connected subgraphs is genereted
        
        Returns
        -------
        
        workingSubgraphsList : list of lists
            each list contains a tuple with the graph node indexes and their 
            similatiry as weigth

        """

        workingSubgraphsList = []

        for subgraph in self.initialSubgraphList:

            newSubgraphList = [subgraph.subgraph(c).copy() for c in nx.connected_components(subgraph)]

            for newSubgraph in newSubgraphList:
                workingSubgraphsList.append(newSubgraph)

        return workingSubgraphsList

    def minimize_edges(self):
        """
        Minimize edges in each subgraph while ensuring constraints are met
        """

        for subgraph in self.workingSubgraphsList:

            weightsList = self.workingSubgraphScoresLists[self.workingSubgraphsList.index(subgraph)]

            # ISSUE ORDER IS ORIGINATED HERE
            # weightsList = sorted(weightsList, key = itemgetter(1))

            # This part has been copied from the original code
            self.nonCycleNodesSet = self.find_non_cyclic_nodes(subgraph)
            self.nonCycleEdgesSet = self.find_non_cyclic_edges(subgraph)
            numberOfComponents = nx.number_connected_components(subgraph)
            self.distanceToActiveFailures = self.count_distance_to_active_failures(subgraph)

            if len(subgraph.edges()) > 2:  # Graphs must have at least 3 edges to be minimzed

                for edge in weightsList:
                    if self.lead_index is not None:
                        # Here the radial option is appplied, will check if the remove_edge is connect to
                        # the hub(lead) compound, if the edge is connected to the lead compound,
                        # then add it back into the graph.
                        if self.lead_index not in [edge[0], edge[1]]:
                            subgraph.remove_edge(edge[0], edge[1])
                            if self.check_constraints(subgraph, numberOfComponents) == False:
                                subgraph.add_edge(edge[0], edge[1], similarity=edge[2], strict_flag=True)
                    elif edge[2] < 1.0:  # Don't remove edges with similarity 1
                        logging.info("Trying to remove edge %d-%d with similarity %f" % (edge[0],edge[1],edge[2]))
                        subgraph.remove_edge(edge[0], edge[1])
                        if self.check_constraints(subgraph, numberOfComponents) == False:
                            subgraph.add_edge(edge[0], edge[1], similarity=edge[2], strict_flag=True)
                        else:
                            logging.info("Removed edge %d-%d" % (edge[0],edge[1]))
                    else:
                        logging.info("Skipping edge %d-%d as it has similarity 1" % (edge[0],edge[1]))

    def add_surrounding_edges(self):
        """
        Add surrounding edges in each subgraph to make sure all nodes are in cycle
        """
        for subgraph in self.workingSubgraphsList:
            subgraph_nodes = subgraph.nodes()
            if self.lead_index in subgraph_nodes:
                # here we only consider the subgraph with lead compound
                self.nonCycleNodesSet = self.find_non_cyclic_nodes(subgraph)
                self.nonCycleEdgesSet = self.find_non_cyclic_edges(subgraph)
                for node in self.nonCycleNodesSet:
                    # for each node in the noncyclenodeset, find the similarity compare to all other surrounding nodes and pick the one with the max score and connect them
                    node_score_list = []
                    for i in range(0, self.dbase.nums()):
                        if i != node and i != self.lead_index:
                            node_score_list.append(self.dbase.strict_mtx[node, i])
                        else:
                            node_score_list.append(0.0)
                    max_value = max(node_score_list)
                    if max_value > self.similarityScoresLimit:
                        max_index = [i for i, x in enumerate(node_score_list) if x == max_value]
                        max_index_final = max_index[0]
                        subgraph.add_edge(node, max_index_final,
                                          similarity=self.dbase.strict_mtx[node, max_index_final], strict_flag=True)
                return subgraph

    def find_non_cyclic_nodes(self, subgraph):
        """
        Generates a list of nodes of the subgraph that are not in a cycle
         
        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for not cycle nodes

        Returns
        -------
        missingNodesSet : set of graph nodes
            the set of graph nodes that are not in a cycle
        
        """

        missingNodesSet = set()

        cycleNodes = []

        cycleList = nx.cycle_basis(subgraph)

        cycleNodes = [node for cycle in cycleList for node in cycle]

        missingNodesSet = set([node for node in subgraph.nodes() if node not in cycleNodes])

        return missingNodesSet

    def find_non_cyclic_edges(self, subgraph):
        """
        Generates a set of edges of the subgraph that are not in a cycle (called
        "bridges" in networkX terminology).
         
        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for not cycle nodes

        Returns
        -------
        missingEdgesSet : set of graph edges
            the set of edges that are not in a cycle
        
        """

        missingEdgesSet = set(nx.bridges(subgraph))

        return missingEdgesSet

    def check_constraints(self, subgraph, numComp):
        """
        Determine if the given subgraph still meets the constraints
        

        Parameters
        ----------
        subgraph : NetworkX subgraph obj
             the subgraph to check for the constraints 
        
        numComp : int
            the number of connected componets

        Returns
        -------
        constraintsMet : bool
           True if all the constraints are met, False otherwise
        """

        constraintsMet = True

        if not self.remains_connected(subgraph, numComp):
            constraintsMet = False

        # The requirement to keep a cycle covering is now optional
        if constraintsMet and self.requireCycleCovering:
            if not self.check_cycle_covering(subgraph):
                constraintsMet = False

        if constraintsMet:
            if not self.check_max_distance(subgraph):
                constraintsMet = False

        if constraintsMet:
            if not self.check_distance_to_active(subgraph):
                constraintsMet = False

        return constraintsMet

    def remains_connected(self, subgraph, numComponents):
        """
        Determine if the subgraph remains connected after an edge has been 
        removed
        
        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for connection after the edge deletition
        
        numComp : int
            the number of connected componets
        
        Returns
        -------
        isConnected : bool
            True if the subgraph is connected, False otherwise
            :param numComponents:

        """

        isConnected = False

        if numComponents == nx.number_connected_components(subgraph): 
            isConnected = True
        else:
            logging.info("Rejecting edge deletion on graph connectivity")

        return isConnected

    def check_cycle_covering(self, subgraph):
        """
        Checks if the subgraph has a cycle covering. Note that this has been extended from
        the original algorithm: we not only care if the number of acyclic nodes has
        increased, but we also care if the number of acyclic edges (bridges) has increased.
        Note that if the number of acyclic edges hasn't increased, then the number of
        acyclic nodes hasn't either, so that test is included in the edges test.
        
        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for connection after the edge deletion

        Returns
        -------
        hasCovering : bool
            True if the subgraph has a cycle covering, False otherwise

        """

        hasCovering = True

        # Have we increased the number of non-cyclic edges?
        if self.find_non_cyclic_edges(subgraph).difference(self.nonCycleEdgesSet): 
            hasCovering = False
            logging.info("Rejecting edge deletion on cycle covering")

        return hasCovering

    def check_max_distance(self, subgraph):
        """
        Check to see if the graph has paths from all compounds to all other 
        compounds within the specified limit

        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for the max distance between nodes

        Returns
        -------
        withinMaxDistance : bool
            True if the subgraph has all the nodes within the specified 
            max distance
        """

        withinMaxDistance = True

        for node in subgraph:
            eccentricity = nx.eccentricity(subgraph, node)
            if eccentricity > self.maxPathLength: 
                withinMaxDistance = False
                logging.info("Rejecting edge deletion on graph diameter for node %d" % (node))

        return withinMaxDistance

    def count_distance_to_active_failures(self, subgraph):
        """
        Count the number of compounds that don't have a minimum-length path to an active
        within the specified limit

        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for the max distance between nodes

        Returns
        -------
        failures : int
            Number of nodes that are not within the max distance to any active node
        """

        failures = 0

        hasActives=False
        for node in subgraph.nodes():
            if (subgraph.nodes[node]["active"]):
                hasActives=True
        if (not hasActives):
            return 0     # No actives, so don't bother checking

        paths = nx.shortest_path(subgraph)
        for node in subgraph.nodes():
            if (not subgraph.nodes[node]["active"]):
                ok=False
                for node2 in subgraph.nodes():
                    if (subgraph.nodes[node2]["active"]):
                        pathlen = len(paths[node][node2]) - 1   # No. edges is 1 less than no. nodes
                        if (pathlen <= self.maxDistFromActive): ok=True
                if (not ok): 
                    failures = failures + 1

        return failures

    def check_distance_to_active(self, subgraph):
        """
        Check to see if we have increased the number of distance-to-active failures

        Parameters
        ---------
        subgraph : NetworkX subgraph obj
            the subgraph to check for the max distance between nodes

        Returns
        -------
        ok : bool
            True if we have not increased the number of failed nodes
        """

        count = self.count_distance_to_active_failures(subgraph)
        failed =  (count > self.distanceToActiveFailures)
        if (failed): logging.info("Rejecting edge deletion on distance-to-actives %d vs %d" % (count,self.distanceToActiveFailures))
        logging.info("Checking edge deletion on distance-to-actives %d vs %d" % (count,self.distanceToActiveFailures))
        return not failed


    def merge_all_subgraphs(self):
        """Generates a single networkx graph object from the subgraphs that have
        been processed

        Returns
        -------
        finalGraph : NetworkX graph obj
            the final graph produced merging all the subgraphs. The produced
            graph may have disconneted parts

        """

        finalGraph = nx.Graph()

        for subgraph in self.workingSubgraphsList:
            finalGraph = nx.union(finalGraph, subgraph)

        return finalGraph

    def connect_subgraphs(self):
        """

        Adds edges to the resultGraph to connect as many components of the final
        graph possible
        
        """

        connectSuccess = self.connect_graph_components_brute_force()

        while connectSuccess:
            connectSuccess = self.connect_graph_components_brute_force()

        # WARNING: The self.workingSubgraphsList at this point is different from
        # the copy self.resultingSubgraphsList made before

        connectSuccess = self.connect_graph_components_brute_force_2()

        while connectSuccess:
            connectSuccess = self.connect_graph_components_brute_force_2()

    def connect_graph_components_brute_force(self):
        """
        Adds edges to the resultGraph to connect all components that can be 
        connected, only one edge is added per component, to form a tree like 
        structure between the different components of the resultGraph
        
        Returns
        -------
        bool
            True if the addition of edges was possible in strict mode, False otherwise

        """

        generator_graph = [self.resultGraph.subgraph(c).copy() for c in nx.connected_components(self.resultGraph)]

        self.workingSubgraphsList = [x for x in generator_graph]

        if len(self.workingSubgraphsList) == 1:
            return False

        edgesToCheck = []
        edgesToCheckAdditionalInfo = []
        numzeros = 0

        for i in range(0, len(self.workingSubgraphsList)):

            nodesOfI = self.workingSubgraphsList[i].nodes()

            for j in range(i + 1, len(self.workingSubgraphsList)):
                nodesOfJ = self.workingSubgraphsList[j].nodes()

                # change the following lines to be compatible with networkx 2.0
                for k in nodesOfI.keys():

                    for l in nodesOfJ.keys():
                        # produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push
                        # this edge into possibleEdgeList """

                        # print 'Molecules (%d,%d)' % (nodesOfI[k],nodesOfJ[l])
                        # I assumed that the score matrix is symmetric. In the Graph part this
                        # does not seems to be true:

                        similarity = self.dbase.loose_mtx[nodesOfI[k]["ID"], nodesOfJ[l]["ID"]]

                        if similarity > 0.0:
                            edgesToCheck.append((nodesOfI[k]["ID"], nodesOfJ[l]["ID"], similarity))
                            edgesToCheckAdditionalInfo.append((nodesOfI[k]["ID"], nodesOfJ[l]["ID"], similarity, i, j))
                        else:
                            numzeros = numzeros + 1

        if len(edgesToCheck) > 0:

            sortedList = sorted(edgesToCheck, key=itemgetter(2), reverse=True)
            sortedListAdditionalInfo = sorted(edgesToCheckAdditionalInfo, key=itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            # self.edgeFile.write("\n" + str(edgeToAdd))
            edgeToAddAdditionalInfo = sortedListAdditionalInfo[0]
            self.edgesAddedInFirstTreePass.append(edgeToAdd)
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False)

            generator_graph = [self.resultGraph.subgraph(c).copy() for c in nx.connected_components(self.resultGraph)]
            self.workingSubgraphsList = [x for x in generator_graph]

            return True

        else:
            return False

    def connect_graph_components_brute_force_2(self):
        """
        Adds a second edge between each of the (former) components of the
        resultGraph to try to provide cycles between (former) components
        
        Returns
        -------
        bool
            True if the addition of edges was possible in loose mode, False otherwise

        """

        if len(self.resultingSubgraphsList) == 1:
            return False

        edgesToCheck = []

        for i in range(0, len(self.resultingSubgraphsList)):

            nodesOfI = self.resultingSubgraphsList[i].nodes()

            for j in range(i + 1, len(self.resultingSubgraphsList)):

                nodesOfJ = self.resultingSubgraphsList[j].nodes()

                for k in nodesOfI.keys():

                    for l in nodesOfJ.keys():

                        # produce an edge from nodesOfI[k] and nodesofJ[l] if
                        # nonzero weights push this edge into possibleEdgeList """

                        # print 'Molecules (%d,%d)' % (nodesOfI[k],nodesOfJ[l])
                        # I assumed that the score matrix is symmetric. In the Graph part
                        # this does not seems to be true: <<<<<<<<<<<<<DEBUG>>>>>>>>>>>>>>>
                        similarity = self.dbase.loose_mtx[nodesOfI[k]["ID"], nodesOfJ[l]["ID"]]

                        if similarity > 0.0:
                            edgesToCheck.append((nodesOfI[k]["ID"], nodesOfJ[l]["ID"], similarity))

        finalEdgesToCheck = [edge for edge in edgesToCheck if edge not in self.edgesAddedInFirstTreePass]

        if len(finalEdgesToCheck) > 0:

            sortedList = sorted(finalEdgesToCheck, key=itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]

            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False)
            self.copyResultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False)

            generator_graph = [self.copyResultGraph.subgraph(c).copy() for c in nx.connected_components(self.copyResultGraph)]
            self.resultingSubgraphsList = [x for x in generator_graph]

            return True

        else:
            return False

    def get_graph(self):
        """

        Returns the final generated NetworkX graph

        """

        return self.resultGraph

    def generate_depictions(self):

        def max_dist_mol(mol):

            max_dist = 0.0
            conf = mol.GetConformer()

            for i in range(0, conf.GetNumAtoms()):

                crdi = np.array([conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z])

                for j in range(i + 1, conf.GetNumAtoms()):
                    crdj = np.array([conf.GetAtomPosition(j).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(j).z])
                    dist = np.linalg.norm(crdi - crdj)

                    if dist > max_dist:
                        max_dist = dist

            return max_dist

        directory_name = tempfile.mkdtemp()

        temp_graph = self.resultGraph.copy()

        if nx.number_of_nodes(temp_graph) <= self.max_images:
            # Draw.DrawingOptions.atomLabelFontSize=30
            # Draw.DrawingOptions.dotsPerAngstrom=100

            for n in temp_graph:

                id_mol = temp_graph.nodes[n]['ID']
                mol = self.dbase[id_mol].getMolecule()
                max_dist = max_dist_mol(mol)

                if max_dist < self.max_mol_size:
                    fname = os.path.join(directory_name, self.dbase[id_mol].getName() + ".png")
                    # 1, modify here to calculate the 2D structure for ligands cannot remove Hydrogens by rdkit
                    # 2, change the graph size to get better resolution
                    try:
                        mol = AllChem.RemoveHs(mol)
                        AllChem.Compute2DCoords(mol)
                        from rdkit.Chem.Draw.MolDrawing import DrawingOptions
                        DrawingOptions.bondLineWidth = 2.5
                        Draw.MolToFile(mol, fname, size=(200, 200), kekulize=False, fitimage=True, imageType='png',
                                       options=DrawingOptions)
                    except:
                        ###### need to ask RDKit to fix this if possible, see the code
                        # issue tracker for more details######
                        logging.info(
                            "Error attempting to remove hydrogens for molecule %s using RDKit. RDKit cannot kekulize the molecule" %
                            self.dbase[id_mol].getName())
                        AllChem.Compute2DCoords(mol)
                        from rdkit.Chem.Draw.MolDrawing import DrawingOptions
                        DrawingOptions.bondLineWidth = 2.5
                        Draw.MolToFile(mol, fname, size=(200, 200), kekulize=False, fitimage=True, imageType='png',
                                       options=DrawingOptions)
                    temp_graph.nodes[n]['image'] = fname
                    # self.resultGraph.nodes[n]['label'] = ''
                    temp_graph.nodes[n]['labelloc'] = 't'
                    temp_graph.nodes[n]['penwidth'] = 2.5
                    # self.resultGraph.node[n]['xlabel'] =  self.resultGraph.nodes[n]['ID']
        for u, v, d in temp_graph.edges(data=True):
            if d['strict_flag'] == True:
                temp_graph[u][v]['color'] = 'blue'
                temp_graph[u][v]['penwidth'] = 2.5
            else:
                temp_graph[u][v]['color'] = 'red'
                temp_graph[u][v]['penwidth'] = 2.5
            if self.edge_labels:
                temp_graph[u][v]['label'] = round(d['similarity'],2)

        nx.nx_agraph.write_dot(temp_graph, self.dbase.options.name + '_tmp.dot')

        cmd = 'dot -Tpng ' + self.dbase.options.name + '_tmp.dot -o ' + self.dbase.options.name + '.png'

        os.system(cmd)
        cmd = 'dot -Teps ' + self.dbase.options.name + '_tmp.dot -o ' + self.dbase.options.name + '.eps'

        os.system(cmd)
        cmd = 'dot -Tpdf ' + self.dbase.options.name + '_tmp.dot -o ' + self.dbase.options.name + '.pdf'

        os.system(cmd)
        os.remove(self.dbase.options.name + '_tmp.dot')
        shutil.rmtree(directory_name, ignore_errors=True)

    # The function to output the score and connectivity txt file

    def layout_info(self):
        # pass the lead compound index if the radial option is on and generate the
        # morph type of output required by FESetup
        if self.lead_index is not None:
            morph_txt = open(self.dbase.options.name + "_morph.txt", "w")
            morph_data = "morph_pairs = "
        with open(self.dbase.options.name + "_score_with_connection.txt", "w") as info_txt:
            all_key_id = self.dbase.dic_mapping.keys()
            data = ["%-10s,%-10s,%-25s,%-25s,%-15s,%-15s,%-15s,%-10s\n" % (
            "Index_1", "Index_2", "Filename_1", "Filename_2", "Str_sim", "Eff_sim", "Loose_sim", "Connect")]
            for i in range(len(all_key_id) - 1):
                for j in range(i + 1, len(all_key_id)):
                    morph_string = None
                    connected = False
                    similarity=0
                    try:
                        edgedata=[d for (u,v,d) in self.resultGraph.edges(data=True) if ((u==i and v==j) or (u==j and v==i))]
                        similarity = edgedata[0]['similarity']
                        connected = True
                    except IndexError:
                        pass
                    Filename_i = self.dbase.dic_mapping[i]
                    Filename_j = self.dbase.dic_mapping[j]
                    MCmap = self.dbase.get_MCSmap(i,j)
                    mapString=""
                    if MCmap is not None:
                        mapString = MCmap
                    # print "Check the filename", Filename_i, Filename_j
                    strict_similarity = self.dbase.strict_mtx[i, j]
                    loose_similarity = self.dbase.loose_mtx[i, j]
                    true_strict_similarity = self.dbase.true_strict_mtx[i, j]
                    if connected:
                        new_line = "%-10s,%-10s,%-25s,%-25s,%-15.5f,%-15.5f,%-15.5f,%-10s,%s\n" % (
                        i, j, Filename_i, Filename_j, true_strict_similarity, strict_similarity, loose_similarity, "Yes",mapString)
                        # generate the morph type, and pick the start ligand based on the similarity
                        if self.lead_index is not None:
                            morph_i = Filename_i.split(".")[0]
                            morph_j = Filename_j.split(".")[0]
                            if i == self.lead_index:
                                morph_string = "%s > %s, " % (morph_i, morph_j)
                            elif j == self.lead_index:
                                morph_string = "%s > %s, " % (morph_j, morph_i)
                            else:
                                # compare i and j with the lead compound, and
                                # pick the one with the higher similarity as the start ligand
                                similarity_i = self.dbase.strict_mtx[self.lead_index, i]
                                similarity_j = self.dbase.strict_mtx[self.lead_index, j]
                                if similarity_i > similarity_j:
                                    morph_string = "%s > %s, " % (morph_i, morph_j)
                                else:
                                    morph_string = "%s > %s, " % (morph_j, morph_i)
                            morph_data += morph_string
                    else:
                        new_line = "%-10s,%-10s,%-25s,%-25s,%-15.5f,%-15.5f,%-15.5f,%-10s,%s\n" % (
                        i, j, Filename_i, Filename_j, true_strict_similarity, strict_similarity, loose_similarity, "No",mapString)
                    data.append(new_line)
            info_txt.writelines(data)
            if self.lead_index is not None:
                morph_txt.write(morph_data)
                

    def write_graph(self, output_no_images, output_no_graph):
        """

        This function writes to a file the final generated NetworkX graph as 
        .dot and the .ps files. The mapping between molecule IDs and compounds
        name is saved as text file


        """

        try:
            self.dbase.write_dic()
            self.layout_info()
        except Exception as e:
            traceback.print_exc()
            raise IOError("%s: %s.txt" % (str(e), self.dbase.options.name))

        try:
            if not output_no_images:
                self.generate_depictions()
            if not output_no_graph:
                nx.nx_agraph.write_dot(self.resultGraph, self.dbase.options.name + '.dot')
        except Exception as e:
            traceback.print_exc()
            raise IOError('Problems during the file generation: %s' % str(e))

        logging.info(30 * '-')

        log = 'The following files have been generated:'
        if not output_no_graph:
            log += f'\n{self.dbase.options.name}.dot\tGraph file'
        if not output_no_images:
            log += f'\n{self.dbase.options.name}.png\tPng file'
        log += f'\n{self.dbase.options.name}.txt\tMapping Text file'

        logging.info(30 * '-')

        return

    ###### Still in developing stage ######

    def draw(self):
        """
        This function plots the NetworkX graph by using Matplotlib
        
        """

        logging.info('\nDrawing....')

        if nx.number_of_nodes(self.resultGraph) > self.max_nodes:
            logging.info('The number of generated graph nodes %d exceed the max number of drawable nodes %s' % (
            nx.number_of_nodes(self.resultGraph), self.max_nodes))
            return

        def max_dist_mol(mol):

            max_dist = 0.0
            conf = mol.GetConformer()

            for i in range(0, conf.GetNumAtoms()):

                crdi = np.array([conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z])

                for j in range(i + 1, conf.GetNumAtoms()):
                    crdj = np.array([conf.GetAtomPosition(j).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(j).z])
                    dist = np.linalg.norm(crdi - crdj)

                    if dist > max_dist:
                        max_dist = dist

            return max_dist

        # Determine the screen resolution by using dxpyinfo and removing massive qt dependency
        command = ('xdpyinfo | grep dimensions')
        p = subprocess.run(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
        width = int(p.stdout.split()[1].split(b'x')[0])
        height = int(p.stdout.split()[1].split(b'x')[1])

        # Canvas scale factor 
        scale_canvas = 0.75

        # Canvas resolution
        max_canvas_size = (int(width * scale_canvas), int(height * scale_canvas))

        fig = plt.figure(1, facecolor='white')

        fig.set_dpi(100)

        fig.set_size_inches(max_canvas_size[0] / fig.get_dpi(), max_canvas_size[1] / fig.get_dpi(), forward=True)

        ax = plt.subplot(111)
        plt.axis('off')

        pos = nx.nx_agraph.graphviz_layout(self.resultGraph, prog="neato")

        strict_edges = [(u, v) for (u, v, d) in self.resultGraph.edges(data=True) if d['strict_flag'] == True]
        loose_edges = [(u, v) for (u, v, d) in self.resultGraph.edges(data=True) if d['strict_flag'] == False]

        node_labels = dict([(u, d['ID']) for u, d in self.resultGraph.nodes(data=True)])

        # Draw nodes
        nx.draw_networkx_nodes(self.resultGraph, pos, node_size=500, node_color='r')
        # Draw node labels
        nx.draw_networkx_labels(self.resultGraph, pos, labels=node_labels, font_size=10)

        if self.edge_labels:
            edge_weight_strict = dict([((u, v,), d['similarity']) for u, v, d in self.resultGraph.edges(data=True) if
                                       d['strict_flag'] == True])
            edge_weight_loose = dict([((u, v,), d['similarity']) for u, v, d in self.resultGraph.edges(data=True) if
                                      d['strict_flag'] == False])

            for key in edge_weight_strict:
                edge_weight_strict[key] = round(edge_weight_strict[key], 2)

            for key in edge_weight_loose:
                edge_weight_loose[key] = round(edge_weight_loose[key], 2)

            # edge strict
            nx.draw_networkx_edge_labels(self.resultGraph, pos, edge_labels=edge_weight_strict, font_color='g')
            # edge loose
            nx.draw_networkx_edge_labels(self.resultGraph, pos, edge_labels=edge_weight_loose, font_color='r')

        # edges strict
        nx.draw_networkx_edges(self.resultGraph, pos, edgelist=strict_edges, edge_color='g')
        # edges loose
        nx.draw_networkx_edges(self.resultGraph, pos, edgelist=loose_edges, edge_color='r')

        if nx.number_of_nodes(self.resultGraph) <= self.max_images:

            trans = ax.transData.transform
            trans2 = fig.transFigure.inverted().transform

            cut = 1.0

            frame = 10
            xmax = cut * max(xx for xx, yy in pos.values()) + frame
            ymax = cut * max(yy for xx, yy in pos.values()) + frame

            xmin = cut * min(xx for xx, yy in pos.values()) - frame
            ymin = cut * min(yy for xx, yy in pos.values()) - frame

            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            h = 20
            w = 20

            mol_size = (200, 200)

            for each_node in self.resultGraph:

                id_mol = self.resultGraph.nodes[each_node]['ID']
                # skip remove Hs by rdkit if Hs cannot be removed
                try:
                    mol = AllChem.RemoveHs(self.dbase[id_mol].getMolecule())
                except:
                    ###### need to ask RDKit to fix this if possible, see the code
                    # issue tracker for more details######
                    mol = self.dbase[id_mol].getMolecule()
                    logging.info(
                        "Error attempting to remove hydrogens for molecule %s using RDKit. RDKit cannot kekulize the molecule" %
                        self.dbase[id_mol].getName())

                # max_dist = max_dist_mol(mol)
                # if max_dist > 7.0:
                #     continue

                AllChem.Compute2DCoords(mol)
                # add try exception for cases cannot be draw
                try:
                    img_mol = Draw.MolToImage(mol, mol_size, kekulize=False)
                except Exception as ex:
                    img_mol = None
                    logging.exception(
                        "This mol cannot be draw using the RDKit Draw function, need to check for more details...")

                xx, yy = trans(pos[each_node])
                xa, ya = trans2((xx, yy))

                nodesize_1 = (300.0 / (h * 100))
                nodesize_2 = (300.0 / (w * 100))

                p2_2 = nodesize_2 / 2
                p2_1 = nodesize_1 / 2

                a = plt.axes([xa - p2_2, ya - p2_1, nodesize_2, nodesize_1])
                # self.resultGraph.nodes[id_mol]['image'] = img_mol
                # a.imshow(self.resultGraph.node[each_node]['image'])
                a.imshow(img_mol)
                a.axis('off')

        # plt.savefig('graph.png', facecolor=fig.get_facecolor())
        # print 'Graph .png file has been generated...'

        plt.show()

        return
