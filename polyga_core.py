#!/usr/bin/env python

"""

PolyGA - version 1.2b - Last modified: 31-DEC-2013

This file contains the core functions of the genetic algorithm (GA) 
implemented in the PolyGA program.

Copyright (C) 2013 Michael Mooney

This file is part of PolyGA.

PolyGA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PolyGA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the file COPYING). 
If not, see <http://www.gnu.org/licenses/>.

Send bug reports, requests, etc. to mooneymi@ohsu.edu

"""

import sys
import os
import random
import math
import time
import copy
import tables
import networkx as nx
import rpy2.robjects as rob
import polyga_utils

class PolyGAInd():
	"""
	
	"""
	
	def __init__(self, connected=False, level='feature', min_group_size=2, max_group_size=2):
		"""
		
		"""
		
		## Define object attributes
		if connected == True:
			self.connected = True
		else:
			self.connected = False
			
		if level == 'gene':
			self.level = 'gene'
		else:
			self.level = 'feature'
			
		if self.connected:
			self.max_size = 2*max_group_size - 1
		else:
			self.max_size = max_group_size
		
		self.min_group_size = min_group_size
		self.max_group_size = max_group_size
		self.group_id = None
		self.fitness = 999.0
		
		## Create empty NetworkX graph
		self.graph = nx.Graph()
			
		return None
	
	
	def __str__(self):
		"""
		
		"""
		
		return str(self.group_id)+'\n'+str(self.connected)+'\n'+str(self.level)+'\n'+str(self.min_group_size)+'\n'+str(self.max_group_size)+'\n'+str(self.max_size)+'\n'+str(self.graph.nodes())+'\n'+str(self.graph.edges())+'\n'+str(nx.get_node_attributes(self.graph, 'features'))+'\n'+str(self.get_size())+'\n'+str(self.fitness)
	
	
	def new(self):
		"""
		
		"""
		
		new_ind = PolyGAInd(self.connected, self.level, self.min_group_size, self.max_group_size)
		
		return new_ind
	
	
	def copy(self):
		"""
		
		"""
		
		new_ind = PolyGAInd(self.connected, self.level, self.min_group_size, self.max_group_size)
		new_ind.graph = self.graph.copy()
		new_ind.fitness = self.fitness
		
		return new_ind
	
	
	def update(self, ind):
		"""
		
		"""
		
		self.connected = ind.connected
		self.level = ind.level
		self.min_group_size = ind.min_group_size
		self.max_group_size = ind.max_group_size
		self.max_size = ind.max_size
		self.group_id = ind.group_id
		self.graph = ind.graph.copy()
		self.fitness = ind.fitness
	
		return self
	
	
	def get_size(self):
		"""
		
		"""
		if self.level == 'feature':
			node_attrs = nx.get_node_attributes(self.graph, 'features')
			size = 0
			for g, f in node_attrs.items():
				size = size + len(f)
			return size
		else:
			return self.graph.number_of_nodes()
	
	def get_genes(self):
		"""
		
		"""
		return self.graph.nodes()
	
	
	def get_features(self):
		"""
		
		"""
		
		return [fid for n, f in nx.get_node_attributes(self.graph, 'features').items() for fid in f.keys()]
	
	
	def get_node_properties(self):
		"""
		
		"""
		nodes = []
		weights = []
		for n, attrs in self.graph.nodes(data=True):
			g = n
			gs = attrs['weight']
			if self.level == 'feature':
				for f, fs in attrs['features'].items():
					nodes.append([g, f])
					weights.append([gs, fs])
			else:
				nodes.append(g)
				weights.append(gs)
		
		return nodes, weights
	
	
	def get_neighbors(self, G, exclude_origins=False, steps=1):
		"""
		
		"""
		
		neighbors = {}
		nodes = self.graph.nodes()
		if self.level == 'gene':
			exclude_origins = True
		
		for n in nodes:
			for k, v in G[n].items():
				if (n,k) in neighbors:
					if v['weight'] > neighbors[(n,k)]:
						neighbors[(n,k)] = v['weight']
				else:
					neighbors[(n,k)] = v['weight']
		
		if exclude_origins:
			neighbors = dict([((n,k), w) for (n,k), w in neighbors.items() if n not in nodes])
		
		while steps > 1:
			new_neighbors = {}
			for n, w in neighbors.items():
				for k, v in G[n].items():
					if (n,k) in new_neighbors:
						if v['weight'] > new_neighbors[(n,k)]:
							new_neighbors[(n,k)] = v['weight']
					else:
						new_neighbors[(n,k)] = v['weight']
				
			if exclude_origins:
				new_neighbors = dict([((n,k), w) for (n,k), w in new_neighbors.items() if n not in nodes])
			
			for (n,k), w in new_neighbors.items():
				if (n,k) in neighbors:
					if w > neighbors[(n,k)]:
						neighbors[(n,k)] = w
				else:
					neighbors[(n,k)] = w
			
			steps = steps - 1
		
		return neighbors
	
	
	def get_neighborhood(self, G, steps=1):
		"""
		
		"""
		
		nodes = self.graph.nodes()
		neighbors = list(nodes)
		
		for n in nodes:
			neighbors.extend([k for k,v in G[n].items()])
		neighbors = list(set(neighbors))
		
		while steps > 1:
			new_neighbors = []
			for n in neighbors:
				new_neighbors.extend([k for k,v in G[n].items()])
				
			neighbors.extend(new_neighbors)
			neighbors = list(set(neighbors))
			steps = steps - 1
		
		neighbors = list(set(neighbors))
		neighborhood = G.subgraph(neighbors).copy()
		
		return neighborhood
	
	
	def is_connected(self):
		"""
		
		"""
		
		if self.connected == False:
			return False
		else:
			try:
				c = nx.is_connected(self.graph)
			except nx.exception.NetworkXPointlessConcept:
				return None
			else:
				return c
	
	
	def largest_connected(self):
		"""
		
		"""
		
		if self.connected == False:
			return None
		else:
			components = nx.connected_components(self.graph)
			if len(components) > 0:
				largest_connected_subgraph = self.graph.subgraph(components[0]).copy()
				self.graph = largest_connected_subgraph
				self.fitness = 999.0
			return None
	
	
	def intersection(self, ind2):
		"""
		
		"""
		#new_ind1 = self.new()
		#new_ind2 = ind2.new()
		intersection_nodes = []
		
		for n1 in self.graph.nodes():
			for n2 in ind2.graph.nodes():
				if self.graph.node[n1] == ind2.graph.node[n2]:
					intersection_nodes.append(n1)
		
		#new_ind1.graph = self.graph.subgraph(intersection_nodes).copy()
		#new_ind2.graph = ind2.graph.subgraph(intersection_nodes).copy()
		
		#return new_ind1, new_ind2
		return intersection_nodes
	
	
	def difference(self, ind2):
		"""
		********NOT FINISHED********
		"""
		#new_ind = self.new()
		nodes = self.graph.nodes()
		intersection_nodes = self.intersection(ind2)
		difference_nodes = list(set(nodes) - set(intersection_nodes))
		
		return difference_nodes
	
	
	def create_connected(self, G, data, gene_weight, int_weight):
		"""
		
		"""
		
		self.connected = True
		
		## Get the HDF5 file handle and data tables
		fh, node_table, edge_table, feature_table = data.get_tables()
		
		## Load the gene-gene interaction network from the HDF5 file
		#G = data.get_graph()
		
		## Select the first gene
		#b1, ts1 = data.gene_scores_wheel
		#gene1 = data.genes[spin_roulette_wheel(b1, ts1)]
		gene1 = self.get_genes()[0]
		
		if self.level == 'feature':
			## Select a feature mapped to the first gene
			features1_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""(gene == gene1)""")]
			b2, ts2 = create_roulette_wheel([fs for f, fs in features1_scores])
			idx2 = spin_roulette_wheel(b2, ts2)
			feature1 = features1_scores[idx2][0]
			fscore1 = features1_scores[idx2][1]
			
			## Append the gene and feature to the graph
			self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={feature1:fscore1})
		else:
			self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={})
		
		rand_group_size = random.randint(self.min_group_size,  self.max_group_size)
		while self.get_size() < rand_group_size:
			## Select the next gene from the group's neighbors
			neighbor_gene_scores = []
			group_genes = self.get_genes()
			for cur_gene in group_genes:
				neighbor_genes = [(ng, G.edge[cur_gene][ng]['weight']) for ng in G.neighbors(cur_gene)]
				
				## Get weighted scores
				for neighbor_gene, int_score in neighbor_genes:
					gene_score = G.node[neighbor_gene]['weight']
					weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
					neighbor_gene_scores.append((cur_gene, neighbor_gene, weighted_score))
			
			#print neighbor_gene_scores
			if self.level == 'feature':
				## Check that non-duplicate neighbor features are available
				group_features = self.get_features()
				all_neighbor_features = []
				for cg, ng, ws in neighbor_gene_scores:
					neighbor_features = [row['feature'] for row in feature_table.where("""gene == ng""")]
					neighbor_features = [f for f in neighbor_features if f not in group_features]
					all_neighbor_features.extend(neighbor_features)
				
				if len(all_neighbor_features) == 0:
					if self.get_size() > self.min_group_size:
						break
					else:
						## Start over by selecting new gene
						self.graph.clear()
						gene1 = data.genes[spin_roulette_wheel(b1, ts1)]
						
						## Select a feature mapped to the first gene
						features1_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""(gene == gene1)""")]
						b2, ts2 = create_roulette_wheel([fs for f, fs in features1_scores])
						idx2 = spin_roulette_wheel(b2, ts2)
						feature1 = features1_scores[idx2][0]
						fscore1 = features1_scores[idx2][1]
						
						## Append the gene and feature to the graph
						self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={feature1:fscore1})
						continue
				else:
					## Select the gene to add based on weighted scores
					new_gene = None
					new_feature = None
					b3, ts3 = create_roulette_wheel([ws for cg, ng, ws in neighbor_gene_scores])
					while new_gene == None or new_feature == None:
						idx3 = spin_roulette_wheel(b3, ts3)
						group_gene = neighbor_gene_scores[idx3][0]
						new_gene = neighbor_gene_scores[idx3][1]
						
						## Get all features mapped to selected gene, and randomly select the feature to add
						neighbor_feature_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""gene == new_gene""")]
						neighbor_feature_scores = [(f, fs) for f, fs in neighbor_feature_scores if f not in group_features]
						if len(neighbor_feature_scores) > 0:
							b4, ts4 = create_roulette_wheel([fs for f, fs in neighbor_feature_scores])
							idx4 = spin_roulette_wheel(b4, ts4)
							new_feature = neighbor_feature_scores[idx4][0]
							new_fscore = neighbor_feature_scores[idx4][1]
					
					## Append the gene and feature to the graph
					if new_gene in self.get_genes():
						self.graph.node[new_gene]['features'].update({new_feature:new_fscore})
						self.graph.add_edge(group_gene, new_gene, weight=G.edge[group_gene][new_gene]['weight'])
					else:
						self.graph.add_node(new_gene, weight=G.node[new_gene]['weight'], features={new_feature:new_fscore})
						self.graph.add_edge(group_gene, new_gene, weight=G.edge[group_gene][new_gene]['weight'])
			else:
				## Remove duplicate genes from the neighbors
				group_genes = self.get_genes()
				neighbor_gene_scores = [(cg, ng, ws) for cg, ng, ws in neighbor_gene_scores if ng not in group_genes]
				
				if len(neighbor_gene_scores) == 0:
					if self.get_size() > self.min_group_size:
						break
					else:
						## Start over by selecting new gene
						self.graph.clear()
						gene1 = data.genes[spin_roulette_wheel(b1, ts1)]
						self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={})
						continue
				else:
					## Select gene to add based on weighted scores
					new_gene = None
					b3, ts3 = create_roulette_wheel([ws for cg, ng, ws in neighbor_gene_scores])
					idx3 = spin_roulette_wheel(b3, ts3)
					group_gene = neighbor_gene_scores[idx3][0]
					new_gene = neighbor_gene_scores[idx3][1]
					
					## Append the gene to the graph
					self.graph.add_node(new_gene, weight=G.node[new_gene]['weight'], features={})
					self.graph.add_edge(group_gene, new_gene, weight=G.edge[group_gene][new_gene]['weight'])
		del G
		del data
		fh.close()
		
		return self
	
	
	def create(self, G, data, gene_weight=0.5, int_weight=0.5):
		"""
		
		"""
		
		if (gene_weight + int_weight != 1) or (gene_weight < 0) or (int_weight < 0):
			gene_weight = 0.5
			int_weight = 0.5
		
		if self.connected:
			return self.create_connected(data, gene_weight, int_weight)
		else:
			## Get the HDF5 file handle and data tables
			fh, node_table, edge_table, feature_table = data.get_tables()
			
			## Load the gene-gene interaction network from the HDF5 file
			#G = data.get_graph(edges=False)
		
			## Select the first gene
			#b1, ts1 = data.gene_scores_wheel
			#gene1 = data.genes[spin_roulette_wheel(b1, ts1)]
			gene1 = self.get_genes()[0]
				
			if self.level == 'feature':
				## Select a feature mapped to the first gene
				features1_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""(gene == gene1)""")]
				b2, ts2 = create_roulette_wheel([fs for f, fs in features1_scores])
				idx2 = spin_roulette_wheel(b2, ts2)
				feature1 = features1_scores[idx2][0]
				fscore1 = features1_scores[idx2][1]
				
				## Append the gene and feature to the graph
				self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={feature1:fscore1})
			else:
				self.graph.add_node(gene1, weight=G.node[gene1]['weight'], features={})
		
			rand_group_size = random.randint(self.min_group_size,  self.max_group_size)
			for j in range(1, rand_group_size):
				## Select the next gene
				if self.level == 'feature':
					group_features = self.get_features()
					new_gene = data.genes[spin_roulette_wheel(b1, ts1)]
					new_features_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""(gene == new_gene)""")]
					new_features_scores = [(f, fs) for f, fs in new_features_scores if f not in group_features]
					
					## Check that the selected new gene has non-duplicate features available
					while (len(new_features_scores) == 0): 
						new_gene = data.genes[spin_roulette_wheel(b1, ts1)]
						new_features_scores = [(row['feature'], row['fscore']) for row in feature_table.where("""(gene == new_gene)""")]
						new_features_scores = [(f, fs) for f, fs in new_features_scores if f not in group_features]
					
					## Select the next feature
					b3, ts3 = create_roulette_wheel([fs for f, fs in new_features_scores])
					idx3 = spin_roulette_wheel(b3, ts3)
					new_feature = new_features_scores[idx3][0]
					new_fscore = new_features_scores[idx3][1]
					
					## Append the gene-feature pair to the group list
					self.graph.add_node(new_gene, weight=G.node[new_gene]['weight'], features={new_feature:new_fscore})
				else:
					group_genes = self.get_genes()
					while True:
						new_gene = data.genes[spin_roulette_wheel(b1, ts1)]
						if not new_gene in group_genes:
							self.graph.add_node(new_gene, weight=G.node[new_gene]['weight'], features={})
							break
			del G
			del data
			fh.close()
			return self
	
	
	def crossover(self, parent2, cross_type, cross_prob):
		"""
		
		"""
		
		# *** DEBUGGING / TESTING ***
		#print self
		#print parent2
		#print "crossing over"
		
		## Check if the two feature groups are identical, if so enter into the next generation unchanged
		exact_match = False
		if self.level == 'feature':
			if set(self.get_features()) == set(parent2.get_features()):
				exact_match = True
		else:
			if set(self.get_genes()) == set(parent2.get_genes()):
				exact_match = True
		if exact_match:
			## Insert the two parents into the population unchanged
			child1 = self.copy()
			child2 = parent2.copy()
		else:
			if cross_type == 'uniform':
				## Swap only features that differ between the two parents, so first identify the matches and mismatches
				child1 = self.new()
				child2 = parent2.new()
				intersection_genes = self.intersection(parent2)
				child1.graph = self.graph.subgraph(intersection_nodes).copy()
				child2.graph = parent2.graph.subgraph(intersection_nodes).copy()
				parent1_genes = self.difference(parent2)
				parent2_genes = parent2.difference(self)
				
				## Perform the crossover operations with the specified frequency
				if len(parent1_genes) < len(parent2_genes):
					for idx, gene in enumerate(parent1_genes):
						rand = random.random()
						if rand < cross_prob:
							rand = random.randint(0, (len(parent2_genes)-1))
							attrs1 = self.graph.node[gene]
							child2.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
							attrs2 = parent2.graph.node[parent2_genes[rand]]
							child1.graph.add_node(parent2_genes[rand], weight=attrs2['weight'], features=attrs2['features'])
							parent2_genes.pop(rand)
						else:
							attrs1 = self.graph.node[gene]
							child1.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
					## Add remainder of larger group to appropriate offspring
					for gene in parent2_genes:
						attrs2 = parent2.graph.node[gene]
						child2.graph.add_node(gene, weight=attrs2['weight'], features=attrs2['features'])
				elif len(parent2_genes) < len(parent1_genes):
					for idx, gene in enumerate(parent2_genes):
						rand = random.random()
						if rand < cross_prob:
							rand = random.randint(0, (len(parent1_genes)-1))
							attrs2 = parent2.graph.node[gene]
							child1.graph.add_node(parent2_genes[rand], weight=attrs2['weight'], features=attrs2['features'])
							attrs1 = self.graph.node[parent1_genes[rand]]
							child2.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
							parent1_genes.pop(rand)
						else:
							attrs2 = parent2.graph.node[gene]
							child2.graph.add_node(gene, weight=attrs2['weight'], features=attrs2['features'])
					## Add remainder of larger group to appropriate offspring
					for gene in parent1_genes:
						attrs1 = self.graph.node[gene]
						child1.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
				else:
					for idx, gene in enumerate(parent1_genes[:-1]):
						rand = random.random()
						if rand < cross_prob:
							#print "Crossing over 3!\n"
							rand = random.randint(0, (len(parent2_genes)-1))
							attrs1 = self.graph.node[gene]
							child2.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
							attrs2 = parent2.graph.node[parent2_genes[rand]]
							child1.graph.add_node(parent2_genes[rand], weight=attrs2['weight'], features=attrs2['features'])
							parent2_genes.pop(rand)
						else:
							attrs1 = self.graph.node[gene]
							child1.graph.add_node(gene, weight=attrs1['weight'], features=attrs1['features'])
					## Add remainders to the appropriate offspring
					attrs1 = self.graph.node[parent1_genes[-1]]
					child1.genes.append(parent1_genes[-1], weight=attrs1['weight'], features=attrs1['features'])
					for gene in parent2_genes:
						attrs2 = parent2.graph.node[gene]
						child2.graph.add_node(gene, weight=attrs2['weight'], features=attrs2['features'])
			elif cross_type == 'single-point':
				child1 = self.copy()
				child2 = parent2.copy()
			
			## Check the group sizes
			
			#if size == 0
			
			
			#if size > max_group_size
		
		# *** DEBUGGING / TESTING ***
		#print self
		#print parent2
		#print child1
		#print child2
		#print len(child1.features)
		#print len(child2.features)
		#print len(child1.genes)
		#print len(child2.genes)
		
		return child1, child2
	
	
	def mutate_remove_connected(self, data, gene_weight=0.5, int_weight=0.5):
		"""
		****** CONSOLIDATE WITH OTHER REMOVE FUNCTION *******
		"""
		
		# *** DEBUGGING / TESTING *** 
		#print "Removing connected..."
		#print self
		
		self.connected = True
		
		## Get the HDF5 file handle and data tables
		fh, gene_table, int_table, feature_table = data.get_tables()
			
		if len(self.graph.edges()) > 1:
			## Randomly select feature to remove from group
			gene_scores = [self.graph.node[n]['weight'] for n in self.graph.nodes()]
						
			## Create a roulette wheel with gene scores and iterate through, the last to be selected will be removed from the group
			new_genes = [None]*len(self.get_genes())
			i = 0
			while i < (len(self.genes)-1):
				b1, ts1 = create_roulette_wheel(gene_scores)
				keep_index = spin_roulette_wheel(b1, ts1)
				if self.genes[keep_index] != None:
					new_genes[keep_index] = self.genes[keep_index]
					self.genes[keep_index] = None
					i = i + 1
			new_genes = [gene for gene in new_genes if gene != None]
			removed_gene = [gene for gene in self.genes if gene != None][0]
			
			## Remove selected gene from graph and get the largest connected subgraph
			self.graph.remove_node(removed_gene)
			
			if self.connected == True:
				self.largest_connected()
				
				## Check that the mutated group is not empty. If it is empty create a new group.
				if len(self.graph.edges()) == 0:
					# *** DEBUGGING / TESTING ***
					#print self
					#print "create new group"
					
					new_group = self.new()
					new_group.create(data, gene_weight, int_weight)
					new_group.group_id = self.group_id
					self.update(new_group)
		
		fh.close()
		return self
	
	
	def mutate_add_connected(self, G, data, gene_weight=0.5, int_weight=0.5):
		"""
		
		"""
		
		# *** DEBUGGING / TESTING *** 
		#print "Adding connected..."
		#print self
		
		self.connected = True
		
		## Get the HDF5 file handle and data tables
		fh, gene_table, int_table, feature_table = data.get_tables()
		
		## Select the next gene from the group's neighbors
		neighbor_genes = self.get_neighbors(G)
		neighbor_gene_scores = []
			
		## Get neighbor gene scores
		for (cur_gene, neighbor_gene), int_score in neighbor_genes.items():
			gene_score = G.node[neighbor_gene]['weight']
			weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
			neighbor_gene_scores.append((cur_gene, neighbor_gene, weighted_score))
			
		if self.level == 'feature':
			## Check that non-duplicate neighbor features are available
			all_neighbor_features = []
			group_features = self.get_features()
			for cg, ng, ws in neighbor_gene_scores:
				neighbor_features = [row['feature'] for row in feature_table.where("""gene == ng""")]
				neighbor_features = [f for f in neighbor_features if f not in group_features]
				all_neighbor_features.extend(neighbor_features)
			
			if len(all_neighbor_features) == 0:
				fh.close()
				return self
		
		if len(neighbor_gene_scores) > 0:
			## Select gene to add based on weighted scores, checking that a duplicate is not added
			new_gene = None
			b3, ts3 = create_roulette_wheel([ws for cg, ng, ws in neighbor_gene_scores])
			if self.level == 'feature':
				new_feature = None
			else:
				new_feature = 1
			while (new_gene == None) or (new_feature == None):
				neighbor_index = spin_roulette_wheel(b3, ts3)
				group_gene = neighbor_gene_scores[neighbor_index][0]
				new_gene = neighbor_gene_scores[neighbor_index][1]
				edge_weight = neighbor_genes[(group_gene, new_gene)]
				
				if self.level == 'feature':
					## Get all features mapped to selected gene, and randomly select the feature to add
					neighbor_feature_scores = [(row['feature'], row['score']) for row in feature_table.where("""gene == new_gene""")]
					neighbor_feature_scores = [(f, s) for f, s in neighbor_feature_scores if f not in group_features]
					if len(neighbor_feature_scores) > 0:
						b4, ts4 = create_roulette_wheel([s for f, s in neighbor_feature_scores])
						idx = spin_roulette_wheel(b4, ts4)
						new_feature = neighbor_feature_scores[idx][0]
						new_fscore = neighbor_feature_scores[idx][1]
			
			## Add new gene and/or feature to graph
			if self.level == 'feature':
				if new_gene in self.get_genes():
					self.graph.node[new_gene]['features'].update({new_feature:new_fscore})
					self.graph.add_edge(group_gene, new_gene, weight=edge_weight)
					self.fitness = 999.0
				else:
					self.graph.add_node(new_gene, weight=G.node[new_gene]['weight'], features={new_feature:new_fscore})
					self.graph.add_edge(group_gene, new_gene, weight=edge_weight)
					self.fitness = 999.0
			else:
				# **** NOT FINISHED ****
				pass
		
		fh.close()
		return self
	
	
	def mutate_alter_connected(self, data, gene_weight=0.5, int_weight=0.5):
		"""
		
		"""
		
		# *** DEBUGGING / TESTING *** 
		#print "Altering connected..."
		#print self
		
		self.connected = True
		
		## Get the HDF5 file handle and data tables
		fh, gene_table, int_table, feature_table = data.get_tables()
		
		## Randomly select gene to alter
		gene_scores = []
		for cur_gene in self.genes:
			gene_scores.extend([row['score'] for row in gene_table.where("""gene == cur_gene""")])
		
		## Create a roulette wheel with gene scores and iterate through, the last to be selected will be removed from the group
		new_genes = [None]*len(self.genes)
		if self.level == 'feature':
			new_features = [None]*len(self.features)
		i = 0
		keep_indices = []
		while i < (len(self.genes)-1):
			b1, ts1 = create_roulette_wheel(gene_scores)
			keep_index = spin_roulette_wheel(b1, ts1)
			if keep_index not in keep_indices:
				new_genes[keep_index] = self.genes[keep_index]
				if self.level == 'feature':
					new_features[keep_index] = self.features[keep_index]
				keep_indices.append(keep_index)
				i = i + 1
		alter_index = next(idx for idx,val in enumerate(self.genes) if idx not in keep_indices)
		alter_gene = self.genes[alter_index]
		
		if self.level == 'feature':
			alter_feature = self.features[alter_index]
		
		## Select a gene pair to be altered
		alter_pairs = [pair for pair in self.gene_pairs if (pair[0] == alter_gene or pair[1] == alter_gene)]
		pair_index = random.randint(0, len(alter_pairs)-1)
		
		for idx, pair in enumerate(self.gene_pairs):
			if pair == alter_pairs[pair_index]:
				pair_index = idx
				break
		
		## Remove the selected pair
		new_gene_pairs = self.gene_pairs[:]
		new_gene_pairs.pop(pair_index)
		if self.level == 'feature':
			new_feature_pairs = self.feature_pairs[:]
			new_feature_pairs.pop(pair_index)
		
		if len(new_gene_pairs) > 0:
			## Update self.genes and self.features
			self.gene_pairs = new_gene_pairs
			if self.level == 'feature':
				self.feature_pairs = new_feature_pairs
				genes = [g for gene_pair in self.gene_pairs for g in gene_pair]
				features = [f for feature_pair in self.feature_pairs for f in feature_pair]
				keep_features = []
				keep_indices = []
				for i, f in enumerate(features):
					if f not in keep_features:
						keep_features.append(f)
						keep_indices.append(i)
				keep_genes = [genes[j] for j in keep_indices]
				self.genes = list(keep_genes)
				self.features = list(keep_features)
				self.size = len(self.features)
			else:
				genes = [g for gene_pair in self.gene_pairs for g in gene_pair]
				keep_genes = []
				for g in genes:
					if g not in keep_genes:
						keep_genes.append(g)
				self.genes = list(keep_genes)
				self.size = len(self.genes)
			
			## Check that the mutated group is connected
			if not self.is_connected():
				self.largest_connected()
		else:
			self.genes = [alter_gene]
			self.gene_pairs = []
			if selflevel == 'feature':
				self.features = [alter_feature]
				self.feature_pairs = []
			
		## Select the next gene from the group's neighbors
		neighbor_gene_scores = []
		if self.level == 'feature':
			## Get loops (self connections)
			for cur_gene, cur_feature in zip(self.genes, self.features):
				loop_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""(gene1 == cur_gene) & (gene2 == cur_gene)""")]
				loop_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""(gene2 == cur_gene) & (gene1 == cur_gene)""")]
				cur_loop_genes = loop_genes1
				cur_loop_genes.extend([pair for pair in loop_genes2 if (pair not in cur_loop_genes and pair not in loop_genes)])
				if len(cur_loop_genes) == 0:
					continue
				else:
					## Get gene scores
					for loop_gene, int_score in cur_loop_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == loop_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, loop_gene, cur_feature, weighted_score))
			
			## Get all neighbors of current group genes
			for cur_gene, cur_feature in zip(self.genes, self.features):
				neighbor_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""gene1 == cur_gene""")]
				neighbor_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""gene2 == cur_gene""")]
				cur_neighbor_genes = neighbor_genes1
				cur_neighbor_genes.extend([pair for pair in neighbor_genes2 if pair not in cur_neighbor_genes])
				neighbor_genes = [(g, s) for g, s in cur_neighbor_genes if g not in self.genes]
				if len(neighbor_genes) == 0:
					continue
				else:
					## Get gene scores
					for neighbor_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == neighbor_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, neighbor_gene, cur_feature, weighted_score))
			
			## Check that non-duplicate neighbor features are available
			all_neighbor_features = []
			for cg, ng, cf, ws in neighbor_gene_scores:
				neighbor_features = [row['feature'] for row in feature_table.where("""gene == ng""")]
				neighbor_features = [f for f in neighbor_features if f not in self.features]
				all_neighbor_features.extend(neighbor_features)
			
			if len(all_neighbor_features) == 0:
				if len(self.feature_pairs) == 0:
					# *** DEBUGGING / TESTING ***
					#print "create new group"
					
					new_group = self.new()
					new_group.create(data, gene_weight, int_weight)
					new_group.group_id = self.group_id
					self.update(new_group)
				else:
					fh.close()
					return self
		else:
			## Get neighbor genes
			for cur_gene in self.genes:
				neighbor_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""gene1 == cur_gene""")]
				neighbor_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""gene2 == cur_gene""")]
				cur_neighbor_genes = neighbor_genes1
				cur_neighbor_genes.extend([pair for pair in neighbor_genes2 if pair not in cur_neighbor_genes])
				neighbor_genes = [(g, s) for g, s in cur_neighbor_genes if g not in self.genes]
				if len(neighbor_genes) == 0:
					continue
				else:
					## Get gene scores
					for neighbor_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == neighbor_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, neighbor_gene, weighted_score))
		
		if len(neighbor_gene_scores) > 0:
			## Select gene to add based on weighted scores, checking that a duplicate is not added
			new_gene = None
			if self.level == 'feature':
				new_feature = None
				b3, ts3 = create_roulette_wheel([ws for cg, ng, cf, ws in neighbor_gene_scores])
			else:
				new_feature = 1
				b3, ts3 = create_roulette_wheel([ws for cg, ng, ws in neighbor_gene_scores])
			while (new_gene == None) or (new_feature == None):
				neighbor_index = spin_roulette_wheel(b3, ts3)
				group_gene = neighbor_gene_scores[neighbor_index][0]
				new_gene = neighbor_gene_scores[neighbor_index][1]
				
				if self.level == 'feature':
					group_feature = neighbor_gene_scores[neighbor_index][2]
					## Get all features mapped to selected gene, and randomly select the feature to add
					neighbor_feature_scores = [(row['feature'], row['score']) for row in feature_table.where("""gene == new_gene""")]
					neighbor_feature_scores = [(f, s) for f, s in neighbor_feature_scores if f not in self.features]
					if len(neighbor_feature_scores) > 0:
						b4, ts4 = create_roulette_wheel([s for f, s in neighbor_feature_scores])
						new_feature = neighbor_feature_scores[spin_roulette_wheel(b4, ts4)][0]
			
			## Add connected pair to the group list
			if self.level == 'feature':
				self.genes.append(new_gene)
				self.features.append(new_feature)
				self.gene_pairs.append([group_gene, new_gene])
				self.feature_pairs.append([group_feature, new_feature])
				self.size = len(self.features)
				self.fitness = 999.0
			else:
				self.genes.append(new_gene)
				self.gene_pairs.append([group_gene, new_gene])
				self.size = len(self.genes)
				self.fitness = 999.0
		else:
			if len(self.gene_pairs) == 0:
				# *** DEBUGGING / TESTING ***
				#print "create new group"
				
				new_group = self.new()
				new_group.create(data, gene_weight, int_weight)
				new_group.group_id = self.group_id
				self.update(new_group)
		
		fh.close()
		return self
	
	
	def mutate_swap(self):
		"""
		
		"""
		
		#print "Swapping..."
		if self.connected:
			return self
		else:
			if len(self.genes) > 1:
				## Randomly select feature to move to the front of the list
				rand = random.randint(0, (len(genes)-1))
				self.genes = [self.genes[rand]]+self.genes[:rand]+self.genes[rand+1:]
				self.size = len(self.genes)
				if self.level == 'feature':
					self.features = [self.features[rand]]+self.features[:rand]+self.features[rand+1:]
					self.size = len(self.features)
			
			return self
	
	
	def mutate_remove(self, data, gene_weight=0.5, int_weight=0.5):
		"""
		
		"""
		
		#print "Removing..."
		if self.connected:
			return self.mutate_remove_connected(data, gene_weight, int_weight)
		else:
			## Get the HDF5 file handle and data tables
			fh, gene_table, int_table, feature_table = data.get_tables()
			
			if len(self.genes) > 1:
				## Randomly select feature to remove from group
				gene_scores = []
				for cur_gene in self.genes:
					gene_scores.extend([row['score'] for row in gene_table.where("""gene == cur_gene""")])
				## Create a roulette wheel with gene scores and iterate through, the last to be selected will be removed from the group
				new_genes = [None]*len(self.genes)
				if self.level == 'feature':
					new_features = [None]*len(self.features)
				i = 0
				while i < (len(self.genes)-1):
					b1, ts1 = create_roulette_wheel(gene_scores)
					keep_index = spin_roulette_wheel(b1, ts1)
					if self.genes[keep_index] != None:
						new_genes[keep_index] = self.genes[keep_index]
						if self.level == 'feature':
							new_features[keep_index] = self.features[keep_index]
						self.genes[keep_index] = None
						i = i + 1
				self.genes = [g for g in new_genes if g != None]
				if self.level == 'feature':
					self.features = [f for f in new_features if f != None]
					self.size = len(self.features)
					self.fitness = 999.0
				else:
					self.size = len(self.genes)
					self.fitness = 999.0
			
			fh.close()
			return self
	
	
	def mutate_add(self, data, gene_weight=0.5, int_weight=0.5, step=1):
		"""
		
		"""
		
		#print "Adding..."
		if (gene_weight + int_weight != 1) or (gene_weight < 0) or (int_weight < 0):
			gene_weight = 0.5
			int_weight = 0.5
			
		if self.connected:
			return self.mutate_add_connected(data, gene_weight, int_weight)
		else:
			## Get the HDF5 file handle and data tables
			fh, gene_table, int_table, feature_table = data.get_tables()
			
			## Get all loops (self connections) for the current group genes
			loop_genes = []
			for cur_gene in self.genes:
				loop_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""(gene1 == cur_gene) & (gene2 == cur_gene)""")]
				loop_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""(gene2 == cur_gene) & (gene1 == cur_gene)""")]
				cur_loop_genes = [pair for pair in loop_genes1 if pair not in loop_genes]
				cur_loop_genes.extend([pair for pair in loop_genes2 if (pair not in cur_loop_genes and pair not in loop_genes)])
				loop_genes.extend(cur_loop_genes)
			
			## Get all neighbors of current group genes
			neighbor_genes = []
			for cur_gene in self.genes:
				neighbor_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""gene1 == cur_gene""")]
				neighbor_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""gene2 == cur_gene""")]
				cur_neighbor_genes = [pair for pair in neighbor_genes1 if pair not in neighbor_genes]
				cur_neighbor_genes.extend([pair for pair in neighbor_genes2 if (pair not in cur_neighbor_genes and pair not in neighbor_genes)])
				neighbor_genes.extend([(g, s) for g, s in cur_neighbor_genes if g not in self.genes])
			
			## Get all genes N steps away from current group genes
			if len(neighbor_genes) > 0:
				if step > 1:
					for i in range(2, step+1):
						new_neighbor_genes1 = []
						new_neighbor_genes2 = []
						for cur_gene, score in neighbor_genes:
							new_neighbor_genes1.extend([(row['gene2'], row['score']) for row in int_table.where("""gene1 == cur_gene""")])
							new_neighbor_genes2.extend([(row['gene1'], row['score']) for row in int_table.where("""gene2 == cur_gene""")])
						new_neighbor_genes1 = [(g, (s+score)/i) for g, s in new_neighbor_genes1 if (g not in self.genes and g not in [ng for ng, ns in neighbor_genes])]
						neighbor_genes.extend(new_neighbor_genes1)
						new_neighbor_genes2 = [(g, (s+score)/i) for g, s in new_neighbor_genes2 if (g not in self.genes and g not in [ng for ng, ns in neighbor_genes])]
						neighbor_genes.extend(new_neighbor_genes2)
			neighbor_genes.extend(loop_genes)
			
			if self.level == 'feature':
				## Check that non-duplicate neighbor features are available
				all_neighbor_features = []
				for ng, s in neighbor_genes:
					neighbor_features = [row['feature'] for row in feature_table.where("""gene == ng""")]
					neighbor_features = [f for f in neighbor_features if f not in self.features]
					all_neighbor_features.extend(neighbor_features)
					if len(all_neighbor_features) > 0:
						break
				
				if len(all_neighbor_features) > 0:
					## Get gene scores
					neighbor_gene_scores = []
					for cur_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == cur_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, weighted_score))
					
					## Select gene to add based on weighted scores, checking that a duplicate is not added
					new_gene = None
					new_feature = None
					while (new_gene == None) or (new_feature == None):
						b1, ts1 = create_roulette_wheel([ws for g, ws in neighbor_gene_scores])
						new_gene = neighbor_gene_scores[spin_roulette_wheel(b1, ts1)][0]
					
						## Get all features mapped to selected gene, and randomly select the feature to add
						neighbor_feature_scores = [(row['feature'], row['score']) for row in feature_table.where("""gene == new_gene""")]
						neighbor_feature_scores = [(f, s) for f, s in neighbor_feature_scores if f not in self.features]
						if len(neighbor_feature_scores) > 0:
							b2, ts2 = create_roulette_wheel([s for f, s in neighbor_feature_scores])
							new_feature = neighbor_feature_scores[spin_roulette_wheel(b2, ts2)][0]
					
					## Add the new gene and feature to the group
					self.genes.append(new_gene)
					self.features.append(new_feature)
					self.size = len(self.features)
					self.fitness = 999.0			
			else:
				## Check that non-duplicate neighbor genes are available
				cleaned_neighbor_genes = [(g, s) for g, s in neighbor_genes if g not in self.genes]
				if len(cleaned_neighbor_genes) > 0:
					## Get gene scores
					neighbor_gene_scores = []
					for cur_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == cur_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, weighted_score))
					
					new_gene = None
					while new_gene == None:
						b1, ts1 = create_roulette_wheel([ws for g, ws in neighbor_gene_scores])
						new_gene = neighbor_gene_scores[spin_roulette_wheel(b1, ts1)][0]
						if new_gene in self.genes:
							new_gene = None
					
					## Add the new gene to the group
					self.genes.append(new_gene)
					self.size = len(self.genes)
					self.fitness = 999.0
			
			fh.close()
			return self
	
	
	def mutate_alter(self, data, gene_weight=0.5, int_weight=0.5, step=1):
		"""
		
		"""
		
		#print "Altering..."
		if (gene_weight + int_weight != 1) or (gene_weight < 0) or (int_weight < 0):
			gene_weight = 0.5
			int_weight = 0.5
			
		if self.connected:
			return self.mutate_alter_connected(data, gene_weight, int_weight)
		else:
			## Get the HDF5 file handle and data tables
			fh, gene_table, int_table, feature_table = data.get_tables()
			
			## Randomly select gene to alter
			gene_scores = []
			for cur_gene in self.genes:
				gene_scores.extend([row['score'] for row in gene_table.where("""gene == cur_gene""")])
			## Create a roulette wheel with gene scores and iterate through, the last to be selected will be removed from the group
			new_genes = [None]*len(self.genes)
			if self.level == 'feature':
				new_features = [None]*len(self.features)
			i = 0
			keep_indices = []
			b1, ts1 = create_roulette_wheel(gene_scores)
			while i < (len(self.genes)-1):
				keep_index = spin_roulette_wheel(b1, ts1)
				if keep_index not in keep_indices:
					new_genes[keep_index] = self.genes[keep_index]
					if self.level == 'feature':
						new_features[keep_index] = self.features[keep_index]
					keep_indices.append(keep_index)
					i = i + 1
			alter_index = next(idx for idx, val in enumerate(self.genes) if idx not in keep_indices)
			alter_gene = self.genes[alter_index]
			
			## Get all loops (self connections) for the current group genes
			loop_genes = []
			loop_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""(gene1 == alter_gene) & (gene2 == alter_gene)""")]
			loop_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""(gene2 == alter_gene) & (gene1 == alter_gene)""")]
			loop_genes = loop_genes1
			loop_genes.extend([pair for pair in loop_genes2 if pair not in loop_genes])
			
			## Get all neighbors of the selected gene
			neighbor_genes = []
			neighbor_genes1 = [(row['gene2'], row['score']) for row in int_table.where("""gene1 == alter_gene""")]
			neighbor_genes2 = [(row['gene1'], row['score']) for row in int_table.where("""gene2 == alter_gene""")]
			neighbor_genes = neighbor_genes1
			neighbor_genes.extend([pair for pair in neighbor_genes2 if pair not in neighbor_genes])
			neighbor_genes = [(g, s) for g, s in neighbor_genes if g not in self.genes]
			
			## Get all genes N steps away from current group genes
			if len(neighbor_genes) > 0:
				if step > 1:
					for i in range(2, step+1):
						new_neighbor_genes1 = []
						new_neighbor_genes2 = []
						for cur_gene, score in neighbor_genes:
							new_neighbor_genes1.extend([(row['gene2'], row['score']) for row in int_table.where("""gene1 == cur_gene""")])
							new_neighbor_genes2.extend([(row['gene1'], row['score']) for row in int_table.where("""gene2 == cur_gene""")])
						new_neighbor_genes1 = [(g, (s+score)/i) for g, s in new_neighbor_genes1 if (g not in self.genes and g not in [ng for ng, ns in neighbor_genes])]
						neighbor_genes.extend(new_neighbor_genes1)
						new_neighbor_genes2 = [(g, (s+score)/i) for g, s in new_neighbor_genes2 if (g not in self.genes and g not in [ng for ng, ns in neighbor_genes])]
						neighbor_genes.extend(new_neighbor_genes2)
			neighbor_genes.extend(loop_genes)
			
			if self.level == 'feature':
				## Check that non-duplicate neighbor features are available
				all_neighbor_features = []
				for ng, s in neighbor_genes:
					neighbor_features = [row['feature'] for row in feature_table.where("""gene == ng""")]
					neighbor_features = [f for f in neighbor_features if f not in self.features]
					all_neighbor_features.extend(neighbor_features)
					if len(all_neighbor_features) > 0:
						break
				
				if len(all_neighbor_features) > 0:
					## Get gene scores
					neighbor_gene_scores = []
					for cur_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == cur_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, weighted_score))
					
					## Select gene to add based on weighted scores, checking that a duplicate is not added
					new_gene = None
					new_feature = None
					while (new_gene == None) or (new_feature == None):
						b1, ts1 = create_roulette_wheel([ws for g, ws in neighbor_gene_scores])
						new_gene = neighbor_gene_scores[spin_roulette_wheel(b1, ts1)][0]
						
						## Get all features mapped to selected gene, and randomly select the feature to add
						neighbor_feature_scores = [(row['feature'], row['score']) for row in feature_table.where("""gene == new_gene""")]
						neighbor_feature_scores = [(f, s) for f, s in neighbor_feature_scores if f not in self.features]
						if len(neighbor_feature_scores) > 0:
							b2, ts2 = create_roulette_wheel([s for f, s in neighbor_feature_scores])
							new_feature = neighbor_feature_scores[spin_roulette_wheel(b2, ts2)][0]
					
					## Add the new gene and feature to the group
					self.genes[alter_index] = new_gene
					self.features[alter_index] = new_feature
					self.size = len(self.features)
					self.fitness = 999.0
				else:
					## Remove the altered gene and feature
					self.genes = [g for g in new_genes if g != None]
					self.features = [f for f in new_features if f != None]
					self.size = len(self.features)
					self.fitness = 999.0
			else:
				## Check that non-duplicate neighbor genes are available
				cleaned_neighbor_genes = [(g, s) for g, s in neighbor_genes if g not in self.genes]
				if len(cleaned_neighbor_genes) > 0:
					## Get gene scores
					neighbor_gene_scores = []
					for cur_gene, int_score in neighbor_genes:
						gene_score = next(row['score'] for row in gene_table.where("""gene == cur_gene"""))
						weighted_score = (int_weight*(int_score/data.max_int_score)) + (gene_weight*(gene_score/data.max_gene_score))
						neighbor_gene_scores.append((cur_gene, weighted_score))
					
					new_gene = None
					while (new_gene == None):
						b1, ts1 = create_roulette_wheel([ws for g, ws in neighbor_gene_scores])
						new_gene = neighbor_gene_scores[spin_roulette_wheel(b1, ts1)][0]
						if new_gene in self.genes:
							new_gene = None
					
					## Add the new gene to the group
					self.genes[alter_index] = new_gene
					self.size = len(self.genes)
					self.fitness = 999.0
				else:
					## Remove the altered gene
					self.genes = [g for g in new_genes if g != None]
					self.size = len(self.genes)
					self.fitness = 999.0
			
			fh.close()
			return self


class PolyGAPop():
	"""
	
	"""
	
	def __init__(self, opts, params, data):
		"""
		
		"""
		
		## Define object attributes
		self.level = opts['analysis']
		self.pop_size = int(params['pop_size'])
		self.generation = 1
		self.max_generations = int(params['generations'])
		self.min_group_size = int(params['min_group_size'])
		self.max_group_size = int(params['max_group_size'])
		self.connected = params['connected']
		self.select_type = params['select_type']
		if self.select_type == 'hybrid':
			self.hybrid_top = int(params['hybrid_top'])
		else:
			self.hybrid_top = None
		self.cross_type = params['cross_type']
		self.cross_prob = float(params['cross_prob'])
		self.mut_prob = float(params['mut_prob'])
		self.elite_num = int(params['elite_num'])
		if 'step_size' in params:
			self.step_size = int(params['step_size'])
		else:
			self.step_size = 1
		self.migrants = int(params['migrants'])
		if 'R_script' in params:
			self.R_script = params['R_script']
		else:
			self.R_script = None
		self.gene_weight = float(params['node_weight'])
		self.int_weight = float(params['edge_weight'])
		self.generation_restart = params['generation_restart']
		if self.generation_restart != None:
			self.generation_restart = int(self.generation_restart)
		self.fitness_restart = params['fitness_restart']
		if self.fitness_restart != None:
			self.fitness_restart = float(self.fitness_restart)
		if 'nprocs' in params and params['nprocs'] > 1:
			self.nprocs = int(params['nprocs'])
		elif 'nprocs' in params and params['nprocs'] <= 1:
			self.nprocs = 1
		else:
			self.nprocs = None
		if 'parallel_nodes' in params:
			self.parallel_nodes = params['parallel_nodes']
		else:
			self.parallel_nodes = None
		self.global_best_fitness = 999.0
		self.local_best_fitness = 999.0
		self.global_best_generation = 1
		self.last_improvement = 1
		
		## Initialize R
		init_R(self.R_script)
		
		## Create Parallel Python (pp) job server
		self.job_server = init_PP(self.nprocs, self.parallel_nodes)
		
		## Create PyTables table to store the GA population
		self.file_name = opts['out']+'_out.h5'
		h5_ga = tables.openFile(self.file_name, mode='w', title='GA Output File - polyGA v1.0')
		group_ga = h5_ga.createGroup("/", 'ga', 'GA Population Data')
		
		## Create PyTables table definition for the GA population
		if self.connected:
			max_size = 2*self.max_group_size - 1
			if self.level == 'feature':
				class ga_pop(tables.IsDescription):
					generation = tables.Int32Col()
					group = tables.Int32Col()
					genes = tables.StringCol(itemsize=32, shape=(max_size,2))
					features = tables.StringCol(itemsize=32, shape=(max_size,2))
					gscores = tables.Float32Col(shape=(max_size,2))
					fscores = tables.Float32Col(shape=(max_size,2))
					escores = tables.Float32Col(shape=(1,max_size))
					size = tables.Int32Col()
					fitness = tables.Float32Col()
			else:
				class ga_pop(tables.IsDescription):
					generation = tables.Int32Col()
					group = tables.Int32Col()
					genes = tables.StringCol(itemsize=32, shape=(max_size,2))
					size = tables.Int32Col()
					fitness = tables.Float32Col()
		else:
			max_size = self.max_group_size
			if self.level == 'feature':
				class ga_pop(tables.IsDescription):
					generation = tables.Int32Col()
					group = tables.Int32Col()
					genes = tables.StringCol(itemsize=32, shape=(1,max_size))
					features = tables.StringCol(itemsize=32, shape=(1,max_size))
					size = tables.Int32Col()
					fitness = tables.Float32Col()
			else:
				class ga_pop(tables.IsDescription):
					generation = tables.Int32Col()
					group = tables.Int32Col()
					genes = tables.StringCol(itemsize=32, shape=(1,max_size))
					size = tables.Int32Col()
					fitness = tables.Float32Col()

		## Create the GA population table and close the HDF5 file
		ga_table = h5_ga.createTable(group_ga, 'groups', ga_pop, "GA Population")
		if self.connected:
			ga_table.attrs.connected = 'yes'
			if self.level == 'feature':
				ga_table.attrs.level = 'feature'
			else:
				ga_table.attrs.level = 'gene'
		else:
			ga_table.attrs.connected = 'no'
			if self.level == 'feature':
				ga_table.attrs.level = 'feature'
			else:
				ga_table.attrs.level = 'gene'
		h5_ga.close()
		
		
		## Create the gene network
		G_full = data.get_graph()
		G_nodes_only = data.get_graph(edges=False)
		
		## Create the GA population
		self.pop = []
		for i in range(self.pop_size):
			self.pop.append(PolyGAInd(self.connected, self.level, self.min_group_size, self.max_group_size))
		
		if self.job_server != None:
			pp_jobs = []
			for ind in self.pop:
				b1, ts1 = data.gene_scores_wheel
				gene1 = data.genes[spin_roulette_wheel(b1, ts1)]
				ind.graph.add_node(gene1)
				if self.connected:
					G = G_full
				else:
					G = G_nodes_only
				pp_jobs.append(self.job_server.submit(PolyGAInd.create, (ind, G, data, self.gene_weight, self.int_weight), (create_roulette_wheel, spin_roulette_wheel), ("random",)))
			
			self.job_server.wait()
			for i, job in enumerate(pp_jobs):
				self.pop[i] = job()
				self.pop[i].group_id = i + 1
				job = None
			pp_jobs = None
		else:
			for i, ind in enumerate(self.pop):
				ind.create(data, self.gene_weight, self.int_weight)
				ind.group_id = i + 1
		
		return None
	
	
	def get_from_file(self, generation=None):
		"""
		
		"""
		
		return None
	
	
	def write_to_file(self):
		"""
		
		"""
		
		## Get PyTables
		fh = tables.openFile(self.file_name, mode='a')
		ga_table = fh.root.ga.groups
		
		## Initialize a feature group
		feature_group = ga_table.row
		
		for ind in self.pop:
			feature_group['generation'] = self.generation
			feature_group['group'] = ind.group_id
			if ind.connected:
				genes = [['NA','NA']]*ind.max_size
				genes[:len(ind.gene_pairs)] = ind.gene_pairs
				feature_group['genes'] = genes
			else:
				genes = ['NA']*ind.max_size
				genes[:len(ind.genes)] = ind.genes
				feature_group['genes'] = genes
			if self.level == 'feature':
				if ind.connected:
					features = [['NA','NA']]*ind.max_size
					features[:len(ind.feature_pairs)] = ind.feature_pairs
					feature_group['features'] = features
				else:
					features = ['NA']*ind.max_size
					features[:len(ind.features)] = ind.features
					feature_group['features'] = features
			feature_group['size'] = ind.size
			feature_group['fitness'] = ind.fitness
			feature_group.append()

		ga_table.flush()
		fh.close()
		return None
	
	
	def fitness(self):
		""" 
		This function uses the Rpy2 module to calculate the fitness (e.g., association p-value) for each 
		feature group in the GA population. The feature groups are passed to R and a list of lists is 
		created. The R script specified in the polyGA configuration file is used to calculate the fitness. 
		For more details about how this R script should be created please see the polyGA documentation.
		
		"""

		## Create list of feature groups in R
		#rob.r("r_ga_pop = list()")
		groups_to_update = []
		r_groups = []
		for ind in [i for i in self.pop if i.fitness == 999]:
			groups_to_update.append(ind.group_id)
			if ind.level == 'feature':
				if ind.connected:
					grp = [f for pair in ind.feature_pairs for f in pair]
					r_grp = rob.StrVector(list(set(grp)))
					r_groups.append(r_grp)
				else:
					r_grp = rob.StrVector(ind.features)
					r_groups.append(r_grp)
			else:
				if ind.connected:
					grp = [g for pair in ind.gene_pairs for g in pair]
					r_grp = rob.StrVector(list(set(grp)))
					r_groups.append(r_grp)
				else:
					r_grp = rob.StrVector(ind.genes)
					r_groups.append(r_grp)
		
		r_code =  "r_ga_pop = c(list(" + "), list(".join([g.r_repr() for g in r_groups]) + "))"
			
		# *** DEBUGGING / TESTING ***
		#print r_code
			
		rob.r(r_code)
		
		# *** DEBUGGING / TESTING ***
		#rob.r("r_fit_vals = c()")
		#rob.r("r_fit_vals = runif("+str(len(cur_gen))+")")
		#test = rob.r("getDoParWorkers()")		
		#print test
		#test2 = rob.r("clusterEvalQ(cl, ls())")
		#print test2
		#test3 = rob.r("clusterEvalQ(cl, gc(TRUE))")
		#print test3
		
		
		## Perform the association tests in R for each feature group
		fit_calc_failed = 0
		try:
			r_code = "r_fit_vals = r_get_fitness(r_ga_pop)"
			rob.r(r_code)
		except rob.rinterface.RRuntimeError as rerr:
			if "cannot allocate vector of size" in str(rerr):
				print "R returned a memory error!"
			else:
				print "R returned a runtime error!"
				print str(rerr)
				fit_calc_failed = 1
		except:
			print "Unknown error returned by R!"
			fit_calc_failed = 1

		if fit_calc_failed == 1:
			print "P-value calculation failed!"
			sys.exit(5)
		
		fit_vals = list(rob.r('r_fit_vals'))
		fit_vals = [f+1e-40 for f in fit_vals]
		
		
		# *** DEBUGGING / TESTING ***
		#print fit_vals
		#fit_vals = [random.random() for i in range(len(groups_to_update))]
		
		groups_fit_vals = zip(groups_to_update, fit_vals)
		groups_fit_vals_dict = dict(groups_fit_vals)
	
		## Insert p-values into GA population table
		for ind in [i for i in self.pop if i.fitness == 999.0]:
			ind.fitness = groups_fit_vals_dict[ind.group_id]
		
		return None
	
	
	def select(self, data):
		""" 
		Selects the parents, based on fitness values, for the next generation of the GA. There are three
		different selection methods: 'truncate' (simply take the most fit half of the current generation), 
		'roulette' (probabilistically select the feature groups based on fitness values--default), and 'hybrid' 
		(take a specified number--using the 'num_top' parameter--of the most fit groups, and probabilistically 
		select the rest).
		
		"""
		
		next_gen = self.generation + 1
		
		## Get current GA population
		#cur_pop = self.pop[:]
		
		## Create new GA population for the next generation
		## Create migrants
		migrants = []
		for i in range(self.migrants):
			migrants.append(PolyGAInd(self.connected, self.level, self.min_group_size, self.max_group_size))
		
		new_pop = []
		if self.job_server != None:
			pp_jobs = []
			for ind in migrants:
				pp_jobs.append(self.job_server.submit(PolyGAInd.create, (ind, data, self.gene_weight, self.int_weight), (create_roulette_wheel, spin_roulette_wheel), ("random",)))
			
			self.job_server.wait()
			for i, job in enumerate(pp_jobs):
				new_ind = job()
				new_ind.group_id = i + 1
				new_pop.append(new_ind)
				job = None
			pp_jobs = None
		else:
			for i, ind in enumerate(migrants):
				ind.create(data, self.gene_weight, self.int_weight)
				ind.group_id = i + 1
				new_pop.append(ind)
		
		## Select parents for next generation
		if self.select_type == 'truncate':
			self.pop.sort(key=lambda x: float(x.fitness))
			num_parents = self.pop_size/2 - self.migrants
			for idx, ind in enumerate(self.pop[:num_parents], self.migrants+1):
				ind.group_id = idx
				new_pop.append(ind)
			self.pop = new_pop
			self.generation = next_gen
		elif self.select_type == 'hybrid':
			self.pop.sort(key=lambda x: float(x.fitness))
			num_parents = self.pop_size/2 - self.migrants
			if self.hybrid_top > num_parents:
				num_top = int(num_parents/2)
			else:
				num_top = self.hybrid_top
			num_parents = num_parents - num_top
			
			## Append the hybrid_top individuals to the new population
			for i in xrange(num_top):
				self.pop[i].group_id = i + 1 + self.migrants
				new_pop.append(self.pop[i])
			
			## Get the fitness values for the current population and create a roulette wheel
			fitvals = []
			for ind in self.pop:
				if -math.log10(ind.fitness) < 0.05:
					fitvals.append(0.05)
				else:
					fitvals.append(-math.log10(ind.fitness))
			
			# *** DEBUGGING / TESTING ***
			#print "parents: "+str(num_parents)
			#print "migrants: "+str(self.migrants)
			#print "top: "+str(num_top)
			
			## Spin the roulette wheel to select the parents
			for i in xrange(num_parents):
				b1, ts1 = create_roulette_wheel(fitvals)
				group_index = spin_roulette_wheel(b1, ts1)
				new_ind = self.pop[group_index].copy()
				new_ind.group_id = i + 1 + num_top + self.migrants
				new_pop.append(new_ind)
			self.pop = new_pop
			self.generation = next_gen
		else:
			num_parents = self.pop_size/2 - self.migrants
			
			## Get the fitness values for the current population and create a roulette wheel
			fitvals = []
			for ind in self.pop:
				if -math.log10(ind.fitness) < 0.05:
					fitvals.append(0.05)
				else:
					fitvals.append(-math.log10(ind.fitness))
			
			for i in xrange(num_parents):
				b1, ts1 = create_roulette_wheel(fitvals)
				group_index = spin_roulette_wheel(b1, ts1)
				new_ind = self.pop[group_index].copy()
				new_ind.group_id = i + 1 + self.migrants
				new_pop.append(new_ind)
			self.pop = new_pop
			self.generation = next_gen
		
		return None
	
	
	def crossover(self):
		"""
		
		"""
		
		num_parents = self.pop_size/2
		
		if self.connected:
			offspring = []
			for idx, ind in enumerate(self.pop, 1):
				new_ind = ind.copy()
				new_ind.group_id = num_parents + idx
				offspring.append(new_ind)
			self.pop.extend(offspring)
		else:
			offspring = [None]*num_parents
			parents = []
			group_ids = range(num_parents)
			while len(group_ids) > 0:
				pair = tuple(random.sample(group_ids, 2))
				parents.append((self.pop[pair[0]], self.pop[pair[1]]))
				group_ids.remove(pair[0])
				group_ids.remove(pair[1])
			
			if self.job_server != None:
				pp_jobs = []
				for ind, parent2 in parents:
					pp_jobs.append(self.job_server.submit(PolyGAInd.crossover, (ind, parent2, self.cross_type, self.cross_prob), (), ("random",)))
			
				self.job_server.wait()
				i = 0
				for job in pp_jobs:
					children = job()
					offspring[i] = children[0]
					offspring[i].group_id = i + 1 + num_parents
					i = i + 1
					offspring[i] = children[1]
					offspring[i].group_id = i + 1 + num_parents
					i = i + 1
					job = None
				pp_jobs = None
			else:
				i = 0
				for ind, parent2 in parents:
					children = ind.crossover(parent2, self.cross_type, self.cross_prob)
					offspring[i] = children[0]
					offspring[i].group_id = i + 1 + num_parents
					i = i + 1
					offspring[i] = children[1]
					offspring[i].group_id = i + 1 + num_parents
					i = i + 1
					
					# *** DEBUGGING / TESTING ***
					#print ind
					#print parent2
					#print "crossed"
					#print children[0]
					#print children[1]
			
			self.pop.extend(offspring)
		
		return None
	
	
	def mutate(self, data):
		"""
		
		"""
		
		num_parents = self.pop_size/2
		self.pop.sort(key=lambda x: float(x.fitness))
		
		## Create mutated GA population
		new_pop = []
		
		## Insert elite individuals into the population
		for idx, ind in enumerate(self.pop[:self.elite_num], 1):
			new_ind = ind.copy()
			new_ind.group_id = idx
			new_pop.append(new_ind)
		
		for idx, ind in enumerate(self.pop[self.elite_num:], self.elite_num+1):
			if random.random() > self.mut_prob:
				new_ind = ind.copy()
				new_ind.group_id = idx
				new_pop.append(new_ind)
			else:
				new_ind = ind.copy()
				new_ind.group_id = idx
				
				## Determine which operations are allowed (Remove, Add, Alter, Swap)
				opts = []
				grp_size = new_ind.size
				if self.cross_type == 'uniform':
					if grp_size < self.max_group_size and grp_size > self.min_group_size:
						opts = ['add', 'remove', 'alter']
					elif grp_size == self.min_group_size and grp_size < self.max_group_size:
						opts = ['add', 'alter']
					elif grp_size == self.max_group_size and grp_size > self.min_group_size:
						opts = ['alter', 'remove']
					elif grp_size == self.max_group_size and grp_size == self.min_group_size:
						opts = ['alter']
					elif grp_size < self.min_group_size:
						opts = ['add']
				elif self.cross_type == 'single-point':
					if grp_size < self.max_group_size and grp_size > self.min_group_size:
						opts = ['add', 'remove', 'alter', 'swap']
					elif grp_size == self.min_group_size and grp_size < self.max_group_size:
						opts = ['add', 'alter', 'swap']
					elif grp_size == self.max_group_size and grp_size > self.min_group_size:
						opts = ['alter', 'remove', 'swap']
					elif grp_size == self.max_group_size and grp_size == self.min_group_size:
						opts = ['alter', 'swap']
					elif grp_size < self.min_group_size:
						opts = ['add']
				
				## Randomly choose the operation to perform
				rand = random.randint(0, (len(opts)-1))
				opt = opts[rand]
				
				# *** DEBUGGING / TESTING ***
				#print opt
				
				if self.job_server != None:
					pp_jobs = []
					if opt == 'swap':
						pp_jobs.append(self.job_server.submit(PolyGAInd.mutate_swap, (new_ind,), (), ("random",)))
					elif opt == 'remove':
						pp_jobs.append(self.job_server.submit(PolyGAInd.mutate_remove, (new_ind, data, self.gene_weight, self.int_weight), (create_roulette_wheel, spin_roulette_wheel), ("random",)))
					elif opt == 'add':
						pp_jobs.append(self.job_server.submit(PolyGAInd.mutate_add, (new_ind, data, self.gene_weight, self.int_weight, self.step_size), (create_roulette_wheel, spin_roulette_wheel), ("random",)))
					elif opt == 'alter':
						pp_jobs.append(self.job_server.submit(PolyGAInd.mutate_alter, (new_ind, data, self.gene_weight, self.int_weight, self.step_size), (create_roulette_wheel, spin_roulette_wheel), ("random",)))
				else:
					if opt == 'swap':
						new_ind.mutate_swap()
						new_pop.append(new_ind)
					elif opt == 'remove':
						new_ind.mutate_remove(data, self.gene_weight, self.int_weight)
						new_pop.append(new_ind)
					elif opt == 'add':
						new_ind.mutate_add(data, self.gene_weight, self.int_weight)
						new_pop.append(new_ind)
					elif opt == 'alter':
						new_ind.mutate_alter(data, self.gene_weight, self.int_weight)
						new_pop.append(new_ind)
					new_ind = None
		
		if self.job_server != None:
			self.job_server.wait()
			for job in pp_jobs:
				new_ind = job()
				
				# *** DEBUGGING / TESTING ***
				if new_ind == None:
					print "MUTANT CREATION FAILED!!!"
					sys.exit()
				
				new_pop.append(new_ind)
				job = None
			pp_jobs = None
		
		new_pop.sort(key=lambda x: int(x.group_id))
		self.pop = new_pop
		
		return None
	
	
	def restart(self):
		"""
		
		"""
		
		gen_restart_tag = False
		fit_restart_tag = False
		
		# *** DEBUGGING / TESTING ***
		#print "gen", self.generation
		#print self.local_best_fitness
		#print "last improvement", self.last_improvement
		#print self.global_best_fitness
		#print self.global_best_generation
		
		## Update fitness improvement variables
		cur_min_fitness = min([ind.fitness for ind in self.pop])
		if cur_min_fitness < self.local_best_fitness:
			self.local_best_fitness = cur_min_fitness
			self.last_improvement = self.generation + 1
			if self.local_best_fitness < self.global_best_fitness:
				self.global_best_fitness = self.local_best_fitness
				self.global_best_generation = self.generation + 1
		
		## Check if search should re-start
		if self.generation_restart != None:
			if (self.generation - self.last_improvement) > self.generation_restart:
				self.last_improvement = self.generation + 1
				# Could also keep best_fitness as the global best, rather than resetting its value
				self.local_best_fitness = 999.0
				gen_restart_tag = True
		if self.fitness_restart != None:
			if self.local_best_fitness < self.fitness_restart:
				# *** DEBUGGING / TESTING ***
				#print "fit restart", self.local_best_fitness
				
				self.last_improvement = self.generation + 1
				# Could also keep best_fitness as the global best, rather than resetting its value
				self.local_best_fitness = 999.0
				fit_restart_tag = True
				
		return gen_restart_tag, fit_restart_tag
	
	
	def restart_PP(self):
		"""
		
		"""
		if self.job_server != None:
			self.job_server = shutdown_PP(self.job_server)
			self.job_server = init_PP(self.nprocs, self.parallel_nodes, True)
		
		return None


def create_roulette_wheel(scores):
	""" 
	The function takes a list of numeric values (scores) and returns a list containing interval boundaries 
	(i.e., the sizes of the pie slices on a roulette wheel), and the total sum of all scores. The intervals 
	are used to probabilistically select items based on their 'scores'.
	"""

	boundaries = []
	score_sum = 0
	for x in scores:
		score_sum = score_sum + x
		boundaries.append(score_sum)
	
	return boundaries, score_sum


def spin_roulette_wheel(boundaries, total_score):
	rand = random.uniform(0, total_score)
	index = next(idx for idx,val in enumerate(boundaries) if val >= rand)
	
	return index

def is_neighbor(neighbors, neighbor_dict, dest_node, visited_nodes):
	"""
	
	"""
	new_neighbors = neighbors.difference(visited_nodes)
	if len(new_neighbors) == 0:
		return False
	elif dest_node in new_neighbors:
		return True
	else:
		next_neighbors = []
		for n in new_neighbors:
			visited_nodes.add(n)
			next_neighbors.extend(neighbor_dict[n])
		return is_neighbor(set(next_neighbors), neighbor_dict, dest_node, visited_nodes)


def is_path(neighbor_dict, node1, node2):
	"""
	
	"""
	
	cur_neighbors = set(neighbor_dict[node1])
	visited_nodes = set([node1])
	
	return is_neighbor(cur_neighbors, neighbor_dict, node2, visited_nodes)


class RLoadError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		print ": \n"+self.value
		return None


def init_R(R_script):
	"""
	
	"""
	
	## Start R and run the initialization script
	print "Initializing R ... "
	print "--------------------------------------------------------"
	print "Sourcing R script: " + R_script
	r_code = "source('" + R_script + "')"
	try:
		rob.r(r_code)
	except rob.rinterface.RRuntimeError as rerr:
		print "   "+str(rerr)
		print "\nThere was an error loading the R script (" + R_script + "): \n"
		raise RLoadError("Please see the polyGA documentation for details about the required format of the R script.\n")
	else:
		r_libs_loaded = rob.r("success")
		if not r_libs_loaded[0]:
			raise RLoadError("The required R libraries failed to load!")
	print "--------------------------------------------------------\n"
	sys.stdout.flush()
	
	return None

def shutdown_R(R_script):
	"""
	
	"""
	if R_script != None:
		rob.r("shutdownR()")
	return None


def init_PP(nprocs, parallel_nodes, quiet=False):
	"""
	
	"""
	
	## Start a Parallel Python (pp) server for parallel processing
	if nprocs != None and nprocs <= 1:
		job_server = None
	else:
		import pp
		
		if parallel_nodes != None and len(parallel_nodes) > 0:
			job_server = pp.Server(ppservers=parallel_nodes)
			wp = str(sum(job_server.get_active_nodes().values()))
			cn = ', '.join(job_server.get_active_nodes().keys())
		elif nprocs != None and nprocs > 1:
			job_server = pp.Server(nprocs)
			wp = str(sum(job_server.get_active_nodes().values()))
			cn = "localhost"
		else:
			job_server = pp.Server()
			wp = str(sum(job_server.get_active_nodes().values()))
			cn = "localhost"
		
		if quiet == False:
			## Print message
			print "Initializing Parallel Python (pp) job server ... "
			print "--------------------------------------------------------"
			print "%19s:  %s" % ("Worker processes", wp)
			print "%19s:  %s" % ("Computing nodes", cn)
			print "--------------------------------------------------------\n"
			sys.stdout.flush()
	
	return job_server


def shutdown_PP(job_server):
	""
	
	""
	
	job_server.destroy()
	return None


def ga_main(opts, params):
	"""
	
	"""
	
	## Create the PolyGA data oject
	#t0 = time.time()
	if opts['data'] != None:
		data = polyga_utils.PolyGAData.from_hdf5(opts['data'])
	else:
		data = polyga_utils.PolyGAData.from_delim(opts['out'], opts['genes'], opts['interactions'], opts['features'])
	
	## Create the GA population
	#t1 = time.time()
	population = PolyGAPop(opts, params, data)
	#t1b = time.time()
	
	print "Searching for Polygenic Associations ...\n"
	sys.stdout.flush()
	
	## Initialize variables association with algorithm progress and search re-start
	gen_restart_count = 0
	fit_restart_count = 0
	if population.max_generations < 10:
		progress_interval = 1
	else:
		progress_interval = population.max_generations/10
		
	## Initialize the migrant parameter used for search re-starts
	migrant_param = population.migrants
	
	## Run the GA for the specified number of generations
	while population.generation < population.max_generations:			
		# *** DEBUGGING / TESTING ***
		#print "generation: "+str(population.generation)
		
		#t2 = time.time()
		population.fitness()
		
		#t3 = time.time()
		population.write_to_file()
		
		## Print progress message
		if population.generation % progress_interval == 0:
			fstr = '%' + str(len(str(population.max_generations))+3) + 's' + ' generations completed ...'
			print fstr % str(population.generation)
			sys.stdout.flush()
		
		## Check whether the search should be re-initializd by increasing the number of migrants
		#t4 = time.time()
		gen_restart_tag, fit_restart_tag = population.restart()
		restart_tag = any((gen_restart_tag, fit_restart_tag))
		if gen_restart_tag:
			gen_restart_count = gen_restart_count + 1
		if fit_restart_tag:
			fit_restart_count = fit_restart_count + 1
		
		if restart_tag == False:
			#t5 = time.time()
			population.select(data)
		else:
			#t5 = time.time()
			population.migrants = population.pop_size/2
			population.select(data)
			population.migrants = migrant_param
		
		#t6 = time.time()
		population.crossover()
		
		#t7 = time.time()
		population.mutate(data)
		
		#t8 = time.time()
		#if population.job_server != None:
		#	population.restart_PP()
		
		#t9 = time.time()
		#if population.generation % 100 == 0:
		#	gc.collect()
		
		# *** DEBUGGING / TESTING ***
		#for ind in population.pop:
		#	print ind
		#t10 = time.time()
		#print "Data load time: "+str(t1 - t0)
		#print "Create population time: "+str(t1b - t1)
		#print "Fitness time: "+str(t3 - t2)
		#print "Data write time: "+str(t4 - t3)
		#print "Restart calc time: "+str(t5 - t4)
		#print "Selection time: "+str(t6 - t5)
		#print "Crossover time: "+str(t7 - t6)
		#print "Mutate time: "+str(t8 - t7)
		#print "PP restart time: "+str(t9 - t8)
		#print "Garbage collection time: "+str(t10 - t9)
		#print "Generation time: "+str(t10 - t2)
		#sys.stdout.flush()
		
	## Calculate fitness for the last population and write to file
	population.fitness()
	population.write_to_file()
	
	## Print progress message
	fstr = '%' + str(len(str(population.max_generations))+3) + 's' + ' generations completed ... Done.'
	print fstr % str(population.max_generations)
	sys.stdout.flush()
	
	## Print the smallest p-value found in the search
	print "\nBest fitness found: "+str(population.global_best_fitness)+"\n"
	
	## Print info about search restarts
	if gen_restart_count + fit_restart_count > 0:
		print "Total number of search re-starts: " + str(gen_restart_count + fit_restart_count)
		print "%30s: %s" % ("Due to generation threshold", str(gen_restart_count))
		print "%30s: %s\n" % ("Due to fitness threshold", str(fit_restart_count))
	
	## Stop any R processes
	shutdown_R(population.R_script)
	
	## Stop Parallel Python (pp) job server
	if population.job_server != None:
		shutdown_PP(population.job_server)
	
	## Print message about output file
	print "Output file ("+population.file_name+") created.\n"
				
	return None


if __name__ == '__main__':
	## Start testing procudure
	pass

