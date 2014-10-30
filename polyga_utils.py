#!/usr/bin/env python

"""

PolyGA - version 1.2b - Last modified: 31-DEC-2013

This file contains the input and output functions for the PolyGA 
program, as well as functions related to error checking.

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
import math
import time
import tables
import networkx as nx
import polyga_core 


class DataFormatError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		print self.value
		return None


class PolyGAData():
	"""
	
	"""
	
	def __init__(self, h5_file, level, gene_rows, int_rows, feature_rows, genes, gene_scores_wheel, max_gene_score, max_int_score):
		"""
		
		"""
		
		self.file_name = h5_file
		self.level = level
		self.num_genes = gene_rows
		self.num_interactions = int_rows
		self.num_feature_pairs = feature_rows
		self.genes = genes
		self.gene_scores_wheel = gene_scores_wheel
		self.max_gene_score = max_gene_score
		self.max_int_score = max_int_score
		self.graph = None
		
		return None
	
	
	@classmethod
	def from_hdf5(cls, h5_file):
		"""
		
		"""
		
		try:
			fh = tables.openFile(h5_file, mode='r')
		except IOError as err:
			print "\nError opening the HDF5 data file containing the gene network!\n"
			raise IOError(str(err))
		else:
			correct_format, level, gene_rows, int_rows, feature_rows = table_format_check(fh, 'input')
			if not correct_format:
				raise DataFormatError("Data did not load successfully!\nTry reloading the data from tab-delimited files.\n")
			else:
				## Print data loading message
				print "Loading data from HDF5 file ..."
				
				node_table = fh.root.gene_network.genes
				edge_table = fh.root.gene_network.interactions
				nodes = [(row['gene'], row['score']) for row in node_table.iterrows()]
				edges = [(row['gene1'], row['gene2'], row['score']) for row in edge_table.iterrows()]
				
				## Create NetworkX graph
				t0 = time.time()
				G = nx.Graph()
				for g, s in nodes:
					G.add_node(g, weight=s)
				
				for g1, g2, s in edges:
					G.add_edge(g1, g2, weight=s)
				
				## Get network attributes
				node_attrs = nx.get_node_attributes(G, 'weight')
				genes = [g for g, w in node_attrs.items()]
				gene_scores = [w for g, w in node_attrs.items()]
				gene_scores_wheel = polyga_core.create_roulette_wheel(gene_scores)
				max_gene_score = node_table.attrs.max_gene_score
				edge_table = fh.root.gene_network.interactions
				max_int_score = edge_table.attrs.max_int_score
				fh.close()
				
				## Print data loading message
				if level == 'feature':
					max_num = max(len(str(gene_rows)), len(str(int_rows)), len(str(feature_rows)))
				else:
					max_num = max(len(str(gene_rows)), len(str(int_rows)))
				print "--------------------------------------------------------"
				fstr = '%' + str(max_num+3) + 's' + ' %s'
				print fstr % (str(gene_rows), "gene annotations loaded.")
				print fstr % (str(int_rows), "gene-gene interactions loaded.")
				if level == 'feature':
					print fstr % (str(feature_rows), "gene-feature pairs loaded.")
				print "--------------------------------------------------------\n"
				
				return cls(h5_file, level, gene_rows, int_rows, feature_rows, genes, gene_scores_wheel, max_gene_score, max_int_score)
	
	
	@classmethod
	def from_delim(cls, out, genes, ints, features=None):
		"""
		
		"""
		
		## Check existence of out file
		out_file = out + '.h5'
		if os.path.exists(out_file):
			raise IOError("The output file ("+out_file+") already exists!\n")
		else:
			h5_file = out_file
			
		## Create PyTables table definitions
		class network_nodes(tables.IsDescription):
			gene = tables.StringCol(itemsize=16)
			score = tables.Float32Col()
		
		class node_attributes(tables.IsDescription):
			gene = tables.StringCol(itemsize=16)
			feature = tables.StringCol(16)
			fscore = tables.Float32Col()
		
		class network_edges(tables.IsDescription):
			gene1 = tables.StringCol(itemsize=16) 
			gene2 = tables.StringCol(itemsize=16)
			score = tables.Float32Col()
		
		## Check existence of input files and read data
		if features == None:
			if not os.path.exists(genes) or not os.path.exists(ints):
				raise IOError("At least one of the input files does not exist!\n")
		else:
			if not os.path.exists(genes) or not os.path.exists(ints) or not os.path.exists(features):
				raise IOError("At least one of the input files does not exist!\n")
		
		## Print data loading message
		print "Loading data from text files ..."
		
		## Create the HDF5 file and group
		h5 = tables.openFile(h5_file, mode='w', title='GA Network Data File - PolyGA v1.2b')
		group_network = h5.createGroup("/", 'gene_network', 'GA Network Data')
		
		## Create PyTables tables
		node_table = h5.createTable(group_network, 'genes', network_nodes, "Gene Annotations")
		edge_table = h5.createTable(group_network, 'interactions', network_edges, "Gene Interactions")
		
		## Create NetworkX Graph
		G = nx.Graph()
		
		### Load the gene data into the graph ###
		
		with open(genes, 'r') as fh:
			## Skip the header line (column names)
			fh.readline()
		
			gene_rows = 0
			for line in fh:
				line = line.rstrip()
				fields = line.split('\t')
				G.add_node(fields[0], weight=float(fields[1]), features={})
				gene_rows = gene_rows + 1
		
		### Load the gene-gene interaction data into the graph ###
		
		with open(ints, 'r') as fh:
			## Skip the header line (column names)
			fh.readline()
			
			int_rows = 0
			for line in fh:
				line = line.rstrip()
				fields = line.split('\t')
				edge_data = G.get_edge_data(fields[0], fields[1])
				if edge_data == None:
					G.add_edge(fields[0], fields[1], weight=float(fields[2]))
				elif fields[2] > edge_data['weight']:
					G.add_edge(fields[0], fields[1], weight=float(fields[2]))
				int_rows = int_rows + 1
		
		if features != None:
			level = 'feature'
			node_table.attrs.level = 'feature'
			feature_table = h5.createTable(group_network, 'features', node_attributes, "Feature Annotations")
			
			### Load the gene-feature map into the graph ###
			feature_rows = 0
			with open(features, 'r') as fh:
				## Skip the header line (column names)
				fh.readline()
				
				for line in fh:
					line = line.rstrip()
					fields = line.split('\t')
					G.node[fields[1]]['features'].update({fields[0]:float(fields[2])})
					feature_rows = feature_rows + 1
		else:
			node_table.attrs.level = 'gene'
			feature_rows = 0
			level = 'gene'
		
		### Network clean-up (check data consistency)
		if level == 'feature':
			for n in G.nodes():
				if not G.node[n].has_key('weight') or G.node[n]['features'] == {}:
					G.remove_node(n)
		else:
			for n in G.nodes():
				if not G.node[n].has_key('weight'):
					G.remove_node(n)
		
		### Write network to PyTables tables
		
		## Write nodes
		net_node = node_table.row
		for g in G.nodes():
			net_node['gene'] = g
			net_node['score'] = G.node[g]['weight']
			net_node.append()
		node_table.flush()
		
		## Create index
		node_table.cols.gene.createIndex()
		
		## Write edges
		net_edge = edge_table.row
		for g1, g2, in G.edges():
			net_edge['gene1'] = g1
			net_edge['gene2'] = g2
			net_edge['score'] = G.edge[g1][g2]['weight']
			net_edge.append()
		edge_table.flush()
		
		## Create indexes
		edge_table.cols.gene1.createIndex()
		edge_table.cols.gene2.createIndex()
		
		## Write features
		if level == 'feature':
			num_features = 0
			net_feature = feature_table.row
			for g in G.nodes():
				for f, fs in G.node[g]['features'].items():
					net_feature['gene'] = g
					net_feature['feature'] = f
					net_feature['fscore'] = fs
					net_feature.append()
					num_features =  num_features + 1
			feature_table.flush()
		
		## Create index
		feature_table.cols.gene.createIndex()
		
		## Get network attributes and set table attributes
		node_attrs = nx.get_node_attributes(G, 'weight')
		edge_attrs = nx.get_edge_attributes(G, 'weight')
		genes = [g for g, w in node_attrs.items()]
		gene_scores = [w for g, w in node_attrs.items()]
		edge_scores = [w for e, w in edge_attrs.items()]
		gene_scores_wheel = polyga_core.create_roulette_wheel(gene_scores)
		max_gene_score = max(gene_scores)
		max_int_score = max(edge_scores)
		
		node_table.attrs.max_gene_score = max_gene_score
		edge_table.attrs.max_int_score = max_int_score
		
		### Close the HDF5 file
		h5.close()
		
		## Print data loading message
		num_genes = G.number_of_nodes()
		num_ints = G.number_of_edges()
		if level == 'feature':
			max_num = max(len(str(num_genes)), len(str(num_ints)), len(str(num_features)))
		else:
			max_num = max(len(str(num_genes)), len(str(num_ints)))
		print "--------------------------------------------------------"
		fstr = '%' + str(max_num+3) + 's' + ' %s'
		print fstr % (str(num_genes), "gene annotations loaded.")
		print fstr % (str(num_ints), "gene-gene interactions loaded.")
		if level == 'feature':
			print fstr % (str(num_features), "gene-feature pairs loaded.")
		print "--------------------------------------------------------\n"
		
		return cls(h5_file, level, num_genes, num_ints, num_features, genes, gene_scores_wheel, max_gene_score, max_int_score)
	
	
	def get_tables(self):
		"""
		
		"""
		
		fh = tables.openFile(self.file_name, mode='r')
		
		node_table = fh.root.gene_network.genes
		edge_table = fh.root.gene_network.interactions
		feature_table = fh.root.gene_network.features
		
		return fh, node_table, edge_table, feature_table
	
	
	def get_graph(self, edges=True):
		"""
		
		"""
		
		fh, node_table, edge_table, feature_table = self.get_tables()
		
		## Load network data into graph
		nodes = [(row['gene'], row['score']) for row in node_table.iterrows()]
		G = nx.Graph()
		for g, s in nodes:
			G.add_node(g, weight=s)
		
		if edges:
			edges = [(row['gene1'], row['gene2'], row['score']) for row in edge_table.iterrows()]
			for g1, g2, s in edges:
				G.add_edge(g1, g2, weight=s)
		
		fh.close()
		
		return G


def is_int(s):
	""" 
	Determines if the specified parameter is an integer, or is a 
	string that can be converted to an integer.
	"""
	 
	try:
		int(s)
		return True
	except ValueError:
		return False


def is_float(s):
	""" 
	Determines if the specified parameter is a floating point number, or is a 
	string that can be converted to a floating point number.
	"""
	
	try:
		float(s)
		return True
	except ValueError:
		return False


def table_format_check(hdf5_fh, filetype, level=None):
	""" 
	This function checks that the input (or results) file is in the correct format for the polyGA program.
	
	Function Parameters:
	
	hdf5_fh		=	An HDF5 file handle for either an input data file or a results file
	filetype	=	Specifies whether the HDF5 file is an input or results file ('input', 'results')
	
	"""
	
	if filetype == 'input':
		## Get the file root
		root_name = hdf5_fh.title
		if not 'GA Network Data File - PolyGA' in root_name:
			print "The HDF5 data file was not created by PolyGA, \nand it may not be formatted correctly!\n"
			return False, None, None, None, None
			
		## Get the file groups
		groups = [group for group in hdf5_fh.walkGroups('/')]
		if len(groups) != 2:
			print "The HDF5 data file is not formatted correctly! Number of groups doesn't match!\n"
			return False, None, None, None, None
		if groups[1]._v_name != 'gene_network':
			print "The HDF5 data file is not formatted correctly! Group name doesn't match!\n"
			return False, None, None, None, None
		else:
			## Check that number of tables match
			nodes = hdf5_fh.listNodes(groups[1])
			if len(nodes) < 2:
				print "The HDF5 data file is not formatted correctly! Number of nodes doesn't match!\n"
				return False, None, None, None, None
			
			## Check that table names match
			try:
				gene_node = hdf5_fh.getNode(groups[1], name='genes', classname='Table')
			except:
				print "The HDF5 data file is not formatted correctly! Table names don't match!\n"
				return False, None, None, None, None
			else:
				level = gene_node.attrs.level
			try:
				int_node = hdf5_fh.getNode(groups[1], name='interactions', classname='Table')
			except:
				print "The HDF5 data file is not formatted correctly! Table names don't match!\n"
				return False, None, None, None, None
			if level == 'feature':
				try:
					feature_node = hdf5_fh.getNode(groups[1], name='features', classname='Table')
				except:
					print "The HDF5 data file is not formatted correctly! Table names don't match!\n"
					return False, None, None, None, None
			
			## Check column names and data types
			if gene_node.coltypes.get('gene') == 'string' and gene_node.coltypes.get('score') == 'float32':
				gene_rows = gene_node.nrows
			else:
				print "The HDF5 data file is not formatted correctly! Gene table data types don't match!\n"
				return False, None, None, None, None
			if int_node.coltypes.get('gene1') == 'string' and int_node.coltypes.get('gene2') == 'string' and int_node.coltypes.get('score') == 'float32':
				int_rows = int_node.nrows
			else:
				print "The HDF5 data file is not formatted correctly! Gene interactions table data types don't match!\n"
				return False, None, None, None, None
			if level == 'feature':
				if feature_node.coltypes.get('gene') == 'string' and feature_node.coltypes.get('feature') == 'string' and feature_node.coltypes.get('fscore') == 'float32':
					feature_rows = feature_node.nrows
				else:
					print "The HDF5 data file is not formatted correctly! Feature table data types don't match!\n"
					return False, None, None, None, None
			else:
				feature_rows = None
			return True, level, gene_rows, int_rows, feature_rows
	elif filetype == 'results':
		## Get the file root
		root_name = hdf5_fh.title
		
		# *** DEBUGGING / TESTING ***
		#print root_name
		
		if not 'GA Results File - PolyGA' in root_name:
			print "The HDF5 output file was not created by polyGA, \nand it may not be formatted correctly!\n"
			return False, None, None
		
		groups = [group for group in hdf5_fh.walkGroups('/')]
		if len(groups) != 2:
			print "The HDF5 results file is not formatted correctly! Number of groups doesn't match!\n"
			return False, None, None
		if groups[1]._v_name != 'ga':
			print "The HDF5 results file is not formatted correctly! Group name doesn't match!\n"
			return False, None, None
		else:
			## Check that table names match
			nodes = hdf5_fh.listNodes(groups[1])
			if len(nodes) != 1:
				print "The HDF5 results file is not formatted correctly! Number of nodes doesn't match!\n"
				return False, None, None
			try:
				ga_pop_node = hdf5_fh.getNode(groups[1], name='groups', classname='Table')
			except:
				print "The HDF5 results file is not formatted correctly! Table names don't match!\n"
				return False, None, None
			
			## Get info about the type of analysis
			level = ga_pop_node.attrs.level
			if ga_pop_node.attrs.connected == 'yes':
				connected = True
			else:
				connected = False
				
			## Check column names and data types
			if level == 'feature':
				if ga_pop_node.coltypes.get('generation') == 'int32' and ga_pop_node.coltypes.get('group') == 'int32' and ga_pop_node.coltypes.get('genes') == 'string' and ga_pop_node.coltypes.get('features') == 'string' and ga_pop_node.coltypes.get('size') == 'int32' and ga_pop_node.coltypes.get('fitness') == 'float32':
					ga_rows = ga_pop_node.nrows
				else:
					print "The HDF5 data file is not formatted correctly! GA table data types don't match!\n"
					return False, None, None
			else:
				if ga_pop_node.coltypes.get('generation') == 'int32' and ga_pop_node.coltypes.get('group') == 'int32' and ga_pop_node.coltypes.get('genes') == 'string' and ga_pop_node.coltypes.get('size') == 'int32' and ga_pop_node.coltypes.get('fitness') == 'float32':
					ga_rows = ga_pop_node.nrows
				else:
					print "The HDF5 data file is not formatted correctly! GA table data types don't match!\n"
					return False, None, None
			return True, level, ga_rows


def process_config(conf):
	"""
	
	"""
	## Read the GA parameters from the configuration file
	params = {}
	if not os.path.exists(conf):
		print "\nThe configuration file does not exist!\n"
		sys.exit()
	else:
		fh = open(conf, 'r')
		for line in fh:
			line = line.rstrip()
			if not (line == '' or line[0] == '#'):
				fields = line.split('=', 1)
				params[fields[0]] = fields[1]
		fh.close()
		
	## Perform checks on GA parameters from the configuration file, assign default values if necessary, and print to the screen the parameters used
	required_params = set(['pop_size', 'generations', 'min_group_size', 'max_group_size', 'select_type', 'migrants', 'cross_type', 'cross_prob', 'mut_prob', 'elite_num', 'node_weight', 'edge_weight', 'R_script', 'generation_restart', 'fitness_restart'])
	optional_params = set(['nprocs', 'hybrid_top', 'connected'])
	supported_params = required_params.union(optional_params)
	param_set = set(params.keys())
	missing_params = required_params.difference(param_set)
	
	print "Configuration file parameters:"
        print "--------------------------------------------------------"
	if len(missing_params) > 0:
		print "\nRequired parameters are missing from the configuration file: \n"
		for p in missing_params:
			print "%21s:  %s" % (p, 'Missing')
		print "\nPlease see the polyGA documentation for the correct format \nof the configuration file.\n"
		sys.exit(4)
	else:
		if params['select_type'] == 'hybrid':
			required_params.add('hybrid_top')
		missing_params = required_params.difference(param_set)
		if len(missing_params) > 0:
			print "\nRequired parameters are missing from the configuration file: \n"
                	for p in missing_params:
				print "%21s:  %s" % (p, 'Missing')
                	print "\nPlease see the polyGA documentation for the correct format \nof the configuration file.\n"
                	sys.exit(4)
		
		if not is_int(params['pop_size']) or int(params['pop_size']) < 1 or int(params['pop_size']) % 4 != 0:
			print "\nGA population size (pop_size) must be an integer and a multiple of 4.\n"
			sys.exit(4)
		if not is_int(params['generations']) or int(params['generations']) < 2:
			print "\nThe number of generations (generations) must be an integer greater than 2.\n"
			sys.exit(4)
		if not is_int(params['min_group_size']) or not is_int(params['max_group_size']):
			print "\nThe group size parameters (min_group_size, max_group_size) must be integers.\n"
			sys.exit(4)
		if int(params['min_group_size']) < 2 or int(params['min_group_size']) > int(params['max_group_size']):
			print "\nThe minimum group size must be >= 2, and must be <= the maximum group size.\n"
			sys.exit(4)
		if not params['select_type'] in ('truncate', 'roulette', 'hybrid'):
			print "\nThe selection type must be one of the following: 'truncate', 'roulette', 'hybrid'.\n"
			sys.exit(4)
		if params['select_type'] == 'hybrid':
			if not is_int(params['hybrid_top']) or (int(params['migrants']) + int(params['hybrid_top'])) >=  (int(params['pop_size'])/2):
				print "\nThe number of fittest groups (hybrid_top) to select with the \n'hybrid' selection method must be an integer. Also, the sum of \nfittest groups and random migrants (hybrid_top + migrants) \nmust be less than half the population size (pop_size/2).\n"
				sys.exit(4)
		if not params['cross_type'] in ('uniform'):
			print "\nOnly 'uniform' crossover (cross_type) is supported at this time.\n"
			sys.exit(4)
		if not is_float(params['cross_prob']) or float(params['cross_prob']) > 1 or float(params['cross_prob']) < 0:
			print "\nThe crossover probability (cross_prob) must be a floating \npoint number between 0 and 1.\n"
			sys.exit(4)
		if not is_float(params['mut_prob']) or float(params['mut_prob']) > 1 or float(params['mut_prob']) < 0:
			print "\nThe mutation probability (mut_prob) must be a floating \npoint number between 0 and 1.\n"
			sys.exit(4)
		if not is_int(params['elite_num']) or int(params['elite_num']) > (int(params['pop_size'])/2):
			print "\nThe number of 'elite' groups must be an integer less \nthan or equal to half the population size (pop_size/2).\n"
		if not is_float(params['node_weight']) or not is_float(params['edge_weight']) or float(params['node_weight']) < 0 or float(params['edge_weight']) < 0:
			print "\nThe node and edge weights (node_weight, edge_weight) \nmust be positive floating point numbers.\n"
			sys.exit(4)
		if float(params['node_weight'])+float(params['edge_weight']) != 1:
			print "\nThe node and edge weights (node_weight, edge_weight) \nmust sum to 1.\n"
			sys.exit(4)
		if 'nprocs' in params and not is_int(params['nprocs']):
			print "\nThe number of parallel processors (nprocs) must be an integer.\n"
			sys.exit(4)
		if 'parallel_nodes' in params:
			if not os.path.exists(params['parallel_nodes']):
				print "\nThe node file ("+params['parallel_nodes']+") does not exist!\n"
				sys.exit(4)
			else:
				with open(params['parallel_nodes'], 'r') as fh:
					for line in fh:
						n = line.rstrip()
					if n in nodes:
						nodes[n] = nodes[n] + 1
					else:
						nodes[n] = 1
				if len(nodes.keys()) > 1:
					pp_servers = tuple(nodes.keys())
				else:
					pp_servers = ()
			params['parallel_nodes'] = pp_servers
		if 'R_script' in params:
			if not os.path.exists(params['R_script']):
				print "\nThe R script ("+params['R_script']+") does not exist!\n"
				sys.exit(4)
		if not is_int(params['migrants']) or int(params['migrants']) >= (int(params['pop_size'])/2) or int(params['migrants']) < 0:
			print "\nThe number of random migrants to create each generation \n(num_migrants) must be an integer greater than 0 and \nless than half the population size (pop_size/2).\n"
			sys.exit(4)
		if not is_int(params['generation_restart']):
			print "\nThe restart threshold (generation_restart) must be an integer."
			print "For more information see the polyGA documentation.\n"
			sys.exit(4)
		if int(params['generation_restart']) < 1:
			params['generation_restart'] = None
		if not is_float(params['fitness_restart']) or float(params['fitness_restart']) >= 1:
			print "\nThe fitness restart threshold (fitness_restart) must be a floating \npoint number less than 1."
			print "For more information see the polyGA documentation.\n"
			sys.exit(4)
		if float(params['fitness_restart']) <= 0:
			params['fitness_restart'] = None
		if 'connected' in params:
			if params['connected'] in ['True', 'true', 'TRUE', 'T']:
				params['connected'] = True
			else:
				params['connected'] = False
		else:
			params['connected'] = False
	
	params_order = {'pop_size':1, 'generations':2, 'connected':3, 'min_group_size':4, 'max_group_size':4, 'select_type':5, 'hybrid_top':6, 'migrants':7, 'cross_type':8, 'cross_prob':9, 'mut_prob':10, 'elite_num':11, 'node_weight':12, 'edge_weight':13, 'R_script':14, 'generation_restart':15, 'fitness_restart':16, 'nprocs':17}
	extra_params = param_set.difference(supported_params)
	for key in sorted(params.keys(), key=lambda x: params_order.get(x)):
		print "%21s:  %s" % (key, params[key])
	
	if len(extra_params) > 0:
		print "\nThe following parameters will not be used: "
		for p in sorted(list(extra_params)):
			print "%21s:  %s" % (p, 'Not Used')
	print "--------------------------------------------------------\n"
	sys.stdout.flush()
	
	return params


class OptionsError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		print self.value
		return None


def process_options(args):
	"""
	
	"""
	
	# *** DEBUGGING / TESTING ***
	#print args
	
	cmd = None
	## Perform checks on the parameters specified at the command-line
	if args['results'] == None and args['data'] == None and (args['interactions'] == None or args['genes'] == None):
		## Print error message
		print "A results file, an HDF5 data file, or tab-delimited data \nfiles must be specified as input. For a gene-level \nanalysis, gene-gene interactions as well as gene scores \nmust be supplied. Additionally, a gene-feature map must \nbe supplied for a feature-level analysis.\n"
		raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.")
	elif args['network']:
		if args['results'] == None or args['data'] == None or args['generation'] == None or args['group'] == None:
			## Print error message
			print "A results file, a data file, a generation and a group number \nmust be specified as input to produce a network plot. \n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.")
		else:
			## Check that results file, data file, and threshold exist
			if not os.path.exists(args['results']):
				print "Results file does not exist!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif not os.path.exists(args['data']):
				print "Data file does not exist!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['generation'] == None:
				print "A generation must be specified.\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['group'] == None:
				print "A group number must be specified.\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			else:
				## Print message with parameters to screen
				print "Command-line parameters: "
				print "--------------------------------------------------------"
				print "%22s %s" % ("Data file: ", args['data'])
				print "%22s %s" % ("Results file: ", args['results'])
				print "%22s %s" % ("Level of analysis: ", args['analysis'])
				print "%22s %s" % ("Network plot: ", 'Yes')
				print "%22s %s" % ("Generation number: ", args['generation'])
				print "%22s %s" % ("Group number: ", args['group'])					
				print "--------------------------------------------------------\n"
				
				## Call the function to export the top association results
				cmd = 'network_plot'
	elif args['cytoscape']:
		if args['results'] == None or args['data'] == None or args['threshold'] == None:
			## Print error message
			print "A results file, a data file, and a fitness threshold \nmust be specified as input to produce cytoscape network files. \n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.")
		else:
			## Check that results file, data file, and threshold exist
			if not os.path.exists(args['results']):
				print "Results file does not exist!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif not os.path.exists(args['data']):
				print "Data file does not exist!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['threshold'] == None:
				print "You must specify a p-value threshold!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['out'] == None:
				print "An output file prefix must be specified.\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			else:
				## Print message with parameters to screen
				print "Command-line parameters: "
				print "--------------------------------------------------------"
				print "%23s %s" % ("Data file: ", args['data'])
				print "%23s %s" % ("Results file: ", args['results'])
				print "%23s %s" % ("Level of analysis: ", args['analysis'])
				print "%23s %s" % ("Cytoscape output: ", 'Yes')
				print "%23s %s" % ("Output file prefix: ", args['out'])
				print "%23s %s" % ("Fitness threshold: ", args['threshold'])					
				print "--------------------------------------------------------\n"
				
				## Call the function to export the top association results
				cmd = 'cytoscape'
	elif args['kegg']:
		if args['results'] == None or args['threshold'] == None:
			## Print error message
			print "A results file and a fitness threshold must be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		else:
			if not os.path.exists(args['results']):
				print "Results file does not exist!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['threshold'] == None:
				print "You must specify a p-value threshold!\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			elif args['out'] == None:
				print "An output file prefix must be specified.\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			else:
				## Print message with parameters to screen
				print "Command-line parameters: "
				print "--------------------------------------------------------"
				print "%23s %s" % ("Results file: ", args['results'])
				print "%23s %s" % ("Level of analysis: ", args['analysis'])
				print "%23s %s" % ("KEGG output: ", 'Yes')
				print "%23s %s" % ("Output file prefix: ", args['out'])
				print "%23s %s" % ("Fitness threshold: ", args['threshold'])					
				print "--------------------------------------------------------\n"
				
				## Call the function to export the top association results
				cmd = 'kegg'
	elif args['results']:
		## Check that results file exists
		if not os.path.exists(args['results']):
			print "Results file does not exist!\n"
			sys.exit(1)
		else:
			## Check the p-value threshold
			default_thresh = False
			if args['threshold'] == None:
				default_thresh = True
				args['threshold'] = 0.00005
			
			## Print message with parameters to screen
			print "Command-line parameters: "
			print "--------------------------------------------------------"
			print "%22s %s" % ("Results file: ", args['results'])
			print "%22s %s" % ("Level of analysis: ", args['analysis'])
			if default_thresh:
				print "%22s %s (default)" % ("Fitness threshold: ", args['threshold'])
			else:
				print "%22s %s" % ("Fitness threshold: ", args['threshold'])
			if args['plot']:
				print "%22s %s" % ("Create plot: ", 'Yes')
			else:
				print "%22s %s" % ("Create plot: ", 'No')
			print "--------------------------------------------------------\n"
			
			## Call the function to export the top association results
			if args['plot']:
				cmd = 'performance_plot'
			else:
				cmd = 'results'
	elif args['data']:
		## Check that input data file exists
		if not os.path.exists(args['data']):
			raise OptionsError("HDF5 file containing the gene \ network annotations does not exist!\n")
		## Check the level of analysis
		if args['analysis'] == None or not args['analysis'] in ['feature', 'gene']:
			print "The level of analysis was not specified or was not \nrecognized. Please indicate the type of search to \nperform ('feature' or 'gene').\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		## Check that GA configuration file exists
		if args['conf'] == None:
			print "A configuration file must be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if not os.path.exists(args['conf']):
			raise OptionsError("Configuration file does not exist!\n")
		## Check that the output file exists
		if args['out'] == None:
			print "An output file prefix must be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if os.path.exists(args['out']+'_out.h5'):
			raise OptionsError("Output file '"+args['out']+"_out.h5' already exists!")
		
		print "Command-line parameters: "
		print "--------------------------------------------------------"
		print "%27s %s" % ("Input data (HDF5) file: ", args['data'])
		print "%27s %s" % ("Configuration file: ", args['conf'])
		print "%27s %s" % ("Level of analysis: ", args['analysis'])
		print "%27s %s" % ("Output file prefix: ", args['out'])
		print "--------------------------------------------------------\n"
		
		## Start the GA
		cmd = 'ga_main'
	else:
		## Check that configuration file exists
		if args['conf'] == None:
			print "A configuration file must be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if not os.path.exists(args['conf']):
			raise OptionsError("Configuration file does not exist!\n")
		## Check the level of analysis
			if args['analysis'] == None or not args['analysis'] in ['feature', 'gene']:
				print "The level of analysis was not specified or was not \nrecognized. Please indicate the type of search to \nperform ('feature' or 'gene').\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		## Check that tab-delimited input files exist
		if args['interactions'] == None:
			print "A tab-delimited file containing gene-gene interactions \nmust be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if not os.path.exists(args['interactions']):
			raise OptionsError("Gene-gene interaction file does not exist!\n")
		if args['analysis'] == 'feature':
			if args['features'] == None:
				print "A tab-delimited file containing gene-feature mappings \nmust be specified.\n"
				raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
			if not os.path.exists(args['features']):
				raise OptionsError("Feature annotation file does not exist!\n")
		else:
			if args['features'] != None:
				if not os.path.exists(args['features']):
					args['features'] = None
		if args['genes'] == None:
			print "A tab-delimited file containing gene scores must \nbe specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if not os.path.exists(args['genes']):
			raise OptionsError("Gene annotation file does not exist!\n")
		## Check that the output file exists
		if args['out'] == None:
			print "An output file prefix must be specified.\n"
			raise OptionsError("Use the '-h' option for more info, and please see the \npolyGA documentation for details.\n")
		if os.path.exists(args['out']+'.h5'):
			raise OptionsError("Output file '"+args['out']+".h5' already exists!\n")
		elif os.path.exists(args['out']+'_out.h5'):
			raise OptionsError("Output file '"+args['out']+"_out.h5' already exists!\n")
		else:
			if args['analysis'] == 'feature':
				fstr = "%28s %s"
			else:
				fstr = "%26s %s"
			print "Command-line parameters: "
			print "--------------------------------------------------------"
			print fstr % ("Gene annotation file: ", args['genes'])
			print fstr % ("Gene interaction file: ", args['interactions'])
			if args['analysis'] == 'feature':
				print fstr % ("Feature annotation file: ", args['features'])
			print fstr % ("Configuration file: ", args['conf'])
			print fstr % ("Level of analysis: ", args['analysis'])
			print fstr % ("Output file prefix: ", args['out'])
			print "--------------------------------------------------------\n"
			
			## Start the GA
			cmd = 'ga_main'
	return cmd, args


def results_table(results, level, threshold=0.00005):
	"""
	Exports a tab-delimited file containing the top association results (with p-value <= threshold).
	
	Function Parameters:
	
	results		=	The filename of the results file
	level		=	The level of the analysis that was performed ('gene' or 'snp')
	threshold	=	The p-value threshold
	"""
	
	if not os.path.exists(results):
		raise OptionsError("Results file does not exist!\n")
	else:
		h5_ga = tables.openFile(results, mode='r')
		
		## Check that these tables exist, otherwise print a file format error message
		correct_format, results_level, ga_rows = table_format_check(h5_ga, 'results')
		if not correct_format:
			h5_ga.close()
			raise DataFormatError("Results file did not load successfully!\nThe file may be corrupted or improperly formatted.\n")
		else:
			if level == 'feature' and results_level == 'gene':
				print("\nWARNING: The results file is from a gene-level analysis!")
				print("Only gene groups will be printed to the output file.\n\n")
			fstr = '%' + str(len(str(ga_rows))+3) + 's' + ' %s'
			print "Loading results from HDF5 file ..."
			print "--------------------------------------------------------"
			if results_level == 'feature':
				print fstr % (str(ga_rows), "feature groups loaded.")
			else:
				print fstr % (str(ga_rows), "gene groups loaded.")
			print "--------------------------------------------------------\n"
				
		prefix = '.'.join(results.split('.')[:-1])
		assoc_file = prefix+'.assoc'
		## Open the file, first checking that it doesn't exist
		if os.path.exists(assoc_file):
			i = 2
			while os.path.exists(assoc_file):
				assoc_file = prefix+'.'+str(i)+'.assoc'
				i = i + 1
		
		## Open output file and write the header line
		fh = open(assoc_file, 'w')
		if results_level == 'feature':
			fh.write("generation\tgroup\tgenes\tfeatures\tp-value (<= "+str(threshold)+")\n")
		else:
			fh.write("generation\tgroup\tgenes\tp-value (<= "+str(threshold)+")\n")
		
		## Access the HDF5 tables
		ga_table = h5_ga.root.ga.groups
		if results_level == 'feature':
			feature_groups = [[row['generation'], row['group'], row['genes'],row['features'],row['fitness']] for row in ga_table.where("""fitness <= threshold""")]
			feature_groups.sort(key=lambda x: float(x[0]))
			unique_feature_groups = []
			top_hits = []
			for idx, feature_group in enumerate(feature_groups):
				if feature_group[1].shape[0] > 1:
					feature_list = [f for pair in feature_group[1].tolist() for f in pair if f != 'NA']
					feature_list =  list(set(feature_list))
				else:
					feature_list = [f for f in feature_group[1].tolist()[0] if f != 'NA']
				feature_list.sort()
				if not feature_list in unique_feature_groups:
					unique_feature_groups.append(feature_list)
					top_hits.append(feature_group)
				
			# *** DEBUGGING / TESTING ***
			#print len(top_hits)		
			
			top_hits.sort(key=lambda x: float(x[4]))
			for feature_group in top_hits:
				if feature_group[3].shape[0] > 1:
					feature_pairs = feature_group[3].tolist()
					gene_pairs = feature_group[2].tolist()
					genes = [g for pair in gene_pairs for g in pair if g != 'NA']
					features = [f for pair in feature_pairs for f in pair if f != 'NA']
					keep_features = []
					keep_indices = []
					for i, f in enumerate(features):
						if f not in keep_features:
							keep_features.append(f)
							keep_indices.append(i)
					keep_genes = [genes[j] for j in keep_indices]
					genes = keep_genes
					features = keep_features
				else:
					genes = feature_group[2].tolist()[0]
					features = feature_group[3].tolist()[0]
				fitness = feature_group[4]
				gene_str = ','.join([g for g in genes if g != 'NA'])
				feature_str = ','.join([f for f in features if f != 'NA'])
				fh.write(str(gene_group[0])+'\t'+str(gene_group[1])+'\t'+gene_str+'\t'+feature_str+'\t'+str(fitness)+'\n')
		else:
			gene_groups = [[row['generation'],row['group'],row['genes'],row['fitness']] for row in ga_table.where("""fitness <= threshold""")]
			unique_gene_groups = []
			top_hits = []
			for idx, gene_group in enumerate(gene_groups):
				if gene_group[2].shape[0] > 1:
					gene_list = [g for pair in gene_group[2].tolist() for g in pair if g != 'NA']
					gene_list =  list(set(gene_list))
				else:
					gene_list = [g for g in gene_group[2].tolist()[0] if g != 'NA']
				gene_list.sort()
				if not gene_list in unique_gene_groups:
					unique_gene_groups.append(gene_list)
					top_hits.append(gene_group)
			
			# *** DEBUGGING / TESTING ***
			#print len(top_hits)		
			
			top_hits.sort(key=lambda x: float(x[3]))
			for gene_group in top_hits:
				if gene_group[2].shape[0] > 1:
					gene_pairs = gene_group[2].tolist()
					genes = [g for pair in gene_pairs for g in pair if g != 'NA']
					keep_genes = []
					for g in genes:
						if g not in keep_genes:
							keep_genes.append(g)
					genes = keep_genes
				else:
					genes = gene_group[2].tolist()[0]
				fitness = gene_group[3]
				gene_str = ','.join([g for g in genes if g != 'NA'])
				fh.write(str(gene_group[0])+'\t'+str(gene_group[1])+'\t'+gene_str+'\t'+str(fitness)+'\n')
			
		print "Significant associations printed to file ("+assoc_file+").\n"
		fh.close()
		h5_ga.close()
		return None


def performance_plot(results, level, filetype='pdf'):
	"""
	Exports a plot of the GA performance (fitness vs. generation).
	
	Function Parameters:
	
	results		=	The filename of the results file
	filetype	=	The image format to be exported (any format supported by the Matplotlib savefig() method)
	"""
	
	try:
		import matplotlib
	except ImportError as imperr:
		print "Warning: The matplotlib module is required to create a \nplot. No plot will be generated.\n"
		print str(imperr)+"\n"
		return None
	else:
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		from matplotlib.font_manager import FontProperties
	
	if not os.path.exists(results):
		raise OptionsError("Results file does not exist!\n")
	else:
		h5_ga = tables.openFile(results, mode='r')
		
		## Check the table format
		correct_format, results_level, ga_rows = table_format_check(h5_ga, 'results')
		if not correct_format:
			raise DataFormatError("Results did not load successfully!\n")
		
		prefix = '.'.join(results.split('.')[:-1])
		plot_file = prefix+'.plot.'+filetype
		## Open the file, first checking that it doesn't exist
		if os.path.exists(plot_file):
			i = 2
			while os.path.exists(plot_file):
				plot_file = prefix+'.plot'+str(i)+'.'+filetype
				i = i + 1
		
		ga_table = h5_ga.root.ga.groups
		max_gen = max(ga_table.cols.generation)
		
		## Get the average and best fitness for each generation
		mean_fitness = []
		best_fitness = []
		if max_gen > 100:
			step = max_gen / 100
			gens = range(0, max_gen+1, step)
			gens[0] = 1
		else:
			gens = range(1, max_gen+1)
		for gen in gens:
			fit_list = [-math.log10(row['fitness']) for row in ga_table.where("""generation == gen""")]
			mean_fitness.append(float(sum(fit_list))/len(fit_list))
			best_fitness.append(max(fit_list))
			#print str(gen)+": "+str(len(pval_list))
		h5_ga.close()
		
		## Set font size
		font = {'size':10}
		matplotlib.rc('font', **font)
		## Create plot 
		fig = plt.figure(figsize=(5, 4))
		## Create an Axes object
		ax = fig.add_subplot(1,1,1)
		## Plot the data
		ax.plot(gens, mean_fitness, 'g.--', label="Average Fitness")
		ax.plot(gens, best_fitness, 'b.--', label="Best Fitness")
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height*0.15, box.width, box.height*0.85])
		ax.legend(loc=9, bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=2, fontsize=10)
		ax.set_title("PolyGA Performance", fontsize=10)
		ax.set_xlabel("Generation")
		ax.set_ylabel("-log10(fitness)")
		## Save the figure
		fig.savefig(plot_file, format=filetype)
		print "Plot of GA performance saved to file ("+plot_file+").\n"
		return None


def cytoscape_output(data, results, level, gen, group_id, threshold):
	"""
	Creates a text file that can be used to create a network in Cytoscape version 2.4 or later.
	
	Function Parameters:
	
	data		=	The filename of the data file
	results		=	The filename of the results file
	level		=	The level of the analysis that was performed ('gene' or 'snp')
	gen			=	The generation of the sub-network to output
	group_id	=	The group ID of the sub-network to output
	
	"""
	
	try:
		import networkx as nx
	except ImportError as imperr:
		print "Warning: The networkx module is required.\n"
		print str(imperr)+"\n"
		return None
	
	if not os.path.exists(results) or not os.path.exists(data):
		raise OptionsError("One of the input files does not exist!\n")
	else:
		h5_ga = tables.openFile(results, mode='r')
		
		## Check that these tables exist, otherwise print a file format error message
		correct_format, results_level, ga_rows = table_format_check(h5_ga, 'results')
		if not correct_format:
			h5_ga.close()
			raise DataFormatError("Results file did not load successfully!\nThe file may be corrupted or improperly formatted.\n")
		else:
			if level == 'feature' and results_level == 'gene':
				print("\nWARNING: The results file is from a gene-level analysis!")
				print("Only gene groups will be printed to the output file.\n\n")
			fstr = '%' + str(len(str(ga_rows))+3) + 's' + ' %s'
			print "Loading results from HDF5 file ..."
			print "--------------------------------------------------------"
			if results_level == 'feature':
				print fstr % (str(ga_rows), "feature groups loaded.")
			else:
				print fstr % (str(ga_rows), "gene groups loaded.")
			print "--------------------------------------------------------\n"
		
		h5_data = tables.openFile(data, mode='r')
		correct_format, level, gene_rows, int_rows, feature_rows = table_format_check(h5_data, 'input')
		if not correct_format:
			raise DataFormatError("Data file did not load successfully!\n")
		else:
			gene_table = h5_data.root.gene_network.genes
			int_table = h5_data.root.gene_network.interactions
			
			nodes = [(row['gene'], row['score']) for row in node_table.iterrows()]
			edges = [(row['gene1'], row['gene2'], row['score']) for row in edge_table.iterrows()]
			h5_data.close()
			
			## Create NetworkX graph
			G = nx.Graph()
			for g, s in nodes:
				G.add_node(g, weight=s)
			for g1, g2, s in edges:
				G.add_edge(g1, g2, weight=s)
			
			## Print data loading message
			if level == 'feature':
				max_num = max(len(str(gene_rows)), len(str(int_rows)), len(str(feature_rows)))
			else:
				max_num = max(len(str(gene_rows)), len(str(int_rows)))
			print "Loading data from HDF5 file ..."
			print "--------------------------------------------------------"
			fstr = '%' + str(max_num+3) + 's' + ' %s'
			print fstr % (str(gene_rows), "gene annotations loaded.")
			print fstr % (str(int_rows), "gene-gene interactions loaded.")
			if level == 'feature':
				print fstr % (str(feature_rows), "gene-feature pairs loaded.")
			print "--------------------------------------------------------\n"
		
		prefix = '.'.join(results.split('.')[:-1])
		out_file = prefix+'.cyto'+'.txt'
		## Open the file, first checking that it doesn't exist
		if os.path.exists(plot_file):
			i = 2
			while os.path.exists(plot_file):
				out_file = prefix+'.cyto'+str(i)+'.txt'
				i = i + 1
		
		## Access the HDF5 tables
		ga_table = h5_ga.root.ga.groups
		
		if gen != 0 and group_id != 0:
			gene_groups = [[row['genes']] for row in ga_table.where("""(generation == gen) & (group == group_id)""")]
		elif threshold != False:
			gene_groups = [[row['genes']] for row in ga_table.where("""fitness <= threshold""")]
		else:
			print "Error: No groups specified!\n"
			sys.exit(1)
			
		if len(gene_groups) == 0:
			print "No groups found!\n"
			sys.exit(1)
			#elif len(gene_groups) > 1:
			#print "More than one group found!\n"
			#sys.exit(1)
		else:
			gene_pair_list = []
			for gene_group in gene_groups:
				gene_pair_list.extend([tuple(pair) for pair in gene_group[0].tolist() if pair != ['NA', 'NA']])
			
			## Create NetworkX graph
			G2 = nx.Graph()
			G2.add_edges_from(gene_pair_list)
			h5_ga.close()
			
			## Get node weights
			for g in G2.nodes():
				G2.node[g]['weight'] = G.node[g]['weight']
			
			## Get edge weights
			for g1, g2 in G2.edges():
				G2.edge[g1][g2]['weight'] = G.edge[g1][g2]['weight']
			
			## Write to file
			fh = open(out_file, 'w')
			fh.write("nodeA\tinteraction\tnodeB\tedge_weight\tnodeA_weight\tnodeB_weight\n")
			for g1, g2 in G2.edges():
				fh.write(g1 + "\tpp\t" + g2 + "\t" + str(G2.edge[g1][g2]['weight'])+"\t" + str(G2.node[g1]['weight']) + "\t" + str(G2.node[g2]['weight']) + "\n")
			fh.close()
	
	return None


def network_plot(data, results, level, gen, group_id, title=None, filetype='pdf'):
	"""
	Exports a plot of a subnetwork.
	
	Function Parameters:
	
	data		=	The filename of the data file
	results		=	The filename of the results file
	level		=	The level of the analysis that was performed ('gene' or 'snp')
	gen			=	The generation of the sub-network to draw
	group_id	=	The group ID of the sub-network to draw
	title		=	Plot title
	filetype	=	Graphics file format
	"""
	
	try:
		import networkx as nx
		import matplotlib
	except ImportError as imperr:
		print "Warning: The matplotlib and networkx modules are required \nto create a plot. No plot will be generated.\n"
		print str(imperr)+"\n"
		return None
	else:
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		from matplotlib.font_manager import FontProperties
	
	if not os.path.exists(results) or not os.path.exists(data):
		raise OptionsError("One of the input files does not exist!\n")
	else:
		h5_ga = tables.openFile(results, mode='r')
		
		## Check that these tables exist, otherwise print a file format error message
		correct_format, results_level, ga_rows = table_format_check(h5_ga, 'results')
		if not correct_format:
			h5_ga.close()
			raise DataFormatError("Results file did not load successfully!\nThe file may be corrupted or improperly formatted.\n")
		else:
			if level == 'feature' and results_level == 'gene':
				print("\nWARNING: The results file is from a gene-level analysis!")
				print("Only gene groups will be printed to the output file.\n\n")
			fstr = '%' + str(len(str(ga_rows))+3) + 's' + ' %s'
			print "Loading results from HDF5 file ..."
			print "--------------------------------------------------------"
			if results_level == 'feature':
				print fstr % (str(ga_rows), "feature groups loaded.")
			else:
				print fstr % (str(ga_rows), "gene groups loaded.")
			print "--------------------------------------------------------\n"
		
		h5_data = tables.openFile(data, mode='r')
		correct_format, level, feature_rows, int_rows = table_format_check(h5_data, 'input')
		if not correct_format:
			raise DataFormatError("Data file did not load successfully!\n")
		else:
			gene_table = h5_data.root.gene_network.genes
			int_table = h5_data.root.gene_network.interactions
			
			nodes = [(row['gene'], row['score']) for row in node_table.iterrows()]
			edges = [(row['gene1'], row['gene2'], row['score']) for row in edge_table.iterrows()]
			h5_data.close()
			
			## Create NetworkX graph
			G = nx.Graph()
			for g, s in nodes:
				G.add_node(g, weight=s)
			gene_rows = G.number_of_nodes()
			for g1, g2, s in edges:
				G.add_edge(g1, g2, weight=s)
			
			## Print data loading message
			if level == 'feature':
				max_num = max(len(str(gene_rows)), len(str(int_rows)), len(str(feature_rows)))
			else:
				max_num = max(len(str(gene_rows)), len(str(int_rows)))
			print "Loading data from HDF5 file ..."
			print "--------------------------------------------------------"
			fstr = '%' + str(max_num+3) + 's' + ' %s'
			print fstr % (str(gene_rows), "gene annotations loaded.")
			print fstr % (str(int_rows), "gene-gene interactions loaded.")
			if level == 'feature':
				print fstr % (str(feature_rows), "gene-feature pairs loaded.")
			print "--------------------------------------------------------\n"
		
		prefix = '.'.join(results.split('.')[:-1])
		plot_file = prefix+'.network.'+filetype
		## Open the file, first checking that it doesn't exist
		if os.path.exists(plot_file):
			i = 2
			while os.path.exists(plot_file):
				plot_file = prefix+'.network'+str(i)+'.'+filetype
				i = i + 1
		
		## Access the HDF5 tables
		ga_table = h5_ga.root.ga.groups
		if results_level == 'feature':
			pass
		else:
			gene_groups = [[row['genes']] for row in ga_table.where("""(generation == gen) & (group == group_id)""")]
			if len(gene_groups) == 0:
				print "No groups found"
				sys.exit()
			elif len(gene_groups) > 1:
				print "More than one group found!"
				sys.exit()
			else:
				gene_group = gene_groups[0]
				gene_pair_list = [tuple(pair) for pair in gene_group[0].tolist() if pair != ['NA', 'NA']]
				
				## Create NetworkX graph
				G2 = nx.Graph()
				G2.add_edges_from(gene_pair_list)
				h5_ga.close()
				
				## Get node weights
				for g in G2.nodes():
					G2.node[g]['weight'] = G.node[g]['weight']
				
				## Get edge weights
				for g1, g2 in G2.edges():
					G2.edge[g1][g2]['weight'] = G.edge[g1][g2]['weight']
				
				## Draw figure
				node_list = G2.nodes()
				edge_list = G2.edges()
				node_size = []
				edge_width = []
				for g in node_list:
					node_size.append(G2.node[g]['weight'])
				for g1, g2 in edge_list:
					edge_width.append(G2.edge[g1][g2]['weight'])
				max_node = max(node_size)
				max_edge = max(edge_width)
				node_size = [300*(float(w)/max_node) for w in node_size]
				edge_width = [10*(float(w)/max_edge) for w in edge_width]
				
				font = {'size':10}
				matplotlib.rc('font', **font)
				fig = plt.figure(figsize=(5,4))
				nx.draw(G, font_size=6, nodelist=node_list, node_size=node_size, edgelist=edge_list, width=edge_width, node_color='#A0CBE2', edge_color='#424242')
				if title != None:
					plot_title = fig.suptitle(title, fontsize=10)
				fig.savefig(plot_file, format=filetype)
	
	return None

