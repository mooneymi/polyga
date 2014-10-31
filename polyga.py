#!/usr/bin/env python

"""

PolyGA - version 1.2b - Last modified: 30-DEC-2013

Usage:	
	python polyga.py -h
	python polyga.py -a feature -c <filename> -i <filename> -g <filename> -s <filename> -o <file prefix>
	python polyga.py -a feature -c <filename> -d <HDF5 filename> -o <file prefix>

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

def polyga_main(cmd, opts):
	"""
	
	"""
	
	if cmd == 'ga_main':
		params = polyga_utils.process_config(opts['conf'])
		polyga_core.ga_main(opts, params)
	elif cmd == 'results':
		polyga_utils.results_table(opts['results'], opts['analysis'], opts['threshold'])
	elif cmd == 'performance_plot':
		polyga_utils.results_table(opts['results'], opts['analysis'], opts['threshold'])
		polyga_utils.performance_plot(opts['results'], opts['analysis'])
	elif cmd == 'network_plot':
		polyga_utils.network_plot(opts['data'], opts['results'], opts['analysis'], opts['generation'], opts['group'])
	elif cmd == 'cytoscape':
		polyga_utils.cytoscape_output(opts['data'], opts['results'], opts['analysis'], opts['generation'], opts['group'], opts['threshold'])
	elif cmd == 'kegg':
		polyga_utils.kegg_path_plot(opts['results'], opts['analysis'], opts['threshold'], opts['kegg'], opts['out'])
	else:
		print "There was an error processing the command-line options!\n"
	
	return None


if __name__ == '__main__':
	try:
		import argparse
		import datetime
		import polyga_utils
		import polyga_core
	except ImportError as imperr:
		print "Error: one or more required modules failed to load.\n"
		raise ImportError(str(imperr))
	else:	
		parser = argparse.ArgumentParser(description="Performs a search to identify polygenic signatures associated with a trait of interest. Candidate feature sets are identified by searching a gene-gene interaction network using a Genetic Algorithm, which is guided by user-supplied expert knowledge.")
		parser.add_argument('-a', '--analysis', required=False, default='feature', choices=['feature', 'gene'], help="The level (either 'gene' or 'feature') at which to perform the search.")
		parser.add_argument('-c', '--conf', required=False, help="Location of the GA configuration file containing the algorithm parameters.")
		parser.add_argument('-d', '--data', required=False, help="Location of data in HDF5 format. This option assumes HDF5 format for output also.")
		parser.add_argument('-i', '--interactions', required=False, help="Location of the gene-gene interactions file.")
		parser.add_argument('-f', '--features', required=False, help="Location of the gene-feature map file.")
		parser.add_argument('-g', '--genes', required=False, help="Location of the gene annotation file.")
		parser.add_argument('-o', '--out', required=False, help="Location and prefix to save the data and results.")
		parser.add_argument('-r', '--results', required=False, help="Location of the polyGA results file. This option will produce a tab-delimited file containing the top hits identified by the algorithm (p-value less than the threshold specified with the '-t' option).")
		parser.add_argument('-t', '--threshold', required=False, type=float, help="The fitness threshold for reported feature groups (default = 5e-5). The results file must be specified with '-r'.")
		parser.add_argument('-p', '--plot', required=False, action='store_true', default=False, help="Will produce a plot (as a PDF) of the GA performance (-log10(fitness) vs. generation). The results file must be specified with '-r'.")
		parser.add_argument('-n', '--network', required=False, action='store_true', default=False, help="Will produce a network plot (as a PDF) of the specified feature group. Two integers, indicating the generation and group number of the group to be plotted, must be provided. Both a data file and a results file must also be specified with the '-d' and '-r' options.")
		parser.add_argument('generation', type=int, help="The generation of the feature group to be plotted. Required only for the network plot '-n' option.")
		parser.add_argument('group', type=int, help="The group number of the feature group to be plotted. Required only for the network plot '-n' option.")
		parser.add_argument('-cyto', '--cytoscape', required=False, action='store_true', default=False, help="Will produce files that can be imported into Cytoscape. This will allow visualization of significant variant interactions within the context of the gene-gene interaction network. Currently this option can only be used for analyses investigating variant groups of size 2. Results and data files must be specified with the '-r' and '-d' options. Also, an output file prefix and a fitness threshold must be specified with the '-o' and '-t' options.")
		parser.add_argument('-k', '--kegg', required=False, help="Location of the KEGG pathway XML file. This option will produce a file that can be used with the 'plotKEGGpathway.r' R script to produce a pathway figure annotated with association results. Currently this option can only be used for analyses investigating variant groups of size 2. A results file, an output file prefix and a fitness threshold must be specified with the '-r', '-o' and '-t' options.")
		args_obj = parser.parse_args()
		args = vars(args_obj)
		
		print "\npolyGA (version 1.0b)\n"
		start_time = datetime.datetime.now()
		print "Start date and time: "+start_time.strftime("%Y-%m-%d %H:%M")+"\n"
		
		cmd, opts = polyga_utils.process_options(args)
		
		polyga_main(cmd, opts)
		
		end_time = datetime.datetime.now()
		print "End date and time: "+end_time.strftime("%Y-%m-%d %H:%M")+"\n"
		elapsed_time = end_time - start_time
		print "Elapsed time: "+str(elapsed_time)+"\n"
