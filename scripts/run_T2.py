#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import T2 
import sys, getopt  
import argparse


def main(argv):
    """    
    Use getopt.getopt to parse command line options.
    Parameters
    ----------
    argv : TYPE
        DESCRIPTION.
    """
    inputfile = ''
    outputfile = ''
    plotdir = '' 
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile=", "plotdir="])
    except getopt.GetoptError:
       print('run_T2.py -i <inputfile> -o <outputfile> -p <saveplotdir>')
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print('run_T2.py -i <inputfile> -o <outputfile> -p <plotdir>')
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
       elif opt in ("-o", "--ofile"):
          outputfile = arg
       elif opt in ("-p", "--plotdir"):
          outputfile = arg
    print('Input file is "', inputfile)
    print('Output file is "', outputfile)
    print('Save plots at "', plotdir) 

if __name__ == '__main__':
    T2.socket.parse_socket(host="127.0.0.1", ports=[12345,12346], selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="tests/outputs/cluster_output", plot_dir="tests/outputs/")

