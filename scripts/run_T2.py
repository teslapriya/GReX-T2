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


def test_T2(mode, outputfile="T2_output.txt", candsfile="giants.cand.", host="127.0.0.1", port=12345, plot=True, plot_dir="./"):
    """
    Test T2 functions. Read in heimdall giants, cluster them, send results to outputfile and save plots. 

    Parameters
    ----------
    mode : String. "file" or "socket". 
    outputfile : default is "T2_output.txt".
    candsfile : default is "giants.cand.".
    host, port : use the same host, port with heimdall -coincidencer HOST:PORT 
    plot : default is True.
    plot_dir : default is "./".
    """
    print("========= begin one T2 test: ", mode, "mode ==============")
    # read in giants 
    if mode == "file":
        tab, data, snrs = T2.cluster_heimdall.parse_candsfile(candsfile, selectcols=['itime', 'dm', 'ibox'])
    elif mode == "socket":
        # a better place for socket trails: test_socket.py 
        tab, data, snrs = T2.cluster_heimdall.parse_socket("127.0.0.1", 12345) # use the same host, port with heimdall -coincidencer HOST:PORT 
    else: # can raise error 
        print("mode is either file or socket.") 
   
    # T2 cluster 
    clusterer, data_labeled = T2.cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10, metric='euclidean', allow_single_cluster=True)
    clsnr = T2.cluster_heimdall.get_peak(data_labeled, snrs) 
    
    # send T2 cluster results to outputfile
    T2.cluster_heimdall.dump_cluster_results(tab, clsnr, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm'])
    
    if plot: 
        T2.cluster_heimdall.plot_giants(tab, plot_dir=plot_dir) # plot giants      
        T2.cluster_heimdall.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir) # plot cluster results  
    
    print("one T2 test completed (", mode, "mode).\n")

      
    

        

if __name__ == '__main__':
    '''
    #my_candsfile = "giants_1.cand" # heimdall output of a 45s-long filterbank with 3 pulses 10s apart, DM600.
    #my_outputfile = "T2_output.txt" # dump data to a jason file. 
    #my_plot_dir = "plots" # save plots there. 
    
    # command line options and arguments 
    #main(sys.argv[1:])
    
    # python run_T2_v4.py "giants.cand" "T2_output.txt" "./plots" 
    #my_candsfile = sys.argv[0]
    #my_outputfile = sys.argv[1]
    #my_plot_dir = sys.argv[2]
    #my_cols, my_metric, my_min_cluster_size, my_min_samples, my_plot, my_save_dir = optional arguments 
    '''
    
    # test 1: file 
    test_T2("file", outputfile="T2_output_file.txt", candsfile="giants_1.cand", plot=True, plot_dir="file_")
    
    # test 2: socket, one gulp only 
    #test_T2("socket", outputfile="T2_output_socket.txt", host="127.0.0.1", port=12345, plot=True, plot_dir="socket_")
    
    # test 3: socket, continuously 
    T2.cluster_heimdall.parse_socket_and_cluster_and_plot(host="127.0.0.1", port=12345, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="T2_output_socket", plot=True, plot_dir="socket_")
    
    '''
    # read in and plot giants
    # from heimdall output file 
    tab, data, snrs = T2.cluster_heimdall.parse_candsfile(my_candsfile, selectcols=['itime', 'dm', 'ibox'])
   
    # from socket 
    # only one gulp?  
    tab, data, snrs = T2.cluster_heimdall.parse_socket("127.0.0.1", 12345) # host, port should be the same with heimdall. 
    print(tab)
    
      
    # plot giants 
    T2.cluster_heimdall.plot_giants(tab, plot_dir=my_plot_dir)  
    
    # cluster 
    clusterer = T2.cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10, metric='euclidean', returnclusterer=True, allow_single_cluster=True)
    data_labeled = T2.cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10, metric='euclidean', returnclusterer=False, allow_single_cluster=True)
    clsnr = T2.cluster_heimdall.get_peak(data_labeled, snrs) 
    T2.cluster_heimdall.dump_cluster_results(tab, clsnr, my_outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm'])
    
    # plot cluster results  
    T2.cluster_heimdall.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=my_plot_dir)
    '''
 

      




