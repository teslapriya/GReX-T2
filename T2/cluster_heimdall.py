#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
from astropy.io import ascii
import numpy as np
import hdbscan
import json 
import seaborn as sns
import socket 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable



def parse_candsfile(candsfile, selectcols=['itime', 'idm', 'ibox', 'ibeam']):
    """ Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output
    (Can add cleaning here, eventually)
    """
    #tab = pd.read_csv(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], delim_whitespace=True) 
    tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], guess=False, fast_reader=False, delimiter="\s")
    tab['ibeam'] = tab['ibeam'].astype(int)
    data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
    snrs = tab['snr']
    # how to use ibeam?

    return tab, data, snrs


def cluster_data(data, min_cluster_size=3, min_samples=5, metric='hamming', returnclusterer=False, allow_single_cluster=True):
    """ Take data from parse_candsfile and identify clusters via hamming metric.
    """

    clusterer = hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_cluster_size,
                                min_samples=min_samples, cluster_selection_method='eom',
                                allow_single_cluster=allow_single_cluster).fit(data) 

    nclustered = np.max(clusterer.labels_ + 1) 
    nunclustered = len(np.where(clusterer.labels_ == -1)[0]) 
    print('Found {0} clustered and {1} unclustered rows'.format(nclustered, nunclustered))
    print('All labels: {0}'.format(np.unique(clusterer.labels_)))
    data_labeled = np.hstack((data, clusterer.labels_[:,None])) 

    if returnclusterer:
        return clusterer
    else:
        return data_labeled


def get_peak(data_labeled, snrs):
    """ Given labeled data, find max snr row per cluster
    """

    clsnr = []
    clusters = data_labeled[:, 3]
    for cluster in np.unique(clusters):
        clusterinds = np.where(cluster == clusters)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        clsnr.append((imaxsnr, maxsnr))

    return clsnr


def dump_cluster_results(tab, clsnr, outputfile, output_cols=['mjds', 'snr', 'ibox', 'dm']):    
    """   
    Takes tab from parse_candsfile and clsnr from get_peak, 
    output columns output_cols into a jason file outputfile. 
    """
    imaxsnr_arr = [clsnr[i][0] for i in range(len(clsnr))] 
    output_dict = {} 
    for i in imaxsnr_arr:
        output_dict[tab['mjds'][i]] = {} 
        for col in output_cols:
            if type(tab[col][i]) == np.int64:
                output_dict[tab['mjds'][i]][col] = int(tab[col][i]) 
            else: 
                output_dict[tab['mjds'][i]][col] = tab[col][i] 

    with open(outputfile, 'w') as f: #encoding='utf-8'
        json.dump(output_dict, f, ensure_ascii=False, indent=4) 
    

def plot_clustered(clusterer, clsnr, snrs, data, tab, cols, plot_dir="./"):
    """
    Plot cluster probabilities in color and greyness. 

    Parameters
    ----------
    clusterer : from cluster_data(). 
    clsnr: from get_peak(). 
    tab, data, snrs: from parse_candsfile(). 
    cols : parameters to plot. 
    plot_dir : where to save output plots. 
    """
    # grey for unclustered noise, pure color for well clustered points. 
    palette = sns.color_palette() 
    cluster_colors = [sns.desaturate(palette[col], sat) if col >= 0 else (0.5, 0.5, 0.5) for col, sat in zip(clusterer.labels_, clusterer.probabilities_)] 
    
    for i in range(len(cols)): 
        fig, ax = plt.subplots() 
        ax.cla() 
        ax.scatter(data[:,i], snrs, s=3, c=cluster_colors) 
        ax.set_xlabel(cols[i]) 
        ax.set_ylabel('snr') 
        ax.set_title('cluster cols:'+str(cols))
        fig.savefig(plot_dir+'snr_'+str(cols[i])+'.pdf') 
        fig.clf()
        plt.close('all')
        
        for j in range(len(cols)): 
            if j > i: 
                fig, ax = plt.subplots() 
                ax.cla()
                ax.scatter(data[:,i],data[:,j], c=cluster_colors)
                ax.set_xlabel(cols[i]) 
                ax.set_ylabel(cols[j]) 
                
                for k in range(len(clsnr)): 
                    imaxsnr = clsnr[k][0]
                    maxsnr = int(clsnr[k][1])
                    ax.scatter(data[:,i][imaxsnr],data[:,j][imaxsnr], s=maxsnr, c='k', marker='*') 
                    ax.text(data[:,i][imaxsnr],data[:,j][imaxsnr], str(maxsnr)) 
                
                ax.set_title('cluster cols:'+str(cols))
                fig.savefig(plot_dir+'cluster_prob_'+cols[i]+'_'+cols[j]+'.pdf') 
                fig.clf()
                plt.close('all') 





    
# Below are functions to plot giants.                       
def plot_dm_hist(tab, nbins=30, plot_dir="./", multibeam=False, data_name = None):
    """
    plot the giants DM histogram

    Parameters
    ----------
    tab : ascii table
        MBHeimdall giants output full table from parse_candsfile or parse_socket.
    nbins : int
        Bin number. The default is 30.
    plot_dir : optional.  Where to save plots. The default is "./".
    multibeam : bool, optional. The default is False.
        If True, plot a histogram for each beam.
    data_name : bool, optional. The default is None.

    Returns
    -------
    Save the plot in plot_dir. 
    """
    fig, ax = plt.subplots()
    ax.cla() #Clear the current axes
    dm_min = 1.
    #tab['ibeam'] = tab['ibeam'].astype(int)
    if multibeam:
        for beam in range(tab['ibeam'].max()+1): #due to zero based indexing
            if beam in tab['ibeam']:
                cands = tab[tab['ibeam'] == beam]
                logbins = np.logspace(np.log10(dm_min),np.log10(cands['dm'].max()),nbins+1)
                vals,edges = np.histogram(cands['dm'],bins=logbins)
                ax.step(edges,np.append(vals,0.),where='post',label=str(beam))
                ax.legend(loc=9,ncol=4,fontsize=8)
    else:
        logbins = np.logspace(np.log10(dm_min),np.log10(tab['dm'].max()),nbins+1)
        ax.hist(tab['dm'],bins=logbins,histtype='step', label=data_name)

    ### Set global plot window parameters
    ax.set_xlim(dm_min,tab['dm'].max())
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set ticks and grid
    ax.grid(axis='x', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    formatter = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    ax.set_xlabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
    ax.set_ylabel('$\\rm Giants count$', size=12)
    ax.set_title('giants dm')
    fig.savefig(plot_dir+'giants_dm_hist_multibeam_'+str(multibeam)+'.pdf') 
    fig.clf()
    plt.close('all')     



def plot_dm_snr(ax, ax_cbar,tab,tsamp=1048e-6):
    """
    

    Parameters
    ----------
    ax_cbar : TYPE
        DESCRIPTION.
    tab : ascii table
        MBHeimdall giants output full table from parse_candsfile or parse_socket.
    tsamp : optional. The default is 1048e-6.

    Returns
    -------
    None.

    """
    #fig, ax = plt.subplots()
    ax.cla() 
    ax_cbar.cla()

    max_snr = 100.

    ax.set_xlim(6.,max_snr)
    ax.set_xscale('log')
    ax.set_xlabel('$\\rm SNR$', size=12)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='y', labelleft=False)

    # Set ticks and grid
    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    xformat = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(xformat)

    xticks     = [6.,8.,10.,20.,40.,100.]
    xticklabel = ['%i' % tick for tick in xticks]
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticklabel)
    ax.xaxis.set_minor_locator(plt.NullLocator())

    if tab['snr'].max() > max_snr:
        print("\n    WARNING: The SNR of a candidate is higher than the maximum SNRs plotted!\n")
        
    colormap = ax.scatter(tab['snr'], tab['dm'], s=tab['ibox'], c=tab['ibox'], cmap='viridis', \
                   vmin=0., vmax=12, marker='o', picker=2, alpha=0.5)
    # Add colorbar
    try: colormap
    except NameError:
        ax_cbar.axis("off")
        pass

    ax_cbar.axis("on")
    cbar = plt.colorbar(colormap,cax=ax_cbar,use_gridspec=True,label='$\\rm Boxcar width\;(index)$')
    cticks = np.array(cbar.get_ticks())
    cbar.ax.set_yticklabels([x[0:5] for x in (2.**cticks * tsamp*1000.).astype('str')])
    cbar.set_alpha(1)
    cbar.draw_all()



def plot_time_dm(ax,tab,duration=None,snr_cut=8,snr_thr=25.,mps=30.,multibeam=False,axrange=True,axlabel=True):
    ax.cla() 
    max_point_size = mps #diameter in points
    
    if axrange:

        if duration != None:
            ax.set_xlim(0.,duration) #duration in seconds
        else:
            ax.set_xlim(0.,tab['mjds'].max()*1.05)
        
        ax.set_ylim(1.,tab['dm'].max()*1.05)
        ax.set_yscale('log')

    if axlabel == True:
        ax.set_xlabel('$\\rm Time\; (sec)$', size=12)
        ax.set_ylabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
    else:
        ax.tick_params(axis='x', labelbottom=False)


    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    yformat = plt.FormatStrFormatter('%i')
    ax.yaxis.set_major_formatter(yformat)

    ax.scatter(tab['mjds'],tab['dm'], \
               s=(max_point_size * tab['snr'] / tab['snr'].mean())**2, \
               c=tab['ibox'], cmap='viridis', \
               vmin=0., vmax=12, marker='o', picker=5, alpha=0.5)


def plot_beam_time(tab, plot_dir="./"):
    fig, ax = plt.subplots()
    ax.cla()
    colormap = ax.scatter(tab['mjds'], tab['ibeam'], s=(tab['snr'] / tab['snr'].min())*5, c=tab['ibox'], marker='o', alpha=0.5)   
    fig.colorbar(colormap, ax=ax, use_gridspec=True,label='$\\rm Boxcar width\;(index)$')
    ax.set_xlabel('$\\rm Beam Number$', size=12)
    ax.set_ylabel('$\\rm Time (s) $', size=12)
    ax.set_title("giants snr for each beam")
    fig.savefig(plot_dir+"gitnat_beam_time.pdf") 



def plot_giants(tab, plot_dir="./"):
    """
    Plot the un-clustere heimdall output giants.out 
    Parameters
    ----------
    tab : full table from parse_candsfile.
    plot_dir : optional. The default is "./".
    """

    plot_dm_hist(tab, nbins=30, plot_dir=plot_dir, multibeam=False, data_name = None)
    plot_dm_hist(tab, nbins=30, plot_dir=plot_dir, multibeam=True, data_name = None)
    plot_beam_time(tab, plot_dir=plot_dir) 
    
    # subplot or just save plots in each function?  
    fig2, ax2 = plt.subplots(1,2, sharey=True)
    plot_time_dm(ax2[0],tab,duration=None,snr_cut=6.5,snr_thr=25.,mps=3.,multibeam=False,axrange=True,axlabel=True) 
    divider = make_axes_locatable(ax2[1])
    cbar_ax = divider.append_axes("right", size="5%", pad=0.1)
    plot_dm_snr(ax2[1], cbar_ax, tab) 
    fig2.savefig(plot_dir+'giants_dm_time_snr.pdf') 
    fig2.clf()
    plt.close('all') 
    
 
def parse_socket(host, port, selectcols=['itime', 'idm', 'ibox', 'ibeam']):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, port: same with heimdall -coincidencer host:port 
    selectcol: list of str.  Select columns for clustering. 
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.bind((host,port))    # assigns the socket with an address
    s.listen(5)             # accept no. of incoming connections

    while True:
        clientsocket, address = s.accept() # stores the socket details in 2 variables
        print(f"Connection from {address} has been established")
        
        # read in heimdall socket output  
        ascii_letter = clientsocket.recv(1)           # recieves an alphabet whose ASCII value is the size of the message 
        
        if len(ascii_letter) > 0:
            size = ord(ascii_letter.decode('utf-8'))      # ord() returns the ASCII value of a character
            #candsfile = clientsocket.recv(size)           # recieving the actual msg        
            candsfile = clientsocket.recv(int(1e10))  # how to make sure size is big enough?
            
            candsfile = candsfile.decode('utf-8')               # decode the bytes msg 
            print(candsfile)
            
            clientsocket.close()   
            
            # do we want both heimdall.cand and giants.out?  From two sockets?  
            # reading in one gulp at a time. 
            # will modify to continue reading in.  
            print("reading candsfile...")
            tab = ascii.read(candsfile, names=['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam'], guess=False, fast_reader=False, delimiter="\t")
            data = np.lib.recfunctions.structured_to_unstructured(tab[selectcols].as_array())  # ok for single dtype (int)
            snrs = tab['snr'] 
            
            print("table has", len(tab), "rows")
            
            
            return tab, data, snrs
                