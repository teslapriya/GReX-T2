import numpy as np
import time

from T2 import cluster_heimdall,plotting

candsfile = '/home/liam/software/GReX-T2/tests/data/giants.cand'
candsfile = '/home/liam/giants.cand'
tab = cluster_heimdall.parse_candsfile(candsfile)

t0=time.time()
tab3 = cluster_heimdall.cluster_dumb(tab, t_window=0.5)
print(time.time()-t0)

t0=time.time()
cluster_heimdall.cluster_data(tab, metric='euclidean', 
                                  allow_single_cluster=True, 
                                  return_clusterer=False)
    
tab2 = cluster_heimdall.get_peak(tab)
print(time.time()-t0)

#clusterer = cluster_heimdall.cluster_data(tab, return_clusterer=True)
#tab2 = cluster_heimdall.get_peak(tab)

#import matplotlib
#matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
max_point_size=5.
plt.scatter((tab['mjds']-tab['mjds'].min())*86400., tab['dm'],
            s=(max_point_size * tab["snr"] / tab["snr"].mean()),
            alpha=0.5)
plt.scatter((tab2['mjds']-tab2['mjds'].min())*86400., tab2['dm'], s=3.5, )
plt.scatter((tab3['mjds']-tab3['mjds'].min())*86400., tab3['dm'], s=3.5, )
plt.xlabel('Time (s)')
plt.ylabel('DM (pc cm**-3)')
plt.legend(['Heimdall candidates','Clustered (HDBSCAN)','Clustered (Dumb)'])
plt.savefig('here.png')
