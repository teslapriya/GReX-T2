import os

import numpy as np
import glob 

import dsautils.dsa_store as ds
import time
from astropy.time import Time

d = ds.DsaStore()
fmt_out = '%5.9f  %d  %0.2f %0.1f %0.3f %0.2f %s\n'
fnout = '/home/ubuntu/injection_list.txt'

if not os.path.exists(fnout):
    f = open(fnout,'w+')
    f.write('# MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno\n')
    f.close()

params = np.genfromtxt('/home/ubuntu/simulated_frb_params.txt')
flist = glob.glob('/home/ubuntu/data/test_inj*.dat')
nfrb = len(flist)
nfrb=25

for zz in range(1):
    for kk in range(17,21):
        for ii in range(nfrb):
            f = open(fnout,'a')
            subbeam = (2*ii+1) % 64
            beam = 64*(kk-17)+subbeam
            print("Injecting into beam %d"%beam)
#            fn = flist[ii]
            fn = '/home/ubuntu/data/test_inj_0012.dat'
            frbno = fn.split('_')[-1][:4]
            ind = np.where(params[:,-1]==float(frbno))[0]
            DM, SNR, Width_fwhm, spec_ind = params[ind][0][0],params[ind][0][1],params[ind][0][2],params[ind][0][3]
            d.put_dict('/cmd/corr/%d'%kk,{'cmd':'inject','val':'%d-%s-'%(subbeam,fn)})
            imjd = Time.now().mjd
            f.write(fmt_out % (imjd, beam, DM, SNR, Width_fwhm, spec_ind, frbno))
            f.close()
            time.sleep(600)

f.close()