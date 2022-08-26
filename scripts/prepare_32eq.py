# I want to analyze the data in Julia, but the models are saved in Tensorflow format XXXX
# This script is called from julia so that it generated virtual traces, which can be further analyzed
import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import numpy as np
import tensorflow as tf
import pandas as pd
import obspy
import h5py



event_list=pd.read_csv(r"/data2/isha/events_list.csv")
ev_codes=np.array(event_list.iloc[:,2])
eq_names=ev_codes.tolist()
mname="256nzi_544nz0_600e"
model=tf.keras.models.load_model("/data/losses_plot/eq32_p_symae_"+mname)

strm=obspy.read("/data2/isha/unique_traces_32/*")

Wtraces=[]
eqlist=[]
Gbins=[]
nevlist=[]
for eq_name in eq_names:
    st1=strm.select(id=eq_name+'*')
    nev_list=np.array([tr.stats.sac.nevid for tr in st1])
    nev_list=np.unique(nev_list)
    # print(nev_list, st1)
    
    for nev in nev_list:
        st2=obspy.Stream()
        for tr in st1:
            if tr.stats.sac.nevid==nev:
                st2.append(tr)
                
        d=np.zeros((len(st2), 801))
        #print(np.shape(d))
        
        for k,tr in enumerate(st2):
            # reading data per bin
            d[k]=tr.data
            # normalizing
            d[k] = d[k] / np.std(d[k])
        
        # append
        # print(np.shape([d]))
        d=np.array([d])
        Wtraces.append(model.nencoder(d))
        Gbins.append(model.symencoder(d))
        eqlist.append(eq_name)
        nevlist.append(nev)


Wall=np.concatenate(Wtraces, axis=1)
nbins=len(Gbins)

def generate_virt(eq_name):
    eq_bins=[]
    for i in range(nbins):
        if(eqlist[i]==eq_name):
            eq_bins.append(i)

    virt=[]
    for ibin in eq_bins:
        G=Gbins[ibin]
        Grepeat = [np.repeat(G, np.shape(Wall)[1], axis=0)]
        Z = [np.concatenate((Grepeat, Wall), axis=2)]
        virt.append(np.array(model.decoder(Z))[0])  # virtual trace
    return virt




# # write encoded arrays to disk
# f = h5py.File("encoded"+mname+".hdf5", 'w')
# for i in range(len(Gbins)):
#     if(str(i) not in f.keys()):
#         g = f.create_group(str(i))
#     else:
#         g = f[str(i)]
    
#     g["Wtraces"] = Wtraces[i]
#     g["Gbins"] = Gbins[i]
#     g["eq"] = eqlist[i]
#     g["nev"] = nevlist[i]

    
# f.attrs["eqlist"] = eqlist
# f.attrs["nevlist"] = nevlist                             
# f.close()
