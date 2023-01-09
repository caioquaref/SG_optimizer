import os
import mat73
import pandas as pd
import numpy as np

"""
Iteration number counter.

if niter := None:  # need to figure this out! how to handle glabal and local variables between executions of this file
    niter = 0
else:
    niter += 1

<---------------------------------------------------------------------------------->
"""
dvpaw = 3.84  # Delta Voltage Pedal at WOT
Hz = 10  # Frequency


gofast_path = 'C:\PWT-Tools\GOFAST\GOFAST_2021_07'
model_path = 'C:\PWT-Tools\GOFAST\Data\GG_Perf_model'
lnc_file = 'Dodge_Hornet_GMET4_256HP_V3.1_NormalLaunch.lnc'

"""
To avoid a GF license error, the call from Python IDE must be made from the 
same folder as GOFAST.EXE, so change the directory temporarily and revert
upon completion.
"""
current_file = os.getcwd()
os.chdir(gofast_path)
full_path = gofast_path + '\GOFAST.exe' + " " + '"' + model_path + '\\' + lnc_file + '"'
#os.system(full_path)
os.chdir(current_file)

#getting data from _output folder
out_dir = lnc_file.replace('.lnc','_output')
out_file = model_path + '\\' + out_dir + '\\' + 'TH_Veh1_Task4_Iter1.mat'

annots = mat73.loadmat(out_file, only_include='tTH')
values = list(zip(annots['tTH']['time']['v'], annots['tTH']['dvdt_veh']['v'], annots['tTH']['v_veh_MPH']['v'], annots['tTH']['s_veh']['v'], annots['tTH']['pos_acc']['v']))  # , annots['tTH']['FEcum_comb']['v']
columns = ['time', 'accel','v_veh_MPH', 'dist', 'acc pedal']  #,'FEcum_comb'
df = pd.DataFrame(values, columns=columns)

df['DELPVS_sim'] = df['acc pedal'] * dvpaw / 100
#df['unadj_MPG'] = df['FEcum_comb'] * 0.6213712 * 3.785412
sec = np.arange(0, df.iat[df.index[-1],0], 1/Hz)
mph = np.interp(sec, df['time'], df['v_veh_MPH'])
delpvs = np.interp(sec, df['time'], df['DELPVS_sim'])

print(mph)