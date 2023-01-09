import pandas as pd
"""
Files that have tables:
mot, alt, gsm, pqm, vsc, con, gbm.

The user should specify which tables columns he wants to change and how.
"""
model = 'C:\PWT-Tools\GOFAST\Data\GG_Perf_model'
file = 'Chrysler Adv I4 GDI TC BZ 1995DI TcFG 4V 4L 263.6062CV 405Nm E5 gdi - 25_33_ZeroMM_256HP_nodfso_perf (TEST).mot' # input from main program/user.

#columns = ['Variable', 'Value', 'Unit']  # standard Gofast file format
df_param = pd.read_csv(model + '\\' + file, sep='\t', on_bad_lines='skip', skiprows=1, header=None) # names=columns,skiprows=81
last_index = df_param.index[-1] + 1
df_param.drop(index=df_param.index[-1], axis=0, inplace=True)
print(df_param)
df_table = pd.read_csv(model + '\\' + file, sep='\t', skiprows=last_index, header=[0,1]) # names=columns,skiprows=81
print(df_table)

if file[-3:] in ['mot', 'pqm', 'alt', 'vsc', 'gsm', 'con', 'gbm']:
    print('yes')
#for i in range(len(variables)):
#    index = df.loc[df['Variable'] == variables[i]].index[0]
#    df.iat[index,1] = new_values[i]