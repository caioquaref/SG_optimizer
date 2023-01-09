from numpy import NaN
import pandas as pd

def ModifyFile(model_path, file_to_mod, variables, new_values):
    """
    Modifies a Gofast file according to user's desired varibles and main program's given values.
    Variables to mod in a given optimization should be specified by the user. Values need to come
    from the main program, changing with each iteration.

    Arguments:
    - model_path  = string specifying path to model;
    - file_to_mod = string specifying the file name to change values (vhc, gsm, mot, etc);
    - variables   = list with variables to mod (depends on file extention);
    - new_values  = list of new value inputs for modded file;

    Currently this function does not support modification of table values. Therefore, 'mot', 'pqm',
    'alt', 'vsc', 'gsm', 'con' and 'gbm' files are not usable.
    
    """

    if file_to_mod[-3:] in ['mot', 'pqm', 'alt', 'vsc', 'gsm', 'con', 'gbm']:
        exit("Cannot modify 'mot', 'pqm', 'alt', 'vsc', 'gsm', 'con', or 'gbm' file type in this program version!")  #figure if out

    else:
        columns = ['Variable', 'Value', 'Unit']  # standard Gofast file format
        df = pd.read_csv(model_path + '\\' + file_to_mod, sep='\t', skiprows=1, names=columns, header=None)

        for i in range(len(variables)):
            index = df.loc[df['Variable'] == variables[i]].index[0]
            df.iat[index,1] = new_values[i]

        df.loc[-1] = ['<Info>', NaN, NaN]  # adding a row with standard <Info> text 
        df.index = df.index + 1  # shifting index
        df = df.sort_index()  # sorting by index

        df.to_csv(model_path + '\\' + file_to_mod, sep='\t', header=False, index=False)

# Random inputs to test function

#model = 'C:\PWT-Tools\GOFAST\Data\GG_Perf_model'
#file = 'ZF 948TE AT9 5.54 FDR3.734 AWD - 7c_Trans (TEST) - Copia.trs' # input from main program/user.
#var = ['trs_zmp', 'trs_zm1', 'trs_zm2']
#values = [44, 88, 132]
#ModifyFile(model, file, var, values)