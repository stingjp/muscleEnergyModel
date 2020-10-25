import os
import pandas as pd
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from collections import OrderedDict

from perimysium import postprocessing as pproc


# create a function that will return the X and Y objects for a given muscle and datatype
def importMuscleActivations(muscleStatesFileList, muscles, grffilelist):
    templist = []
    for file in musclestatesfilelist:
        # load the file in to a dataframe with its info
        tempsubj = file[0:7]
        tempcond = file[8:22]
        temptrial = file[23:30]
        tempexp = file[0:4]
        tempfile = os.path.join(musclestatepaths,file)
        tempdata = pd.read_csv(tempfile, skiprows=19, delimiter='\t') # 19 rows for muscle driven 20 for emg
        # print(file)
        # select the data that we want, 
        # normalized along the time axis as % gait cycle
        # somehow select heel strike and string along
        # temptime = tempdata['time']
        muscles = ['bflh','bfsh','gaslat','gasmed','glmax1','glmax2','glmax3','glmed1','glmed2','glmed3',
            'psoas','iliacus','semimem','semiten','soleus','tibant','recfem','vaslat','vasmed'] 
        keys = tempdata.columns
        # print(keys)
        for each in keys:
            if 'activation' not in each and each != 'time':
                tempdata.drop([each], axis=1, inplace=True)

        keys2 = tempdata.columns
        for each in keys2:
            if not any(string in each for string in muscles) and each != 'time':
                tempdata.drop([each], axis=1, inplace=True)

        keys3 = tempdata.columns
        # print(len(keys3))
        # print(keys3)

        muscleactdict = {'subjectname':[tempsubj], 'condname':[tempcond],
                        'trialname':[temptrial], 'experimentname':[tempexp]}
        

        # working
        grffile = file[0:30] + 'grf.mot'
        # print(grffile)
        temptime = tempdata['time'].values
        starttime = temptime[0]
        endtime = temptime[-1]
        rightstrike, leftstrike, rightoff, leftoff = pproc.gait_landmarks_from_grf(os.path.join(grfpaths, grffile),
                                        right_grfy_column_name='ground_force_r_vy',
                                        left_grfy_column_name='ground_force_l_vy',
                                        threshold=1e-5,
                                        do_plot=False,
                                        min_time=starttime,
                                        max_time=endtime,
                                        plot_width=6,
                                        show_legend=True)

        # print('start time %f' % starttime)
        # print('end time %f' % endtime)
        # print(rightstrike)
        # print(rightoff)
        throw, temptime = scipy.signal.resample(x=temptime, t=temptime, num=100)
        # print(temptime)
        # print(len(rightoff))
        # print(len(rightstrike))
        # get new index for movements
        if len(rightoff) == 1:
            # we have a right toeoff
            # get index at that time
            offidx = np.argmax(temptime > rightoff[0])
            # print('offidx: %f' % offidx)
            useoff = True
        else:
            # get the heel strike index
            strikeidx = np.argmax(temptime > rightstrike[0])
            # print('strikeidx: %f' % strikeidx)
            useoff = False

        for each in keys3:
            # print(each)
            temptimethrow = temptime
            temp = tempdata[each].values
            temp = scipy.signal.resample(x=temp, num=100)
            
            # handle actual vectors now
            if each == 'time':
                temptime2 = temptime
                length = temptime[-1] - temptime[0]
                start = temptime[0]
                end = temptime[-1]
                # change to % gait cycle
                for idx, val in enumerate(temptime2):
                    temptime2[idx] = (val-start)/length*100
                
                muscleactdict[each] = [temptime2]
                # print(temptime2)
            else:
                # now to shift the values depending on the indexes found
                if useoff:
                    # use the toeoff index
                    # print(offidx)
                    if offidx == 60: 
                        newtemp = temp
                    elif offidx < 60:
                        numfrom60 = 60 - offidx
                        rearbump = temp[-numfrom60:]
                        # print(rearbump.shape)
                        interm = temp[0:-numfrom60]
                        # print(interm.shape)
                        rearbump = np.append(rearbump, interm)
                        newtemp = rearbump
                        # print(rearbump.shape)
                    elif offidx > 60:
                        numfrom60 = offidx - 60
                        # take this number off the front
                        frontbump = temp[0:numfrom60]
                        # print(frontbump.shape)
                        interm = temp[numfrom60:]
                        # print(interm.shape)
                        interm = np.append(interm, frontbump)
                        # print(interm.shape)
                        newtemp = interm
                else:
                    if strikeidx == 0:
                        pass
                    else:
                        # use the heel strike
                        numfrom0 = strikeidx 
                        frontbump = temp[0:numfrom0]
                        interm = temp[numfrom0:]
                        interm = np.append(interm, frontbump)
                        newtemp = interm

                # plt.figure()
                # plt.plot(temptime, temp)
                # plt.plot(temptime, newtemp)
                # plt.show()
                muscleactdict[each] = [newtemp]


        # print(muscleactdict)
        # import time
        # time.sleep(10)
        tempmuscle_df = pd.DataFrame(muscleactdict)
        templist.append(tempmuscle_df)


    # time.sleep(30)
    # # combine each of the dataframes into one
    act_df = pd.concat(templist, ignore_index=True)
    act_df.sort_values(by=['subjectname', 'condname','trialname'], inplace=True)
    act_df.drop([40], inplace=True)
    print('\nMake sure that you actually want to drop this one entry from subject 14!\n')
    act_df.reset_index(inplace=True, drop=True)

    return act_df

################################################################

def no_top_right(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')



def plot_activations(experimentnames, conditionnames):
    fontsize=8
    n_rows = 11
    n_cols = len(conditionnames)
    # for emg_name, musc_names in emg_sensor_map.items():
        # print(emg_name)
    
    def plot_cond(i_col, condition):
        i_row = 0
        for emg_name, musc_names in emg_sensor_map.items():
            # setup subplot
            ax = plt.subplot2grid((n_rows, n_cols), (i_row, i_col))

            if n_cols == 2 and i_row == 0:
                ax.set_title(loadcond_nice_name[condition],
                            fontsize=fontsize, fontweight='bold')
            ax.tick_params(direction='out', axis='both')
            no_top_right(ax)

            ax.set_ylim((0,1))
            ax.set_xlim((0,100))

            if i_row + 1 == n_rows:
                ax.set_xlabel('time (% gait cycle)', fontsize=fontsize)
                ax.set_xticks(np.array([0,20,40,60,80,100]))
            else:
                ax.axes.get_xaxis().set_ticklabels([])
                ax.set_xticks(np.array([0,20,40,60,80,100]))
            ax.set_yticks([0,1])
            if i_col == 0:
                pass
            else:
                ax.axes.get_yaxis().set_ticklabels([])


            if i_col == 0:
                if n_cols == 1:
                    xpos = -25
                else:
                    xpos = -28
                plt.text(xpos, 0.5, sensor_long_names[emg_name], ha='center',
                        va='center', fontsize=fontsize)

            # here is where I could figure out how to get the EMG values

            # plot simulation values
            handles = []
            labels = []
            for im, musc_prefix in enumerate(musc_names):
                # print(im)
                # print(musc_prefix)
                # need to get the muscle mean curve and the std
                # print(condition)
                # print(sim_df_map[musc_prefix])
                if condition == 'dembnoloadfree':
                    musc_act_mean = noload_act_df[sim_df_map[musc_prefix]].values.mean(axis=0)
                    musc_act_std = noload_act_df[sim_df_map[musc_prefix]].values.mean(axis=0)
                    timecurve = noload_act_df['time'].values.mean(axis=0)
                    # print(musc_act_mean.shape)
                    # print(musc_act_std.shape)
                    # import time
                    # time.sleep(10)
                elif condition == 'dembloadedfree':
                    musc_act_mean = loaded_act_df[sim_df_map[musc_prefix]].values.mean(axis=0)
                    musc_act_std = loaded_act_df[sim_df_map[musc_prefix]].values.mean(axis=0)
                    timecurve = loaded_act_df['time'].values.mean(axis=0)
                if len(musc_names) == 1:
                    label = 'simulation'
                else:
                    label = musc_prefix
                    if musc_prefix in sim_legend_map:
                        label = sim_legend_map[musc_prefix]

                # plot the curve

                ax.plot(timecurve, musc_act_mean, lw=3, label=label, color=colors[im])
                ax.fill_between(timecurve, musc_act_mean+musc_act_std, 
                            musc_act_mean-musc_act_std, facecolor=colors[im],
                            alpha=0.5)
            if i_col == (n_cols - 1):
                first_legend = plt.legend(frameon=False, fontsize=fontsize,
                                    bbox_to_anchor=(1, 1.1), loc='upper left')
                plt.gca().add_artist(first_legend)
                # if i_row == 0:
                #     pl.legend((emgplot,), ('EMG',), loc='upper right',
                #             bbox_to_anchor=(1, 1.1),
                #             frameon=False, fontsize=fontsize)



            i_row += 1
            
    ## TODO: setup a loop for the experiments to create these figures for each experiment


    # working with only demb
    fig = plt.figure(figsize=(5.2, 8.5))
    plot_cond(0, conditionnames[0])
    plot_cond(1, conditionnames[1])
    fig.subplots_adjust(left=0.16, right=0.81, top=0.97, bottom=0.068)
    fig.tight_layout()
    fig.savefig('demb_muscleactivity_withstd.png', dpi=600)

    # plt.show()

    return



if __name__ == '__main__':
    ## import all the muscle activations 
    # set paths
    repobasedir = os.getcwd()
    musclestatepaths = os.path.join(repobasedir,'..\\muscleDrivenResults\\')
    musclestatesfilelist = os.listdir(musclestatepaths)
    grfpaths = os.path.join(repobasedir,'..\\expfiles\\grffiles\\')
    grffilelist = os.listdir(grfpaths)

    ####
    blue = (0.2980392156862745, 0.4470588235294118, 0.6901960784313725)
    green = (0.3333333333333333, 0.6588235294117647, 0.40784313725490196)
    red = (0.7686274509803922, 0.3058823529411765, 0.3215686274509804)
    colors = [blue, green, red]



    sensor_names = ['SOL','GAS','TA','MH','BF','VL','VM','RF','GMAX','GMED','HIP']
    sensor_long_names = {
            'SOL': 'soleus',
            'GAS': 'gastrocnemius',
            'TA': 'tibialis\nanterior',
            'MH': 'medial\nhamstrings',
            'BF': 'biceps\nfemoris',
            'VL': 'vastus\nlateralis',
            'VM': 'vastus\nmedialis',
            'RF': 'rectus\nfemoris',
            'GMAX': 'gluteus\nmaximus',
            'GMED': 'gluteus\nmedius',
            'HIP': 'iliopsoas'}

    emg_sensor_map = OrderedDict()
    emg_sensor_map['GMAX'] = ['glmax3', 'glmax2', 'glmax1']
    emg_sensor_map['GMED'] = ['glmed3', 'glmed2', 'glmed1']
    emg_sensor_map['MH'] = ['semimem', 'semiten']
    emg_sensor_map['BF'] = ['bfsh', 'bflh']
    emg_sensor_map['RF'] = ['recfem']
    emg_sensor_map['VL'] = ['vaslat']
    emg_sensor_map['VM'] = ['vasmed']
    emg_sensor_map['GAS'] = ['gaslat', 'gasmed']
    emg_sensor_map['SOL'] = ['soleus']
    emg_sensor_map['TA'] = ['tibant']
    emg_sensor_map['HIP'] = ['psoas', 'iliacus']

    all_loadconds = [
                'noload/slow',
                'noload/free',
                'loaded/free',
                'loaded/matched'
                ]

    sim_legend_map = {
            'bfsh': 'short head',
            'bflh': 'long head',
            'gaslat': 'lateral',
            'gasmed': 'medial',
            'glmax3': 'posterior',
            'glmax2': 'intermed.',
            'glmax1': 'anterior',
            'glmed3': 'posterior',
            'glmed2': 'intermed.',
            'glmed1': 'anterior',
            'semimem': 'semimem.',
            'semiten': 'semiten',
            'psoas': 'psoas',
            'iliacus': 'iliacus'
            }

    sim_df_map = {
            'bfsh': '/forceset/bfsh_r/activation',
            'bflh': '/forceset/bflh_r/activation',
            'gaslat': '/forceset/gaslat_r/activation',
            'gasmed': '/forceset/gasmed_r/activation',
            'glmax3': '/forceset/glmax3_r/activation',
            'glmax2': '/forceset/glmax2_r/activation',
            'glmax1': '/forceset/glmax1_r/activation',
            'glmed3': '/forceset/glmed3_r/activation',
            'glmed2': '/forceset/glmed2_r/activation',
            'glmed1': '/forceset/glmed1_r/activation',
            'semimem': '/forceset/semimem_r/activation',
            'semiten': '/forceset/semiten_r/activation',
            'vaslat': '/forceset/vaslat_r/activation',
            'vasmed': '/forceset/vasmed_r/activation',
            'tibant': '/forceset/tibant_r/activation',
            'recfem': '/forceset/recfem_r/activation',
            'soleus': '/forceset/soleus_r/activation',
            'iliacus': '/forceset/iliacus_r/activation',
            'psoas': '/forceset/psoas_r/activation'
    }

    loadcond_nice_name = {'dembnoloadfree': 'no load',
                          'dembloadedfree': 'loaded (no assistance)'}
    

    # TODO figure out how to use that structure to my advantage
    muscles = ['bflh','bfsh','gaslat','gasmed','glmax1','glmax2','glmax3','glmed1','glmed2','glmed3',
        'psoas','iliacus','semimem','semiten','soleus','tibant','recfem','vaslat','vasmed'] 

    conditions = ['dembnoloadfree', 'dembloadedfree'] 
    experiments = ['demb']
    # get the activations dataframe
    act_df = importMuscleActivations(musclestatesfilelist, muscles, grffilelist)    
    noload_act_df = act_df[act_df['condname'] == 'dembnoloadfree']
    loaded_act_df = act_df[act_df['condname'] == 'dembloadedfree']

    # print(noload_act_df.columns)

    # print(noload_act_df)

    # timecurve_subj10 = noload_act_df['time'].loc[21]
    # solcurve_subj10 = noload_act_df['/forceset/soleus_r/activation'].loc[21]
    # # print(timecurve.shape)
    # # print(solcurve.shape)
    # timecurve_subj11 = noload_act_df['time'].loc[27]
    # solcurve_subj11 = noload_act_df['/forceset/soleus_r/activation'].loc[27]

    


    # meancurve = act_df['/forceset/psoas_r/activation'].values.mean(axis=0)
    # stdcurve = act_df['/forceset/psoas_r/activation'].values.std(axis=0)

    # print(meancurve.shape)
    # print(stdcurve.shape)

    # test = act_df['/forceset/psoas_r/activation'].loc['']
    # print(test.shape)
    # print(test)

    # import time
    # time.sleep(10)
    # call the plot function
    plot_activations(experiments, conditions)