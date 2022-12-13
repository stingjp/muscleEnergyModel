'''
Jon Stingel
07/15/2022
Script for running matlab batches from the command line interface. 

'''

import os
from shutil import copyfile
from shutil import copy
from shutil import copytree
import time
# setting up the paths
setup_matlab_path = "G:/Shared drives/Exotendon/musclemodel/muscleEnergyModel/"
# output_fpath = os.path.join(setup_matlab_path, 'matlab_log.txt')

print('Going through each subject...')
## set up all the paths that we will need
repodir = os.path.join('G:\\Shared drives\\Exotendon\\musclemodel\\muscleEnergyModel\\')
resultsbasedir = os.path.join(repodir,'..\\results\\')
# analysisbasedir = os.path.join(repodir,'..\\analysis\\')
# metabolicsbasedir = os.path.join(repodir,'..\\metabolicsResults\\')
# targetmuscleresults = os.path.join(metabolicsbasedir,'muscleInverse\\')
# targetmuscleEMGresults = os.path.join(metabolicsbasedir,'muscleInverseWithEMG\\')



## set up all the subject conditions and trials that we will need
subjects = ['wals024','wals077','wals088','wals112','wals127','wals128',
            'welk001','weld002','welk003','welk004','welk005','welk007','welk008','welk009','welk010','welk011','welk013',
            'jack001','jack002','jack003','jack004','jack005','jack006','jack007','jack008',
            'demb005','demb007','demb009','demb010','demb011','demb012','demb014',
            'sild001','sild001b','sild002','sild003','sild004','sild005','sild007_standing','sild007',
            'sild009','sild010','sild012','sild013','sild001','sild016','sild017','sild018','sild020',
            'sild022','sild024','sild025','sild026','sild027','sild028','sild029','sild030','sild031',
            'sild032','sild033','sild034','sild035']

# conditions
walsconditions = ['walsslack','walslow','walsmed','walshigh','walsmax']
welkconditions = ['welknatural','welkexo','welknaturalslow','welknaturalnatural',
                  'welknaturalexo','welkexonatural','welkexoexo','welkexofast']
jackconditions = ['jackpower1','jackpower2','jackpower3','jackpower4','jackpower5','jackpower6',
                  'jacktau1','jacktau2','jacktau3','jacktau4','jacktau5']
dembconditions = ['dembnoloadfree','dembnoloadslow','dembloadedfree','dembloadedmatched']
sildconditions = ['sildbw0','sildbw5','sildbw10','sild10w0','sild10w5','sild10w10',
                  'sild20w0','sild20w5','sild20w10','sild30w0','sild30w5','sild30w10',
                  'sildbwrun0','sild10wrun0','sild20wrun0','sild30wrun0']

# trials
dembtrials = {'demb005dembnoloadfree':['trial02','trial03','trial06'],
              'demb005dembloadedfree':['trial01','trial03','trial06'],
              'demb007dembnoloadfree':['trial01','trial02','trial03'],
              'demb007dembloadedfree':['trial01','trial02','trial05'],
              'demb009dembnoloadfree':['trial01','trial02','trial04'],
              'demb009dembloadedfree':['trial04','trial05','trial08'],
              'demb010dembnoloadfree':['trial02','trial04','trial05'],
              'demb010dembloadedfree':['trial03','trial04','trial05'],
              'demb011dembnoloadfree':['trial03','trial04','trial05'],
              'demb011dembloadedfree':['trial04','trial05','trial06'],
              'demb012dembnoloadfree':['trial05','trial06','trial07'],
              'demb012dembloadedfree':['trial03','trial04','trial05'],
              'demb014dembnoloadfree':['trial02','trial04','trial05'],
              'demb014dembloadedfree':['trial02','trial03','trial04'],
              
              'demb005dembnoloadslow':[],
              'demb005dembloadedmatched':[],
              'demb007dembnoloadslow':[],
              'demb007dembloadedmatched':[],
              'demb009dembnoloadslow':[],
              'demb009dembloadedmatched':[],
              'demb010dembnoloadslow':[],
              'demb010dembloadedmatched':[],
              'demb011dembnoloadslow':[],
              'demb011dembloadedmatched':[],
              'demb012dembnoloadslow':[],
              'demb012dembloadedmatched':[],
              'demb014dembnoloadslow':[],
              'demb014dembloadedmatched':[]
              }

welktrials = {'welk001welknatural':['trial01','trial02','trial03','trial04'],
              'welk001welkexo':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk001welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk001welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk001welkexofast':['trial01','trial02','trial03','trial04'],
              'welk002welknatural':['trial01','trial02','trial03','trial04'],
              'welk002welkexo':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk002welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk002welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk002welkexofast':['trial01','trial02','trial03','trial04'],
              'welk003welknatural':['trial01','trial02','trial03','trial04'],
              'welk003welkexo':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk003welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk003welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk003welkexofast':['trial01','trial02','trial03','trial04'],
              'welk004welknatural':['trial01','trial02','trial03','trial04'],
              'welk004welkexo':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk004welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk004welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk004welkexofast':['trial01','trial02','trial03','trial04'],

              'welk005welknatural':['trial01','trial02','trial03','trial04'],
              'welk005welkexo':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk005welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk005welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk005welkexofast':['trial01','trial02','trial03','trial04'],

              'welk007welknatural':['trial01','trial02','trial03','trial04'],
              'welk007welkexo':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk007welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk007welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk007welkexofast':['trial01','trial02','trial03','trial04'],

              'welk008welknatural':['trial01','trial02','trial03','trial04'],
              'welk008welkexo':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk008welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk008welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk008welkexofast':['trial01','trial02','trial03','trial04'],

              'welk009welknatural':['trial01','trial02','trial03','trial04'],
              'welk009welkexo':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk009welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk009welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk009welkexofast':['trial01','trial02','trial03','trial04'],

              'welk010welknatural':['trial01','trial02','trial03','trial04'],
              'welk010welkexo':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk010welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk010welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk010welkexofast':['trial01','trial02','trial03','trial04'],

              'welk011welknatural':['trial01','trial02','trial03','trial04'],
              'welk011welkexo':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk011welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk011welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk011welkexofast':['trial01','trial02','trial03','trial04'],

              'welk013welknatural':['trial01','trial02','trial03','trial04'],
              'welk013welkexo':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk013welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk013welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk013welkexofast':['trial01','trial02','trial03','trial04']

              }



## navigate to the results directory
os.chdir(resultsbasedir)


## loop through all the directorie7
# TODO loop through everything and actually gather data!
# going for a datastructure


#### 
# scratch space overwrite for subsets of subj and conditions
# subjects = ['welk002','welk003','welk005','welk007','welk008','welk009','welk010','welk013']
# subjects = ['welk005','welk007','welk008']

subjects = ['welk002'] 
welkconditions = ['welknatural'] 
trials = ['trial03'] # nat 1,2 done

command = "analyzeSubject_python"
# command = "analyzeSubject_setup"


########################################## 
## needs changed from the old script to the new function
for subj in subjects:
    if subj[0:4] == 'wals':
        # make the directories
        tempdir = os.path.join(resultsbasedir,subj)
        # TODO: copy and edit for these subjects
    
    elif subj[0:4] == 'welk':
        # make the directories
        tempdir = os.path.join(resultsbasedir,subj) 
        # loop each condition
        for cond in welkconditions:
            tempdir_2 = os.path.join(tempdir, cond)
            # loop through each trial
            for keys in trials:
                # get the trial dir
                trialdir = os.path.join(tempdir_2, keys)                
                os.chdir(trialdir)
                print(trialdir)
                ## now do stuff
                # print(os.path.isfile('matlab_log.txt'))
                if os.path.isfile('matlab_log.txt'):
                    os.remove('matlab_log.txt')
                    
                status = os.system('matlab -nosplash -nodesktop -logfile matlab_log.txt -wait -r "try, '
                        "run('%s'); disp('SUCCESS'); "
                        'catch ME; disp(getReport(ME)); exit(2), end, exit(0);"\n' 
                        % command
                        )

                # exception catch for failed matlab commands
                if status != 0:
                    # print 'Non-zero exist status. Continuing....'
                    raise Exception('Non-zero exit status.')

                # Wait until output mat file exists to finish the action
                while True:
                    time.sleep(3.0)
                    # output_fpath = os.path.join(trialdir,'\\matlab_log.txt')
                    mat_exists = os.path.isfile('matlab_log.txt')
                    mat_exists2 = os.path.isfile('matlab_log2.txt')

                    if mat_exists or mat_exists2:
                        print('it exists')
                        break
                    print('still in the loop')





    elif subj[0:4] == 'jack':
        # make the directories
        tempdir = os.path.join(resultsbasedir,subj) 
        # TODO: copy and edit for these subjects

    elif subj[0:4] == 'demb':
        # make the directories
        tempdir = os.path.join(resultsbasedir,subj)
        # loop each condition
        for cond in dembconditions:
            tempdir_2 = os.path.join(tempdir, cond)
            # loop through each trial
            for keys in dembtrials[subj+cond]:
                # get the trial dir
                trialdir = os.path.join(tempdir_2, keys)

    elif subj[0:4] == 'sild':
        # make the directories
        tempdir = os.path.join(resultsbasedir,subj) 
        # TODO: copy and edit for these subjects



## TODO:
## copy over the experimental data from a csv file to the same format 
# as the simulations. 


print('\n...end')
print('Ran batch simulations - need to check the: \nmetabolic costs, \nkinematics, \nand any Issues for the subject.')
