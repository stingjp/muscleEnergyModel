'''
Jon stingel
08/12/2020
File to create all the results file structure, as they are 
used in the matlab scripts - must be present
'''

import os
from shutil import copyfile
from shutil import copy
from shutil import copytree


print('Creating all the file directories...\n')

repodir = os.getcwd()
resultbasedir = os.path.join(repodir,'..\\results\\')
# print(repodir)
# print(resultbasedir)

templateanalysis = os.path.join(repodir,'templates\\moco\\analyzeSubject.m')
# print(templateanalysis)
templategrf = os.path.join(repodir,'templates\\moco\\grf_walk.xml')
geometrydir = os.path.join(repodir,'Geometry')

os.chdir(resultbasedir)
# print(os.getcwd())

# subjects
subjects = ['wals024','wals077','wals088','wals112','wals127','wals128',
            'welk001','weld002','welk003','welk004',
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

## TODO
# have to figure out the trials for each of the conditions, 
# and then get the expdata folder and copy script



for subj in subjects:
    if subj[0:4] == 'wals':
        # make the directories
        tempdir = os.path.join(resultbasedir,subj) 
        try:
            os.mkdir(tempdir)
        except:
            # print('\nDirectory exists.')
            pass
        # loop the conditions
        for cond in walsconditions:
            try:
                os.mkdir(os.path.join(tempdir, cond))
            except:
                # print('\nDirectory exists.')
                pass
    elif subj[0:4] == 'welk':
        # make the directories
        tempdir = os.path.join(resultbasedir,subj) 
        try:
            os.mkdir(tempdir)
        except:
            # print('\nDirectory exists.')
            pass
        # loop the conditions
        for cond in welkconditions:
            try:
                os.mkdir(os.path.join(tempdir, cond))
            except:
                # print('\nDirectory exists.')
                pass
    elif subj[0:4] == 'jack':
        # make the directories
        tempdir = os.path.join(resultbasedir,subj) 
        try:
            os.mkdir(tempdir)
        except:
            # print('\nDirectory exists.')
            pass
        # loop the conditions
        for cond in jackconditions:
            try:
                os.mkdir(os.path.join(tempdir, cond))
            except:
                # print('\nDirectory exists.')
                pass
    elif subj[0:4] == 'demb':
        # make the directories
        tempdir = os.path.join(resultbasedir,subj) 
        try:
            os.mkdir(tempdir)
        except:
            # print('\nDirectory exists.')
            pass
        # loop the conditions
        for cond in dembconditions:
            try:
                os.mkdir(os.path.join(tempdir, cond))
            except:
                # print('\nDirectory exists.')
                pass
            # do the trial directories
            tempdir_2 = os.path.join(tempdir, cond)
            for keys in dembtrials[subj+cond]:
                try:
                    os.mkdir(os.path.join(tempdir_2, keys))
                except:
                    # print('\nTrial directory exists.')
                    pass
                trialdir = os.path.join(tempdir_2, keys)
                targetfile = os.path.join(trialdir, 'analyzeSubject.m')
                targetfile2 = os.path.join(trialdir, 'grf_walk.xml')
                targetgeometry = os.path.join(trialdir,'Geometry')

                copy(templateanalysis, targetfile)
                copy(templategrf, targetfile2)
                
                # if not os.path.exists(targetfile):
                #     copyfile(templateanalysis, targetfile)
                # if not os.path.exists(targetfile2):
                #     copyfile(templategrf, targetfile2)
                try:
                    os.mkdir(os.path.join(trialdir, 'expdata'))
                    copytree(geometrydir, targetgeometry)
                except:
                    pass
    elif subj[0:4] == 'sild':
        # make the directories
        tempdir = os.path.join(resultbasedir,subj) 
        try:
            os.mkdir(tempdir)
        except:
            # print('\nDirectory exists.')
            pass
        # loop the conditions
        for cond in sildconditions:
            try:
                os.mkdir(os.path.join(tempdir, cond))
            except:
                # print('\nDirectory exists.')
                pass

# end 






print('\n...end')
print('\nNow copy the experimental data that is necessary for each of the subjects, conditions, and trials into the respective "expdata" folder in the results directories we just made.')