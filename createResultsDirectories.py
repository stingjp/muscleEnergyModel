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
from shutil import rmtree
import distutils


print('Creating all the file directories...\n')

repodir = os.getcwd()
resultbasedir = os.path.join(repodir,'..\\results\\')
# print(repodir)
# print(resultbasedir)

templateanalysis = os.path.join(repodir,'templates\\moco\\analyzeSubject.m')
# print(templateanalysis)
templategrf = os.path.join(repodir,'templates\\moco\\grf_walk.xml')
templategrf_welk = os.path.join(repodir,'templates\\moco\\welk_grf_walk.xml')
geometrydir = os.path.join(repodir,'Geometry\\')
idtemplate = os.path.join(repodir,'templates\\moco\\idguitesting.xml')
template_RRA_dir = os.path.join(repodir, 'RRAfiles\\')
# backpackdir = os.path.join(geometrydir, 'backpack0.vtp')


os.chdir(resultbasedir)
# print(os.getcwd())

# subjects
subjects = ['wals024','wals077','wals088','wals112','wals127','wals128',
            'welk001','welk002','welk003','welk004','welk005','welk006','welk007','welk008','welk009',
            'welk010','welk011','welk012','welk013','welk014',
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
              'welk002welknatural':['trial01','trial02','trial03','trial04'],
              'welk003welknatural':['trial01','trial02','trial03','trial04'],
              'welk004welknatural':['trial01','trial02','trial03','trial04'],
              'welk005welknatural':['trial01','trial02','trial03','trial04'],
              'welk006welknatural':['trial01','trial02','trial03','trial04'],
              'welk007welknatural':['trial01','trial02','trial03','trial04'],
              'welk008welknatural':['trial01','trial02','trial03','trial04'],
              'welk009welknatural':['trial01','trial02','trial03','trial04'],
              'welk010welknatural':['trial01','trial02','trial03','trial04'],
              'welk011welknatural':['trial01','trial02','trial03','trial04'],
              'welk012welknatural':['trial01','trial02','trial03','trial04'],
              'welk013welknatural':['trial01','trial02','trial03','trial04'],
              'welk014welknatural':['trial01','trial02','trial03','trial04'],
              'welk001welkexo':['trial01','trial02','trial03','trial04'],
              'welk002welkexo':['trial01','trial02','trial03','trial04'],
              'welk003welkexo':['trial01','trial02','trial03','trial04'],
              'welk004welkexo':['trial01','trial02','trial03','trial04'],
              'welk005welkexo':['trial01','trial02','trial03','trial04'],
              'welk006welkexo':['trial01','trial02','trial03','trial04'],
              'welk007welkexo':['trial01','trial02','trial03','trial04'],
              'welk008welkexo':['trial01','trial02','trial03','trial04'],
              'welk009welkexo':['trial01','trial02','trial03','trial04'],
              'welk010welkexo':['trial01','trial02','trial03','trial04'],
              'welk011welkexo':['trial01','trial02','trial03','trial04'],
              'welk012welkexo':['trial01','trial02','trial03','trial04'],
              'welk013welkexo':['trial01','trial02','trial03','trial04'],
              'welk014welkexo':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk006welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk012welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk014welknaturalslow':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk006welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk012welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk014welknaturalnatural':['trial01','trial02','trial03','trial04'],
              'welk001welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk002welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk003welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk004welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk005welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk006welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk007welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk008welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk009welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk010welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk011welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk012welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk013welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk014welknaturalexo':['trial01','trial02','trial03','trial04'],
              'welk001welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk002welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk003welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk004welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk005welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk006welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk007welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk008welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk009welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk010welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk011welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk012welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk013welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk014welkexonatural':['trial01','trial02','trial03','trial04'],
              'welk001welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk002welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk003welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk004welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk005welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk006welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk007welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk008welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk009welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk010welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk011welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk012welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk013welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk014welkexoexo':['trial01','trial02','trial03','trial04'],
              'welk001welkexofast':['trial01','trial02','trial03','trial04'],
              'welk002welkexofast':['trial01','trial02','trial03','trial04'],
              'welk003welkexofast':['trial01','trial02','trial03','trial04'],
              'welk004welkexofast':['trial01','trial02','trial03','trial04'],
              'welk005welkexofast':['trial01','trial02','trial03','trial04'],
              'welk006welkexofast':['trial01','trial02','trial03','trial04'],
              'welk007welkexofast':['trial01','trial02','trial03','trial04'],
              'welk008welkexofast':['trial01','trial02','trial03','trial04'],
              'welk009welkexofast':['trial01','trial02','trial03','trial04'],
              'welk010welkexofast':['trial01','trial02','trial03','trial04'],
              'welk011welkexofast':['trial01','trial02','trial03','trial04'],
              'welk012welkexofast':['trial01','trial02','trial03','trial04'],
              'welk013welkexofast':['trial01','trial02','trial03','trial04'],
              'welk014welkexofast':['trial01','trial02','trial03','trial04'],
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
            # do the trial directories
            tempdir_2 = os.path.join(tempdir, cond)
            for keys in welktrials[subj+cond]:
                try:
                    os.mkdir(os.path.join(tempdir_2, keys))
                except:
                    pass
                trialdir = os.path.join(tempdir_2, keys)
                targetfile = os.path.join(trialdir, 'analyzeSubject.m')
                targetfile2 = os.path.join(trialdir, 'grf_walk.xml')
                targetfile3 = os.path.join(trialdir, 'idguitesting.xml')
                targetgeometry = os.path.join(trialdir, 'Geometry')
                targetRRA = os.path.join(trialdir, 'RRAfiles')
                copy(templateanalysis, targetfile)
                copy(templategrf_welk, targetfile2)
                copy(idtemplate, targetfile3)
                try:
                    os.mkdir(os.path.join(trialdir, 'expdata'))
                except:
                    pass
                try:
                    # distutils.dir_util.copy_tree(geometrydir, targetgeometry)
                    copytree(geometrydir, targetgeometry)
                except:
                    rmtree(targetgeometry)
                    copytree(geometrydir, targetgeometry)
                try:
                    # distutils.dir_util.copy_tree(template_RRA_dir, targetRRA)
                    copytree(template_RRA_dir, targetRRA)
                except:
                    rmtree(targetRRA)
                    copytree(template_RRA_dir, targetRRA)

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
                targetfile3 = os.path.join(trialdir, 'idguitesting.xml')
                targetgeometry = os.path.join(trialdir,'Geometry')
                # targetbackpack = os.path.join(targetgeometry,'backpack0.vtp')

                copy(templateanalysis, targetfile)
                copy(templategrf, targetfile2)
                copy(idtemplate, targetfile3)
                
                try:
                    os.mkdir(os.path.join(trialdir, 'expdata'))
                except:
                    pass
                try:
                    # distutils.dir_util.copy_tree(geometrydir, targetgeometry)
                    copytree(geometrydir, targetgeometry)
                except:
                    rmtree(targetgeometry)
                    copytree(geometrydir, targetgeometry)
                try:
                    # distutils.dir_util.copy_tree(template_RRA_dir, targetRRA)
                    copytree(template_RRA_dir, targetRRA)
                except:
                    rmtree(targetRRA)
                    copytree(template_RRA_dir, targetRRA)


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
