# Allow using the osimpipeline git submodule. 
import sys
sys.path.insert(1, 'osimpipeline')      # adds the osimpipeline folder to the directory we are accessing
import os

import yaml                       #!!! can we talk about this again, what is it importing
with open('config.yaml') as f:
    config = yaml.load(f)
if not 'opensim_home' in config:
    raise Exception("You must define the field `opensim_home` in config.yaml "
            "to point to the root of your OpenSim 3.3 installation.")
sys.path.insert(1, os.path.join(config['opensim_home'], 'sdk', 'python'))
sys.path.insert(1, 'perimysium')  #!!! should this be in the config.yaml???

DOIT_CONFIG = {                   #!!! setting up a dict for what? what is verbosity?
        'verbosity': 2,
        'default_tasks': None,
        }

# Settings for plots.
import matplotlib
matplotlib.use('TkAgg')
if matplotlib.__version__[0] == '1':
    raise Exception("Must have matplotlib version 2 to avoid "
            "incorrect bar plots.")
import matplotlib.pyplot as plt
plt.rc('font', family='Helvetica, Arial, sans-serif', size=8)
plt.rc('errorbar', capsize=1.5)
plt.rc('lines', markeredgewidth=1)
plt.rc('legend', fontsize=8)

import osimpipeline as osp
print ">>>>>>>>>>right above for the weird conversion warning<<<<<<<<<<< \n"

# This line is necessary for registering the tasks with python-doit.
from osimpipeline.vital_tasks import *
from osimpipeline.mrs_tasks import *
# import osimpipeline.calibrate_tasks as ct

# Custom tasks for this project.
from tasks import *

# create a dictionary of the muscles in the model - slightly different naming in model vs here
#!!! why are the names different?
muscle_name_map = {
        'addbrev': 'add_brev',
        'addlong': 'add_long',
        'addmagProx': 'add_mag1',
        'addmagMid': 'add_mag2',
        'addmagDist': 'add_mag3',
        'addmagIsch': 'add_mag4',
        'bflh': 'bifemlh',
        'bfsh': 'bifemsh',
        'edl': 'ext_dig',
        'ehl': 'ext_hal',
        'fdl': 'flex_dig',
        'fhl': 'flex_hal',
        'gaslat': 'lat_gas',
        'gasmed': 'med_gas',
        'glmax1': 'glut_max1',
        'glmax2': 'glut_max2',
        'glmax3': 'glut_max3',
        'glmed1': 'glut_med1',
        'glmed2': 'glut_med2',
        'glmed3': 'glut_med3',
        'glmin1': 'glut_min1',
        'glmin2': 'glut_min2',
        'glmin3': 'glut_min3',
        'perbrev': 'per_brev',
        'perlong': 'per_long',
        'piri': 'peri',
        'recfem': 'rect_fem',
        'sart': 'sar',
        'tibant': 'tib_ant',
        'tibpost': 'tib_post',
        'vasint': 'vas_int',
        'vaslat': 'vas_lat',
        'vasmed': 'vas_med'
        }


### why do some of the marker vectors get defined here vs down in the adjustment section??
## placement of markers on the model -> !! what frame or how are these determined?
# upper extremity                                     #!!! where are all of these values coming from?? are these the scaled positions for the marker positions??
RAC = np.array([0.02,0.410,0.140])
LAC = np.array([0.02,0.410,-0.140])
MidAC = (RAC+LAC)/2.0

# lower extremity
z_hip_shift = 0.0                                 #!!! what exactly is this?
RIC = np.array([-0.069,0.09,0.135-z_hip_shift])
LIC = np.array([-0.069,0.09,-0.135+z_hip_shift])
MidIC = (RIC+LIC)/2.0
RWST = np.array([-0.029,0.066,0.14-z_hip_shift])  #!!! why do we add to one
LWST = np.array([-0.029,0.066,-0.14+z_hip_shift]) #!!! why do we subt the other
RASIS = np.array([0.005,0.018,0.135-z_hip_shift])
LASIS = np.array([0.005,0.018,-0.135+z_hip_shift])
MidASIS = (RASIS+LASIS)/2.0
RPSIS = np.array([-0.155,0.035,0.035])
LPSIS = np.array([-0.155,0.035,-0.035])
MidPSIS = (RPSIS+LPSIS)/2.0
MidPELV = (MidPSIS+MidASIS)/2.0
RHJC = np.array([-0.0563,-0.0785,0.07726])
RHJC_P = np.array([-0.0563,MidIC[1],0.07726])
RKJC = np.array([0.001731,-0.002389,-0.008452])
RAJC = np.array([0.0,0.0,0.0])
# interesting that the talus marker (ankle joint) is the origin


# Projected foot markers
RAJC_P = np.array([0.04877,-0.0051348,-0.00792004])
# RHL_P = np.array([-0.0249995,-0.0051348,-0.00500004])
# RHL_P = RHL_P.tolist()
RLT_P = np.array([-0.0188,-0.0031348,0.04892])
RMT_P = np.array([0.0112,-0.0031348,-0.05108,])
RMidT_P = (RLT_P+RMT_P)/2.0
RLT_P = RLT_P.tolist()                  #!!! why do only these ones convert to a list
RMT_P = RMT_P.tolist()
RMidT_P = RMidT_P.tolist()

"""
Library of the markers:

          Anatomical Markers:
    Pelvis:
RAC/LAC:      bony prominence on top of the shoulder. Follow the superior face of clavicle until you find the most superior prominence.
C7:           On large spinous process of 7th cervical vert. 
RIC/LIC:      over the Iliac crest in line with the greater trochanter marker.
RASIS/LASIS:  prominent anterior end of the iliac crest. Follow the anterior part of the iliac crest until you reach the most anterior rounded point.
RPSIS/LPSIS:  prominent posterior end of the iliac crest. Follow the posterior part of the iliac crest until you reach the most posterior rounded point. 
RWST/LWST:    placed over the iliac crest (between ASIS and IC). follow the iliac crest and place thest markers between the ASIS and IC markers.
RTC/LTC:      Put the landmark in the most prominent lateral part of the great trochanter. ask the subject to internally and externally rotate the femur, so you can feel the great trochanter moving. 
    Thigh:
RLK:          locate a tubercle near the center of the lateral condyle of the femur. have the subject flex the knee
RMK:          medial condyle shows a small tubercle. 
    Shank:
RLA:          middle of the lateral malleolus. place index and thumb on the posterior and anterior edge of the lateral malleolus, respectively. Then put the marker in the middle of these fingers.
    Foot:
RMA:          middle of the medial malleolus. place your index and thumb on the posterior and anterior edge of the lateral malleolus, respectively. Put the marker in the middle of these fingers.
RLA:          middle of the lateral malleolus. place your index and thumb on the posterior and anterior edge of the lateral malleolus, respectively. put the marker in the middle of these fingers.
RLH:          placed in the heel counter of the boot. placed on the lateral side of the strap for the Harvard suit. 
RHL:          marker over top of the heel counter. do not place over the strap (Harvard)
RLT:          on the superior and lateral aspect of the fifth metatarsal head
RMT:          on the superior and medial aspect of the first metatarsal head
RDTip:        on the superior aspect of the first metatarsal head
          Tracking Markers:
    Right thigh: -> placed in the lateral/posterior and widest part of the thigh
RSLT:         superior, lateral marker
RSMT:         superior, medial marker
RIMT:         inferior, medial marker
RILT:         inferior, lateral marker
    Right Shank: -> placed just above the boot in the lateral aspect of the shank
RSPS:         superior, posterior
RSAS:         superior, anterior
RIAS:         inferior, anterior
RIPS:         inferior, posterior

        FOR HARVARD:
    Cable markers:
RDAC/LDAC:    placed at the distal end of the bowden cable where the inner cable attaches to the boot.
RPAC/LPAC:    placed oat the proximal end of the bowden cable where the inner cable exits the outer sheath
    Suit Markers:
LMS1,LMS2,LMS3,LMS4
LLS1,LLS2,LLS3,LLS4
LMS1,LLS1
LCW

    ikts.add_ikmarkertask('C7', True, 5.0)
    ikts.add_ikmarkertask_bilateral('AC', True, 10.0)
    ikts.add_ikmarkertask('MidAC', True, 5.0)
    ikts.add_ikmarkertask_bilateral('IC', True, 10.0)
    ikts.add_ikmarkertask_bilateral('WST', True, 5.0)
    ikts.add_ikmarkertask_bilateral('ASIS', True, 20.0)
    ikts.add_ikmarkertask_bilateral('PSIS', True, 20.0)
    ikts.add_ikmarkertask('MidIC', True, 5.0)
    ikts.add_ikmarkertask('MidASIS', True, 5.0)
    ikts.add_ikmarkertask('MidPSIS', True, 5.0)
    ikts.add_ikmarkertask('MidPELV', True, 5.0)
    ikts.add_ikmarkertask('RTC', True, 15.0)
    ikts.add_ikmarkertask('LTC', False, 0.0)
    ikts.add_ikmarkertask('MidTC', False, 0.0)
    ikts.add_ikmarkertask('RSLT', False, 0.0)
    ikts.add_ikmarkertask('RSMT', False, 0.0)
    ikts.add_ikmarkertask('RIMT', False, 0.0)
    ikts.add_ikmarkertask('RILT', False, 0.0)
    ikts.add_ikmarkertask('RLK', True, 15.0)
    ikts.add_ikmarkertask('RMK', True, 15.0)
    ikts.add_ikmarkertask('RSPS', False, 0.0)
    ikts.add_ikmarkertask('RSAS', False, 0.0)
    ikts.add_ikmarkertask('RIPS', False, 0.0)
    ikts.add_ikmarkertask('RIAS', False, 0.0)
    ikts.add_ikmarkertask('RLA', True, 15.0)
    ikts.add_ikmarkertask('RMA', True, 15.0)
    # ikts.add_ikmarkertask('RHL', True, 5.0)
    ikts.add_ikmarkertask('RLH', True, 5.0)
    # ikts.add_ikmarkertask('RMH', True, 5.0)
    ikts.add_ikmarkertask('RLT', True, 5.0)
    ikts.add_ikmarkertask('RMT', True, 5.0)
    ikts.add_ikmarkertask('RDTip', True, 2.0)
    ikts.add_ikmarkertask('RHJC', False, 0.0)
    ikts.add_ikmarkertask('RHJC_P', False, 0.0)
    ikts.add_ikmarkertask('RKJC', True, 15.0)
    ikts.add_ikmarkertask('RAJC', True, 15.0)
    ikts.add_ikmarkertask('RAJC_P', True, 5.0)
    # ikts.add_ikmarkertask('RHL_P', True, 5.0)
    ikts.add_ikmarkertask('RLT_P', True, 5.0)
    ikts.add_ikmarkertask('RMT_P', True, 5.0)
    ikts.add_ikmarkertask('RMidT_P', True, 5.0)
    ikts.add_ikmarkertask('RDAC', False, 0.0)
    ikts.add_ikmarkertask('RPAC', False, 0.0)


"""




#!!! what is actually going on here, is this where the actual scaling is taking place, here is where you are ajusting the values??
marker_map = {       #!!! just an explanation of how this is used would be good
    'delete_marker': ['RCAL'],

    'adjust_marker': [('C7','C7',[-0.055,0.465,0.0017],'torso',True),
                      ('RAC','RAC',RAC,'torso',True),
                      ('LAC','LAC',LAC,'torso',True),
                      ('RASI','RASIS',RASIS,'pelvis',True),
                      ('LASI','LASIS',LASIS,'pelvis',True),
                      ('RPSI','RPSIS',RPSIS,'pelvis',True),
                      ('LPSI','LPSIS',LPSIS,'pelvis',True),
                      ('RHJC','RHJC',RHJC,'pelvis',True),
                      ('RTH1','RSLT',[0.070, -0.150, 0.040],'femur_r',False),
                      ('RTH2','RSMT',[0.070,-0.150,-0.010],'femur_r',False),
                      ('RTH3','RIMT',[0.070,-0.200,-0.010],'femur_r',False),
                      ('RLFC','RLK',[0,-0.404,0.05],'femur_r',True),
                      ('RMFC','RMK',[0,-0.404,-0.05],'femur_r',True),
                      ('RKJC','RKJC',RKJC,'tibia_r',True),
                      ('RTB1','RSPS',[-0.025,-0.195,0.050],'tibia_r',False),
                      ('RTB2','RSAS',[0.025,-0.195,0.050],'tibia_r',False),
                      ('RTB3','RIAS',[0.025,-0.245,0.050],'tibia_r',False),
                      ('RAJC','RAJC',RAJC,'talus_r',True),
                      ('RLMAL','RLA',[-0.005,-0.3888,0.053],'tibia_r',True),
                      ('RMMAL','RMA',[0.006,-0.3888,-0.038],'tibia_r',True),
                      # ('RCAL','RHL',[-0.025,0.01,-0.005],'calcn_r',True),
                      ('RMT5','RLT',[0.160,0.02,0.05],'calcn_r',True),
                      ('RTOE','RMT',[0.190,0.010,-0.050],'calcn_r',True),
                     ],

    'add_marker':    [('MidAC',MidAC,'torso',True),
                      ('RIC',RIC,'pelvis',True),
                      ('LIC',LIC,'pelvis',True),
                      ('MidIC',MidIC,'pelvis',True),
                      ('RWST',RWST,'pelvis',True),
                      ('LWST',LWST,'pelvis',True),
                      ('MidASIS',MidASIS,'pelvis',True),
                      ('MidPSIS',MidPSIS,'pelvis',True),
                      ('MidPELV',MidPELV,'pelvis',True),
                      ('RHJC_P',RHJC_P,'pelvis',True),
                      ('RTC',[-0.010,-0.015,0.085],'femur_r',True),
                      ('LTC',[-0.066276,-0.09349,-0.16226],'pelvis',True),
                      ('MidTC',[-0.066276,-0.09349,0.0],'pelvis',True),
                      ('RILT',[0.070,-0.200,0.040],'femur_r',False),
                      ('RIPS',[-0.025,-0.245,0.050],'tibia_r',False),
                      ('RLH',[0.030,0.020,0.050],'calcn_r',True),
                      # ('RMH',[0.100,0.0225,-0.050],'calcn_r',False),
                      ('RAJC_P',RAJC_P,'calcn_r',True),
                      # ('RHL_P',RHL_P,'calcn_r',True),
                      ('RDTip',[0,0.022,-0.018],'toes_r',True),
                      ('RLT_P',RLT_P,'toes_r',True),
                      ('RMT_P',RMT_P,'toes_r',True),
                      ('RMidT_P',RMidT_P, 'toes_r',True),
                      ('RDAC',[-0.025,0.05,-0.005],'calcn_r',False),
                      ('RPAC',[-0.08377,-0.26695,0.00292],'tibia_r',False),
                     ],

        }

# simple_musculature = True   #!!! are you just choosing this on your own based on the emg that you have or what?
simple_musculature = False
print ('\nGo to dodo.py to figure out what model you need with the scaling weight issue.\n')
if simple_musculature: 
    ## switching to try and see what the weight changes with the torso will do with the metabolics
    # generic_model = 'Rajagopal2015_18musc_muscle_names_probed_iliacus_oneleg.osim'  # 'Rajagopal2015_9musc_right_leg.osim'
    # modified_model = 'Rajagopal2015_18musc_muscle_names_probed_iliacus_oneleg.osim' # 'Rajagopal2015_9musc_right_leg_modified.osim' 
    generic_model = 'Rajagopal2015_18musc_muscle_names_probed_iliacus_noarmweighted.osim'
    modified_model = 'Rajagopal2015_18musc_muscle_names_probed_iliacus_noarmweighted_modified.osim'

    muscle_names = [    #!!! are these in any particular order or nah??
            'glut_max2_r',
            'psoas_r',
            'semimem_r',
            'rect_fem_r',
            'bifemsh_r',
            'vas_int_r',
            'med_gas_r',
            'soleus_r',
            'tib_ant_r'
            ]
else:
    # going to change these models bc of the issue with scaling and the weight of the model. 
    # generic_model = 'Rajagopal2015_right_leg.osim'
    # modified_model = 'Rajagopal2015_right_leg_modified.osim' 
    generic_model = 'Rajagopal2015_noarms.osim'
    modified_model = 'Rajagopal2015_noarms_modified.osim'
    #!!!muscle_names TODO -> did I do this right?? -> I did not
    # muscle_names = [
    #         'add_brev_r',
    #         'add_long_r',
    #         'add_mag3_r',
    #         'add_mag4_r',
    #         'add_mag2_r',
    #         'add_mag1_r',
    #         'bifemlh_r',
    #         'bifemsh_r',
    #         'ext_dig_r',
    #         'ext_hal_r',
    #         'flex_dig_r',
    #         'flex_hal_r',
    #         'lat_gas_r',
    #         'med_gas_r',
    #         'glut_max1_r',
    #         'glut_max2_r',
    #         'glut_max3_r',
    #         'glut_med1_r',
    #         'glut_med2_r',
    #         'glut_med3_r',
    #         'glut_min1_r',
    #         'glut_min2_r',
    #         'glut_min3_r',
    #         'grac_r',
    #         'iliacus_r',
    #         'per_brev_r',
    #         'per_long_r',
    #         'peri_r',
    #         'psoas_r',
    #         'rect_fem_r',
    #         'sar_r',
    #         'semimem_r',
    #         'semiten_r',
    #         'soleus_r',
    #         'tfl_r',
    #         'tib_ant_r',
    #         'tib_post_r',
    #         'vas_int_r',
    #         'vas_lat_r',
    #         'vas_med_r'
    #         ]
    muscle_names = [
            'addbrev_r',
            'addlong_r',
            'addmagDist_r',
            'addmagIsch_r',
            'addmagMid_r',
            'addmagProx_r',
            'bflh_r',
            'bfsh_r',
            'edl_r',
            'ehl_r',
            'fdl_r',
            'fhl_r',
            'gaslat_r',
            'gasmed_r',
            'glmax1_r',
            'glmax2_r',
            'glmax3_r',
            'glmed1_r',
            'glmed2_r',
            'glmed3_r',
            'glmin1_r',
            'glmin2_r',
            'glmin3_r',
            'grac_r',
            'iliacus_r',
            'perbrev_r',
            'perlong_r',
            'piri_r',
            'psoas_r',
            'recfem_r',
            'sart_r',
            'semimem_r',
            'semiten_r',
            'soleus_r',
            'tfl_r',
            'tibant_r',
            'tibpost_r',
            'vasint_r',
            'vaslat_r',
            'vasmed_r'
            ]

study = osp.Study('anklehipexosuit',
        modified_model,
        'Rajagopal2015_reserve_actuators.xml')   #!!! so the reserve actuators just get associated with a degree of freedom?? and then you apply forces witht them?
# from this point on, when we use 'study' it is referring to the large Study(object) titled 'anklehipexosuit'

# muscle names
study.muscle_names = muscle_names

# 'Default' (act/exc squared) or 'Met' (Umberger's metabolic cost model)
study.costFunction = 'Met'
#>>>input
# Shift experimental exoskeleton device torque peaks to line up with ID moments
# 'true' or 'false'                   !!! what is going on here? what does this do??
study.shift_exo_peaks = 'false' # switch to true for modified slack kinematic trials


##need to understand what these markers are used for!!! TODO
# Model markers to compute errors for after IK
error_markers = ['C7','RASIS','LASIS','RPSIS','LPSIS','RIC','LIC','RWST',
                  'LWST','RTC','RSLT','RSMT','RIMT','RILT','RLK','RSPS', # 'RSAS',
                  'RIAS','RIPS','RLA','RLH','RLT']

# error_markers = ['RTC','RLK','RLA','RLH','RLT','RDTip']      # can you really just not add in some of the markers??

study.error_markers = error_markers
subj_nums = ['024','077','088','112','126','127','128']
study.all_subjects = ['subject' + num for num in subj_nums]

##need to go through what the tasks are doing!!! TODO
study.add_task(TaskModifyExperimentalData)
# -> check out the method call to tasks.py, this has functions with these names #

#Setup model: rename muscles, add HBL base marker set
study.add_task(TaskGenericModelSetup, generic_model,
        muscle_name_map=muscle_name_map,
        marker_map=marker_map)


## TODO: figure out if this is the right place for this
# adding in stuff for passive calibration
study.param_dict = dict()
other_param_muscle_list = ['med_gas_r','glut_max2_r','rect_fem_r',
                     'semimem_r','soleus_r','tib_ant_r','vas_int_r']
reduced_param_muscle_list = ['semimem_r','psoas_r']
param_muscle_list = muscle_names # note that the psoas and bifem are not in the one nick sent to me
study.param_dict['optimal_fiber_length'] = param_muscle_list # param_muscle_list
# study.param_dict['tendon_slack_length'] =  param_muscle_list
study.param_dict['muscle_strain'] =        param_muscle_list

study.cost_dict = dict()


'''
So we set up the general structure of the study
then we go to the individual subjects, which each have several methods within them. 
We pass the study to these methods because we are going to modify and add to it for each subject
'''

subjects = config['subjects']         # takes the subject list from the config file
for subj in subjects:
    if subj == 24:
        import subject024                 #!!! what is actually being imported?? - is it the subjectXX.py file stuff?
        subject024.add_to_study(study)    #!!! this goes to the add_to_study(study) object - 
    if subj == 77:
        import subject077
        subject077.add_to_study(study)
    if subj == 88:
        import subject088
        subject088.add_to_study(study)
    if subj == 112:
        import subject112
        subject112.add_to_study(study)
    if subj == 126:
        import subject126
        subject126.add_to_study(study)
    if subj == 127:
        import subject127
        subject127.add_to_study(study)
    if subj == 128:
        import subject128
        subject128.add_to_study(study)

### TODO NEXT -> after figuring out the whole subject import stuff

# Copy modified model and associated actuator files to results directory
study.add_task(TaskCopyGenericModelFilesToResults)  #!!! do these names actually point to anything???

# Copy data files for all study subjects
study.add_task(TaskCopyMotionCaptureData)           #!!! same this what does this name actually carry with it??

cycles = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# cycles = [2, 4, 5, 6, 7, 8, 9]
cycle_nums = ['cycle%02i' % i for i in cycles]              #!!! i don't know what this is??

print "\nhere are the cycle nums"
print cycle_nums

# these are good for getting the met rate for each
study.add_task(TaskAggregateMetabolicRateSlackKinematics,
    cycle_nums=cycle_nums)
study.add_task(TaskAggregateMetabolicRateAssistedKinematics,
    cycle_nums=cycle_nums)

# plots the metabolics trends for subjects and conditions
study.add_task(TaskPlotMetabolicReductions)

# gives the average ankle kin for each subject on one plot
study.add_task(TaskAverageAnkleKinematicsPlot)

# gather the stride length data
study.add_task(TaskAggregatePlotStrideLengthsAssistedKinematics,
  cycle_nums=cycle_nums)
print "\nstill need to do the plot function for stride lengths"



study.add_task(TaskAggregateMomentsExperiment,
  cycle_nums=cycle_nums)
study.add_task(TaskPlotMoments, study.tasks[-1])


study.add_task(TaskAggregateMuscleActivity, cycle_nums=cycle_nums)
study.add_task(TaskPlotMuscleData, study.tasks[-1])
study.add_task(TaskAggregateNormalizedMuscleDynamics, cycle_nums=cycle_nums)
study.add_task(TaskPlotMuscleData, study.tasks[-1])
task_count = 4
for subject in study.subjects:
    study.add_task(TaskAggregateMomentsExperiment, cycle_nums=cycle_nums,
        subject=subject.name)
    study.add_task(TaskPlotMoments, study.tasks[-1], subject=subject.name)
    study.add_task(TaskAggregateMuscleActivity, cycle_nums=cycle_nums,
        subject=subject.name)
    study.add_task(TaskPlotMuscleData, study.tasks[-1],
        subject=subject.name)
    study.add_task(TaskAggregateNormalizedMuscleDynamics, 
        cycle_nums=cycle_nums, subject=subject.name)
    study.add_task(TaskPlotMuscleData, study.tasks[-1],
        subject=subject.name)
    task_count += 6

mods = ['low', 'med', 'high', 'max']
mod_count = 6 
for mod in mods:
    study.add_task(TaskAggregateMomentsMod, mod_name=mod, 
        cycle_nums=cycle_nums)
    study.add_task(TaskPlotMoments, study.tasks[-1], mod_name=mod)
    study.add_task(TaskAggregateMuscleActivity, cycle_nums=cycle_nums,
        mod_name=mod)
    task_count += 3
    study.add_task(TaskPlotMuscleData, study.tasks[-1], mod_name=mod,
        agg_tasks_to_compare=[study.tasks[-task_count]], 
        mod_names_to_compare=['experiment'])
    study.add_task(TaskAggregateNormalizedMuscleDynamics, 
        cycle_nums=cycle_nums, mod_name=mod)
    task_count += 2
    study.add_task(TaskPlotMuscleData, study.tasks[-1], mod_name=mod,
        agg_tasks_to_compare=[study.tasks[-task_count + 2]], 
        mod_names_to_compare=['experiment'])
    task_count += 1
    for isubj, subject in enumerate(study.subjects):
        study.add_task(TaskAggregateMomentsMod, mod_name=mod, 
            cycle_nums=cycle_nums, subject=subject.name)
        study.add_task(TaskPlotMoments, study.tasks[-1], mod_name=mod,
            subject=subject.name)

        agg_task_last = -mod_count - 6*len(study.subjects) - 1
        mod_name_last = 'high'
        agg_task_exp = -task_count + 6*isubj + 3
        mod_name_exp = 'experiment'
            
        study.add_task(TaskAggregateMuscleActivity, cycle_nums=cycle_nums,
            mod_name=mod, subject=subject.name)
        study.add_task(TaskPlotMuscleData, study.tasks[-1], mod_name=mod,
            agg_tasks_to_compare=[study.tasks[agg_task_last],
                                  study.tasks[agg_task_exp]], 
            mod_names_to_compare=[mod_name_last, mod_name_exp], 
            subject=subject.name)
        study.add_task(TaskAggregateNormalizedMuscleDynamics, 
            cycle_nums=cycle_nums, mod_name=mod, subject=subject.name)
        study.add_task(TaskPlotMuscleData, study.tasks[-1], mod_name=mod,
            agg_tasks_to_compare=[study.tasks[agg_task_last],
                                  study.tasks[agg_task_exp]], 
            mod_names_to_compare=[mod_name_last, mod_name_exp], 
            subject=subject.name)
        task_count += 6


