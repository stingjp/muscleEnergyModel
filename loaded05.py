import os

import osimpipeline as osp
import tasks
import h5py
import helpers

def scale_setup_fcn(util, mset, sset, ikts):

    # Torso
    #m = util.Measurement('torso_Z', mset)
    #m.add_markerpair('RAC', 'LAC')
    #m.add_bodyscale('torso','Z')

    #m = util.Measurement('torso_Y', mset)
    #m.add_markerpair('MidIC','MidAC')
    #m.add_bodyscale('torso','Y')

    # Manual scale for torso
    # m = util.Scale('torso', 1.0, 0.964239, 1.115168, sset)  
    m = util.Measurement('torso_X', mset)
    m.add_markerpair('R.Clavicle','R.Shoulder')
    m.add_markerpair('L.Clavicle','L.Shoulder')
    m.add_markerpair('MidClavicle','S2')
    m.add_bodyscale('torso','X')

    m = util.Measurement('torso_Y', mset)
    m.add_markerpair('MidASIS','MidClavicle')
    m.add_markerpair('MidPSIS','MidClavicle')
    m.add_markerpair('MidASIS','S2')
    m.add_markerpair('MidPSIS','S2')
    m.add_bodyscale('torso','Y')

    m = util.Measurement('torso_Z', mset)
    m.add_markerpair('R.Shoulder','L.Shoulder')
    m.add_markerpair('R.Clavicle','L.Clavicle')
    m.add_bodyscale('torso','Z')

    # Pelvis
    m = util.Measurement('pelvis_X', mset)
    m.add_markerpair('MidASIS','MidPSIS')
    m.add_bodyscale('pelvis', 'X')

    m = util.Measurement('pelvis_Y', mset)
    m.add_markerpair('L.ASIS','L_HJC')
    m.add_markerpair('R.ASIS','R_HJC')
    m.add_markerpair('L.PSIS','L_HJC')
    m.add_markerpair('R.PSIS','R_HJC')
    m.add_bodyscale('pelvis', 'Y')

    m = util.Measurement('pelvis_Z', mset) # RWST and LWST
    m.add_markerpair('R.ASIS','L.ASIS')
    m.add_markerpair('R_HJC','L_HJC')
    m.add_bodyscale('pelvis', 'Z')

    # Right Femur
    m = util.Measurement('femur_r_Y', mset)
    m.add_markerpair('L.Knee','L_HJC')
    m.add_markerpair('R.Knee','R_HJC')
    m.add_bodyscale('femur_r','Y')
    m.add_bodyscale('femur_r','X')
    m.add_bodyscale('femur_r','Z')
    m.add_bodyscale('femur_l','Y')
    m.add_bodyscale('femur_l','X')
    m.add_bodyscale('femur_l','Z')


    # m = util.Measurement('femur_r_Z', mset)
    # m.add_markerpair('RLK','RMK')
    # m.add_bodyscale('femur_r','Z')
    # m.add_bodyscale('femur_r','X')

    # Right Tibia
    m = util.Measurement('tibia_r_Y', mset)
    m.add_markerpair('L.Knee','L.Ankle')
    m.add_markerpair('R.Knee','R.Ankle')
    m.add_bodyscale('tibia_r','Y')
    m.add_bodyscale('tibia_r','X')
    m.add_bodyscale('tibia_r','Z')
    m.add_bodyscale('tibia_l','Y')
    m.add_bodyscale('tibia_l','X')
    m.add_bodyscale('tibia_l','Z')

    # m = util.Measurement('tibia_r_Z', mset)
    # m.add_markerpair('RMA','RLA')
    # m.add_bodyscale('tibia_r','Z')
    # m.add_bodyscale('tibia_r','X')

    # Right Foot
    # '''
    m = util.Measurement('foot_r_X', mset)
    m.add_markerpair('L.HeelGround','L.MT5Ground')
    m.add_markerpair('R.HeelGround','R.MT5Ground')
    m.add_markerpair('L.HeelGround','L.ToeGround')
    m.add_markerpair('R.HeelGround','R.ToeGround')
    m.add_bodyscale('talus_r','X')
    m.add_bodyscale('calcn_r','X')
    m.add_bodyscale('toes_r','X')
    m.add_bodyscale('talus_l','X')
    m.add_bodyscale('calcn_l','X')
    m.add_bodyscale('toes_l','X')

    m = util.Measurement('foot_r_Y', mset)
    m.add_markerpair('L.AnkleJoint','L.AnkleJointGround') #RLH and RLA
    m.add_markerpair('R.AnkleJoint','R.AnkleJointGround')
    m.add_bodyscale('talus_r','Y')
    m.add_bodyscale('calcn_r','Y')
    m.add_bodyscale('toes_r','Y')
    m.add_bodyscale('talus_l','Y')
    m.add_bodyscale('calcn_l','Y')
    m.add_bodyscale('toes_l','Y')

    m = util.Measurement('foot_r_Z', mset)
    m.add_markerpair('L.MT5Ground','L.ToeGround')
    m.add_markerpair('R.MT5Ground','R.ToeGround')
    m.add_bodyscale('talus_r','Z')
    m.add_bodyscale('calcn_r','Z')
    m.add_bodyscale('toes_r','Z')
    m.add_bodyscale('talus_l','Z')
    m.add_bodyscale('calcn_l','Z')
    m.add_bodyscale('toes_l','Z')
    # '''
        # m = util.Scale('talus_r',   0.953832,   1.5,    1.13861,    sset)
        # m = util.Scale('calcn_r',   0.953832,   1.5,    1.13861,    sset)
        # m = util.Scale('toes_r',    0.953832,   1.5,    1.13861,    sset)

        # # for the case of scaling both legs and then subtracting it
        # # Right Femur
        # m = util.Measurement('femur_l_Y', mset)
        # m.add_markerpair('RTC','RKJC')
        # m.add_bodyscale('femur_l','Y')

        # m = util.Measurement('femur_l_Z', mset)
        # m.add_markerpair('RLK','RMK')
        # m.add_bodyscale('femur_l','Z')
        # m.add_bodyscale('femur_l','X')

        # # Right Tibia
        # m = util.Measurement('tibia_l_Y', mset)
        # m.add_markerpair('RKJC','RAJC')
        # m.add_bodyscale('tibia_l','Y')

        # m = util.Measurement('tibia_l_Z', mset)
        # m.add_markerpair('RMA','RLA')
        # m.add_bodyscale('tibia_l','Z')
        # m.add_bodyscale('tibia_l','X')

        # # Right Foot
        # m = util.Measurement('foot_l_X', mset)
        # m.add_markerpair('RLH','RLT')
        # m.add_bodyscale('talus_l','X')
        # m.add_bodyscale('calcn_l','X')
        # m.add_bodyscale('toes_l','X')

        # m = util.Measurement('foot_l_Y', mset)
        # m.add_markerpair('RAJC','RAJC_P') #RLH and RLA
        # m.add_bodyscale('talus_l','Y')
        # m.add_bodyscale('calcn_l','Y')
        # m.add_bodyscale('toes_l','Y')

        # m = util.Measurement('foot_l_Z', mset)
        # m.add_markerpair('RMT','RLT')
        # m.add_bodyscale('talus_l','Z')
        # m.add_bodyscale('calcn_l','Z')
        # m.add_bodyscale('toes_l','Z')

    # IK step marker tasks
    ikts.add_ikmarkertask_bilateral('.Shoulder', True, 5.0)
    ikts.add_ikmarkertask_bilateral('.Clavicle', True, 3.0)
    ikts.add_ikmarkertask_bilateral('.Biceps', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.Elbow', True, 5.0)
    ikts.add_ikmarkertask_bilateral('.MElbow', True, 5.0)    
    ikts.add_ikmarkertask_bilateral('.Forearm', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.Wrist', True, 5.0)
    ikts.add_ikmarkertask_bilateral('.ASIS', True, 20.0)
    ikts.add_ikmarkertask('S2', True, 5.0)
    ikts.add_ikmarkertask_bilateral('.PSIS', True, 10.0)
    ikts.add_ikmarkertask_bilateral('_HJC', True, 30.0)
    ikts.add_ikmarkertask('R.TH1', False, 0.0)
    ikts.add_ikmarkertask('R.TH2', False, 0.0)
    ikts.add_ikmarkertask('R.TH3', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.Knee', True, 15.0)
    ikts.add_ikmarkertask_bilateral('.MKnee', True, 15.0)
    ikts.add_ikmarkertask('R.SH1', False, 0.0)
    ikts.add_ikmarkertask('R.SH2', False, 0.0)
    ikts.add_ikmarkertask('R.SH3', False, 0.0)
    ikts.add_ikmarkertask('R.SH4', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.Ankle', True, 25.0)
    ikts.add_ikmarkertask_bilateral('.MAnkle', True, 5.0)
    ikts.add_ikmarkertask_bilateral('.Toe', True, 1.0)
    ikts.add_ikmarkertask_bilateral('.MT5', True, 1.0)
    ikts.add_ikmarkertask_bilateral('.Heel', True, 1.0)
    ikts.add_ikmarkertask('L.TH1', False, 0.0)
    ikts.add_ikmarkertask('L.TH2', False, 0.0)
    ikts.add_ikmarkertask('L.TH3', False, 0.0)
    ikts.add_ikmarkertask('L.TH4', False, 0.0)
    ikts.add_ikmarkertask('L.SH1', False, 0.0)
    ikts.add_ikmarkertask('L.SH2', False, 0.0)
    ikts.add_ikmarkertask('L.SH3', False, 0.0)
    ikts.add_ikmarkertask('MidPSIS', True, 10.0)
    ikts.add_ikmarkertask('MidASIS', True, 10.0)
    ikts.add_ikmarkertask('MidClavicle', True, 10.0)
    ikts.add_ikmarkertask('MidPSIS', True, 10.0)
    ikts.add_ikmarkertask('MidPSIS', True, 10.0)
    ikts.add_ikmarkertask_bilateral('.AnkleJoint', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.KneeJoint', False, 0.0)
    ikts.add_ikmarkertask_bilateral('.HeelGround', True, 2.0)
    ikts.add_ikmarkertask_bilateral('.MT5Ground', True, 2.0)
    ikts.add_ikmarkertask_bilateral('.ToeGround', True, 2.0)
    ikts.add_ikmarkertask_bilateral('.AnkleJointGround', False, 0.0)

    # '''
    # IK step marker tasks old
        # ikts.add_ikmarkertask('C7', True, 5.0)
        # ikts.add_ikmarkertask_bilateral('AC', True, 3.0)
        # ikts.add_ikmarkertask('MidAC', True, 5.0)
        # ikts.add_ikmarkertask_bilateral('IC', True, 5.0)
        # ikts.add_ikmarkertask_bilateral('WST', True, 2.0)
        # ikts.add_ikmarkertask_bilateral('ASIS', True, 10.0)
        # ikts.add_ikmarkertask_bilateral('PSIS', True, 10.0)
        # ikts.add_ikmarkertask('MidIC', True, 5.0)
        # ikts.add_ikmarkertask('MidASIS', True, 5.0)
        # ikts.add_ikmarkertask('MidPSIS', True, 5.0)
        # ikts.add_ikmarkertask('MidPELV', True, 18.0)
        # ikts.add_ikmarkertask('RTC', True, 1.0)
        # ikts.add_ikmarkertask('LTC', False, 0.0)
        # ikts.add_ikmarkertask('MidTC', False, 0.0)
        # ikts.add_ikmarkertask('RSLT', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RSMT', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RIMT', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RILT', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RLK', True, 20.0)
        # ikts.add_ikmarkertask('RMK', True, 20.0)
        # ikts.add_ikmarkertask('RSPS', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RSAS', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RIPS', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RIAS', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RLA', True, 10.0)
        # ikts.add_ikmarkertask('RMA', True, 10.0)
        # # ikts.add_ikmarkertask('RHL', True, 5.0)
        # ikts.add_ikmarkertask('RLH', True, 2.0)
        # # ikts.add_ikmarkertask('RMH', True, 5.0)
        # ikts.add_ikmarkertask('RLT', True, 3.0)
        # ikts.add_ikmarkertask('RMT', True, 5.0)
        # ikts.add_ikmarkertask('RDTip', True, 5.0)
        # ikts.add_ikmarkertask('RHJC', False, 0.0)
        # ikts.add_ikmarkertask('RHJC_P', False, 0.0)
        # ikts.add_ikmarkertask('RKJC', True, 20.0)
        # ikts.add_ikmarkertask('RAJC', True, 20.0)
        # ikts.add_ikmarkertask('RAJC_P', True, 20.0)
        # # ikts.add_ikmarkertask('RHL_P', True, 5.0)
        # ikts.add_ikmarkertask('RLT_P', True, 10.0)
        # ikts.add_ikmarkertask('RMT_P', True, 10.0)
        # ikts.add_ikmarkertask('RMidT_P', True, 10.0)
        # ikts.add_ikmarkertask('RDAC', False, 0.0)               # dont use tracking markers in scaling
        # ikts.add_ikmarkertask('RPAC', False, 0.0)               # dont use tracking markers in scaling
        # # '''

def add_to_study(study):
    subject = study.add_subject(5, 112.4317)

    # make sure that the matfile is saved as version 7.3 or later
    ## TODO: write up a helper to make sure that the .mat file is v7.3 or later
    # mat_file_path = os.path.join(study.config['data_path'],'raw','matlab','s024.mat')
    # mat_file_path = os.path.join(study.config['data_path'],'segmented','subject024','subject024.mat')
    # mat_file = h5py.File(mat_file_path,'r')

    # model_check_path = os.path.join(study.config['data_path'],'loadedWalking_amySilder','assistloadwalk_simulations_of_experiments','experiments','subject05','noload')
    
    # check if we have a scaled model, or set up a scaling task for this subject
    # if os.path.isfile(os.path.join(model_check_path,'subject05_new4.osim')):
        # print('\nModel has been scaled and is ready for IK')
    # else:
    ## static trial for scaling
    static = subject.add_condition('static')
    static_trial = static.add_trial(1, omit_trial_dir=True)

    # `os.path.basename(__file__)` should be `subject024.py`.
    scale_setup_task = subject.add_task(osp.TaskScaleSetup,
            init_time=0.9,
            final_time=3.8,
            mocap_trial=static_trial,
            edit_setup_function=scale_setup_fcn,
            addtl_file_dep=['dodo.py', os.path.basename(__file__)])


    subject.add_task(osp.TaskScale,
            scale_setup_task=scale_setup_task,
            #scale_max_isometric_force=True,
            )

    subject.add_task(tasks.TaskScaleMuscleMaxIsometricForce)
    subject.scaled_model_fpath = os.path.join(subject.results_exp_path,
        '%s_scaled_Fmax.osim' % subject.name)


    # set up the noloaded condition with all trials. 
    noload = subject.add_condition('noload')
    noload_trial1 = noload.add_trial(num=1,
                                    metadata='C:/Users/JP/code/repos/Stanford/delplab/data/muscleEnergyModel/loadedWalking_amySilder/assistloadwalk_simulations_of_experiments/experiments/subject05/noload/free/metadata.yaml',
                                    omit_trial_dir=False,
                                    primary_leg=None,
                                    model_to_adjust_fpath=None,
                                    interval=None,
                                    gait_events=None)
    
    # noload_trial2 = noload.add_trial(2)
    print('did it work?')


    print('\nwell hello there\n')
    import pdb
    pdb.set_trace()

    # loaded condition
    loaded = subject.add_condition('loaded')
    loaded_trial = loaded.add_trial(1, # could add in the gait events later
                omit_trial_dir=True)



    ## slack condition
    slack = subject.add_condition('slack')
    mat_file_temp = mat_file['stp_000']['non_tme_r']['x']['raw']

    gait_events = dict()
    gait_events['right_strikes'] = mat_file_temp[:,0].tolist()
    gait_events['stride_times']= (mat_file_temp[:,1000]-mat_file_temp[:,0]).tolist()
    # gait_events['right_toeoffs'] = []
    # gait_events['left_strikes'] = []
    # gait_events['left_toeoffs'] = []

    slack_trial = slack.add_trial(1,
            gait_events=gait_events,
            omit_trial_dir=True,
            )

    """

    ## low condition
    low = subject.add_condition('low')
    mat_file_temp = mat_file['stp_025']['non_tme_r']['x']['raw']

    gait_events = dict()
    gait_events['right_strikes'] = mat_file_temp[:,0].tolist()
    gait_events['stride_times']= (mat_file_temp[:,1000]-mat_file_temp[:,0]).tolist()
    # gait_events['right_toeoffs'] = []
    # gait_events['left_strikes'] = []
    # gait_events['left_toeoffs'] = []

    low_trial = low.add_trial(1,
            gait_events=gait_events,
            omit_trial_dir=True,
            )

    ## med condition
    med = subject.add_condition('med')
    mat_file_temp = mat_file['stp_050']['non_tme_r']['x']['raw']

    gait_events = dict()
    gait_events['right_strikes'] = mat_file_temp[:,0].tolist()
    gait_events['stride_times']= (mat_file_temp[:,1000]-mat_file_temp[:,0]).tolist()
    # gait_events['right_toeoffs'] = []
    # gait_events['left_strikes'] = []
    # gait_events['left_toeoffs'] = []

    med_trial = med.add_trial(1,
            gait_events=gait_events,
            omit_trial_dir=True,
            )

    ## high condition
    high = subject.add_condition('high')
    mat_file_temp = mat_file['stp_075']['non_tme_r']['x']['raw']

    gait_events = dict()
    gait_events['right_strikes'] = mat_file_temp[:,0].tolist()
    gait_events['stride_times']= (mat_file_temp[:,1000]-mat_file_temp[:,0]).tolist()
    # gait_events['right_toeoffs'] = []
    # gait_events['left_strikes'] = []
    # gait_events['left_toeoffs'] = []

    high_trial = high.add_trial(1,
            gait_events=gait_events,
            omit_trial_dir=True,
            )

    ## max condition
    max = subject.add_condition('max')
    mat_file_temp = mat_file['stp_100']['non_tme_r']['x']['raw']

    gait_events = dict()
    gait_events['right_strikes'] = mat_file_temp[:,0].tolist()
    gait_events['stride_times']= (mat_file_temp[:,1000]-mat_file_temp[:,0]).tolist()
    # gait_events['left_strikes'] = []
    # gait_events['left_toeoffs'] = []

    max_trial = max.add_trial(1,
            gait_events=gait_events,
            omit_trial_dir=True,
            )

    ## add in the max isometric force calibration:




    # inverse kinematics/dynamics
    # helpers.generate_ik_id_tasks(slack_trial)
    # helpers.generate_ik_id_tasks(low_trial)
    # helpers.generate_ik_id_tasks(med_trial)
    # helpers.generate_ik_id_tasks(high_trial)
    # helpers.generate_ik_id_tasks(max_trial)

    # subject.add_task(tasks.TaskAvgPlot)



    # # muscle redundancy solver
    # # ========================
    # # first 3 go through the setup for the slack trial
    # # then the next 4 are the helper verision of those setups and post for other conditions
    
    # '''
    # Available cost functions:
    # Met - metabolic cost over the course of the trail
    # Default - activations squared. 
    # '''

    # # be sure that I have the parameters in the param_dict that I want for this subject
    
    # these are the slack kinematics on the slack case
    # mrs_setup_tasks = slack_trial.add_task_cycles(osp.TaskMRSDeGrooteSetup, param_dict=study.param_dict,
    #     cost='Met', use_filtered_id_results=True)
    # slack_trial.add_task_cycles(osp.TaskMRSDeGroote,
    #     setup_tasks=mrs_setup_tasks)
    # slack_trial.add_task_cycles(osp.TaskMRSDeGrootePost,
    #     setup_tasks=mrs_setup_tasks)
    
    
    # ## attempt to put everything together in one
    # mrs_setup_tasks = low_trial.add_task_cycles(osp.TaskMRSDeGrooteSetup, param_dict=study.param_dict,
    #     cost='Met', use_filtered_id_results=True)
    # # first the assisted kinematics without assistance
    # low_trial.add_task_cycles(osp.TaskMRSDeGroote,
    #     setup_tasks=mrs_setup_tasks)
    # low_trial.add_task_cycles(osp.TaskMRSDeGrootePost,
    #     setup_tasks=mrs_setup_tasks)
    # # now with assistance
    # helpers.generate_mrsmod_task(low_trial, mrs_setup_tasks, 'low')
    
    # mrs_setup_tasks = med_trial.add_task_cycles(osp.TaskMRSDeGrooteSetup, param_dict=study.param_dict,
    #     cost='Met', use_filtered_id_results=True)
    # # assisted kinematics with no assistance
    # med_trial.add_task_cycles(osp.TaskMRSDeGroote,
    #     setup_tasks=mrs_setup_tasks)
    # med_trial.add_task_cycles(osp.TaskMRSDeGrootePost,
    #     setup_tasks=mrs_setup_tasks)
    # # assisted kinematics with assistance
    # helpers.generate_mrsmod_task(med_trial, mrs_setup_tasks, 'med')

    # mrs_setup_tasks = high_trial.add_task_cycles(osp.TaskMRSDeGrooteSetup, param_dict=study.param_dict,
    #     cost='Met', use_filtered_id_results=True)
    # # assisted kinematics without assistance
    # high_trial.add_task_cycles(osp.TaskMRSDeGroote,
    #     setup_tasks=mrs_setup_tasks)
    # high_trial.add_task_cycles(osp.TaskMRSDeGrootePost,
    #     setup_tasks=mrs_setup_tasks)
    # # assisted kinematics with assistance
    # helpers.generate_mrsmod_task(high_trial, mrs_setup_tasks, 'high')

    # mrs_setup_tasks = max_trial.add_task_cycles(osp.TaskMRSDeGrooteSetup, param_dict=study.param_dict,
    #     cost='Met', use_filtered_id_results=True)
    # # assisted kinematics without assistance
    # max_trial.add_task_cycles(osp.TaskMRSDeGroote,
    #     setup_tasks=mrs_setup_tasks)
    # max_trial.add_task_cycles(osp.TaskMRSDeGrootePost,
    #     setup_tasks=mrs_setup_tasks)
    # # assisted kinematics with assistance
    # helpers.generate_mrsmod_task(max_trial, mrs_setup_tasks, 'max')
    """