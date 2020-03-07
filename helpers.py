import osimpipeline as osp
import tasks

def generate_mrsmod_task(trial, mrs_setup_tasks, force_level):

    if force_level == 'low': exo_force_level = 1
    elif force_level == 'med': exo_force_level = 2 
    elif force_level == 'high': exo_force_level = 3
    elif force_level == 'max': exo_force_level = 4  
          
    mrsflags = [
        "study='AnkleHipExosuit'",
        "exo_force_level=%s" % exo_force_level, # shift_exo_peaksforce_level,
        "data_path='%s'" % trial.study.config['data_path'],
        "subject='%s'" % trial.subject.name,
        "condition='%s'" % trial.condition.name,
        "shift_exo_peaks=%s" % trial.study.shift_exo_peaks,
        ]

    mrsmod_tasks = trial.add_task_cycles(
        osp.TaskMRSDeGrooteMod,
        force_level,
        'AnkleHipExosuit: apply low exosuit force condition',
        mrsflags,
        setup_tasks=mrs_setup_tasks
        )

    trial.add_task_cycles(tasks.TaskMRSDeGrooteModPost,
        setup_tasks=mrsmod_tasks)

def generate_ik_id_tasks(trial):

    ik_setup_task = trial.add_task(osp.TaskIKSetup)
    trial.add_task(osp.TaskIK, ik_setup_task)
    trial.add_task(osp.TaskIKPost, ik_setup_task,
        error_markers=trial.study.error_markers, side='r')

    trial.add_task(tasks.TaskFilterGRFs)
    trial.add_task(tasks.TaskUpdateGroundReactionColumnLabels)
    id_setup_task = trial.add_task(osp.TaskIDSetup, ik_setup_task)
    trial.add_task(osp.TaskID, id_setup_task)

    # Create Butterworth filter to process net joint moments from ID
    from scipy.signal import butter
    import math     
    
    grf_sample_freq = 2160 # Hz
    cutoff = 200 # Hz

    nyq = 0.5 * grf_sample_freq
    # 1080
    norm_cutoff = cutoff / nyq
    # .185

    b, a = butter(4, 0.15, btype='low', analog=False) # old is 4, 0.15
    butter_polys = (b, a)

    # need to find out where the best and right place to do this filtering and editing

    trial.add_task(osp.TaskIDPost, id_setup_task, 
        butter_polys=butter_polys,  # default is false, this turns off the joint moment filtering
        plot_primary_leg_only=True)

    
    