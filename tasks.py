import os
import itertools
import pyBTK as btk
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

import helpers

import opensim as osm
import osimpipeline as osp
from osimpipeline import utilities as util
from osimpipeline import postprocessing as pp
from osimpipeline import task as task
from matplotlib import colors as mcolors



colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

class TaskAvgPlot(task.SubjectTask):
    """Create a plot of the average kinematics for all conditions run for
    a subject
    TODO:
        fix side selection
        select for dofs
    """
    REGISTRY = []
    def __init__(self, subject):
        super(TaskAvgPlot, self).__init__(subject)
        self.name = '%s_avg_kinematics_plot' % (self.subject.name)
        self.doc = "Create a plot of average kinematics accross conditions."
        self.results_plot_path = os.path.join(
                self.study.config['results_path'], 'experiments',
                self.subject.rel_path)
        self.add_action([], [self.results_plot_path], self.create_plot)

    def create_plot(self, file_dep, target):
        fig = pl.figure(figsize=(7, 10))
        dims = (4, 2)

        cond_list=[]
        f_paths=[]
        base_f_path = self.results_plot_path
        for cond in self.subject.conditions:
            if cond.name != 'static':
                f_paths.append(os.path.join(base_f_path, cond.name, 'ik',
                'avg_joint_angles.txt'))
                cond_list.append(cond.name)
        self.cond_list=cond_list
        fig = pp.plot_lower_limb_kinematics_mod(self,f_paths, side = 'r', avg = True)
        f_name='%s\\%s_avg_joint_angles.pdf' % (base_f_path,self.subject.name)
        fig.savefig(f_name)
        pl.close(fig)

class TaskAverageAnkleKinematicsPlot(task.StudyTask):
    """Create a plot of the average kinematics for all conditions run
    across all subjects in the study.
    TODO:
        fix side selection
        select for dofs
    """
    REGISTRY = []
    def __init__(self, study):
        super(TaskAverageAnkleKinematicsPlot, self).__init__(study)
        self.name = '%s_average_ankle_kinematics_plot' % (self.study.name)
        self.doc = "Create a plot of average kinematics accross conditions."
        self.analysis_png_fpath_ank = os.path.join(
                self.study.config['analysis_path'], 
                'average_ankle_kinematics.png')
        self.analysis_pdf_fpath_ank = os.path.join(
                self.study.config['analysis_path'], 
                'average_ankle_kinematics.pdf')
        self.analysis_png_fpath_hip = os.path.join(
                self.study.config['analysis_path'], 
                'average_hip_kinematics.png')
        self.analysis_pdf_fpath_hip = os.path.join(
                self.study.config['analysis_path'], 
                'average_hip_kinematics.pdf')
        self.analysis_png_fpath_knee = os.path.join(
                self.study.config['analysis_path'], 
                'average_knee_kinematics.png')
        self.analysis_pdf_fpath_knee = os.path.join(
                self.study.config['analysis_path'], 
                'average_knee_kinematics.pdf')
        self.analysis_png_fpath_24_ank = os.path.join(
                self.study.config['analysis_path'], 
                '24_average_ank_kinematics.png')
        self.analysis_pdf_fpath_24_ank = os.path.join(
                self.study.config['analysis_path'], 
                '24_average_ank_kinematics.pdf') 
        self.analysis_png_fpath_77_ank = os.path.join(
                self.study.config['analysis_path'], 
                '77_average_ank_kinematics.png')
        self.analysis_pdf_fpath_77_ank = os.path.join(
                self.study.config['analysis_path'], 
                '77_average_ank_kinematics.pdf')
        self.analysis_png_fpath_88_ank = os.path.join(
                self.study.config['analysis_path'], 
                '88_average_ank_kinematics.png')
        self.analysis_pdf_fpath_88_ank = os.path.join(
                self.study.config['analysis_path'], 
                '88_average_ank_kinematics.pdf')
        self.analysis_png_fpath_112_ank = os.path.join(
                self.study.config['analysis_path'], 
                '112_average_ank_kinematics.png')
        self.analysis_pdf_fpath_112_ank = os.path.join(
                self.study.config['analysis_path'], 
                '112_average_ank_kinematics.pdf') 
        self.analysis_png_fpath_127_ank = os.path.join(
                self.study.config['analysis_path'], 
                '127_average_ank_kinematics.png')
        self.analysis_pdf_fpath_127_ank = os.path.join(
                self.study.config['analysis_path'], 
                '127_average_ank_kinematics.pdf') 
        self.analysis_png_fpath_128_ank = os.path.join(
                self.study.config['analysis_path'], 
                '128_average_ank_kinematics.png')
        self.analysis_pdf_fpath_128_ank = os.path.join(
                self.study.config['analysis_path'], 
                '128_average_ank_kinematics.pdf') 
                

        self.add_action([], 
                        [self.analysis_png_fpath_ank,
                         self.analysis_pdf_fpath_ank], 
                        self.create_plot)
        # self.add_action([], 
        #                 [self.analysis_png_fpath_knee,
        #                  self.analysis_pdf_fpath_knee], 
        #                 self.create_plot)
        # self.add_action([], 
        #                 [self.analysis_png_fpath_hip,
        #                  self.analysis_pdf_fpath_hip], 
        #                 self.create_plot)
        self.add_action([], 
                        [self.analysis_png_fpath_24_ank,
                         self.analysis_pdf_fpath_24_ank,
                         self.analysis_png_fpath_77_ank,
                         self.analysis_pdf_fpath_77_ank,
                         self.analysis_png_fpath_88_ank,
                         self.analysis_pdf_fpath_88_ank,
                         self.analysis_png_fpath_112_ank,
                         self.analysis_pdf_fpath_112_ank,
                         self.analysis_png_fpath_127_ank,
                         self.analysis_pdf_fpath_127_ank,
                         self.analysis_png_fpath_128_ank,
                         self.analysis_pdf_fpath_128_ank], 
                        self.create_plot)


    def create_plot(self, file_dep, target):

        fig = pl.figure(figsize=(6, 3))
        ax = fig.add_subplot(1, 1, 1)

        colors = {'slack': 'grey',
                  'low' : (0.662, 0.047, 0.913),
                  'med' : (0.094, 0.756, 0.905),
                  'high' : (0.094, 0.905, 0.266), 
                  'max' : (0.964, 0.607, 0.074),
                  }
        if len(target) <= 2:
            ankle_kin_conds = dict()
            for cond in self.study.subjects[0].conditions:
                if cond.name != 'static' and cond.name != 'scale':
                    ankle_kin = np.zeros(101).transpose()
                    for subject in self.study.subjects:
                        # kin_fpath = os.path.join(self.study.config['results_path'],
                        #     'experiments', subject.rel_path,  cond.name, 'ik', 
                        #     'avg_joint_angles.txt') # should this be turned to txt or pdf???!!!
                        kin_fpath = os.path.join(self.study.config['results_path'],
                            'experiments', subject.rel_path,  cond.name, 'ik', 
                            'avg_joint_angles.csv')
                        kin = pp.storage2numpy(kin_fpath)
                        if '_knee_' in target[0]:
                            ankle_kin += kin['knee_angle_r']
                        if '_hip_' in target[0]:
                            ankle_kin += kin['hip_flexion_r']
                        if '_ankle_' in target[0]:
                            ankle_kin += kin['ankle_angle_r']

                    ankle_kin_conds[cond.name] = ankle_kin / len(self.study.subjects)

                    per_gc = range(101)
                    ax.plot(per_gc, ankle_kin_conds[cond.name], 
                        color=colors[cond.name],label=cond.name)

            # ax.grid()
            ax.legend()
            ax.set_axisbelow(True)
            ax.set_xlabel('percent gait cycle')
            ax.set_ylabel('angle (degrees)')
            ax.set_yticks(np.linspace(-30,20,11))
            ax.set_yticklabels(np.linspace(-30,20,11))
            fig.tight_layout()
            fig.savefig(target[0], dpi=600)
            fig.savefig(target[1])
            pl.close(fig)

        if len(target) > 2:
            k = 0
            for i, subject in enumerate(self.study.subjects):
                # new figure for each subject
                print(i)
                print(subject)
                sub_fig = pl.figure(figsize=(6, 3))
                sub_ax = sub_fig.add_subplot(1, 1, 1)
                ankle_kin_conds_test = dict()

                for cond in self.study.subjects[0].conditions:
                    if cond.name != 'static' and cond.name != 'scale':
                        ankle_kin = np.zeros(101).transpose()

                        # kin_fpath = os.path.join(self.study.config['results_path'],
                        #     'experiments', subject.rel_path,  cond.name, 'ik', 
                        #     'avg_joint_angles.txt') # should this be turned to txt or pdf???!!!
                        kin_fpath = os.path.join(self.study.config['results_path'],
                            'experiments', subject.rel_path,  cond.name, 'ik', 
                            'avg_joint_angles.csv')
                        kin = pp.storage2numpy(kin_fpath)
                        if '_knee_' in target[0]:
                            ankle_kin += kin['knee_angle_r']
                        if '_hip_' in target[0]:
                            ankle_kin += kin['hip_flexion_r']
                        if '_ank_' in target[0]:
                            ankle_kin += kin['ankle_angle_r']

                        ankle_kin_conds_test[cond.name] = ankle_kin #/ len(self.study.subjects)

                        per_gc = range(101)
                        sub_ax.plot(per_gc, ankle_kin_conds_test[cond.name], 
                            color=colors[cond.name],label=cond.name)

                sub_ax.grid()
                sub_ax.legend()
                sub_ax.set_axisbelow(True)
                sub_ax.set_xlabel('percent gait cycle')
                sub_ax.set_ylabel('angle (degrees)')
                sub_ax.set_yticks(np.linspace(-30,20,11))
                sub_ax.set_yticklabels(np.linspace(-30,20,11))
                sub_fig.tight_layout()
                sub_fig.savefig(target[k], dpi=600)
                sub_fig.savefig(target[k+1])
                pl.close(fig)
                k += 2 

class TaskModifyExperimentalData(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study):

        super(TaskModifyExperimentalData, self).__init__(study)                 #!!! idk what this is
        import pdb
        # pdb.set_trace()
        self.name = '_'.join([study.name, 'modify_exp_data'])
        self.doc = ('Modify experimental data in C3D files (or any files) to prepare for '
            'OpenSim pipeline.')
        self.preprocess_path = os.path.join(self.study.config['preprocess_path'],study.name,'preprocess')
        self.repo_path = self.study.config['repo_path']
        self.cond_map = self.study.config['cond_map']
        self.mocap_path = os.path.join(self.study.config['motion_capture_data_path'],study.name,'mocap')

        self.actions += [#self.modify_data_files,
                         self.copy_data_files_to_osim]                          #!!! does this call these, or not?



    def copy_data_files_to_osim(self):
        """Copy modified data files into "osim" directory under main data
        directory. "osim" is the motion_capture_data_path where files will
        be copied from to the "results" folder for simulations.
        """
        import pdb
        pdb.set_trace()
        import shutil
        for subj in self.study.all_subjects:
            for k, v in self.cond_map.items():
                for ext in ['.trc', '.mot']:
                    mocap_fpath = os.path.join(self.mocap_path, subj, v + ext)
                    prepro_fpath = os.path.join(self.preprocess_path,
                        subj, v + ext)

                    to_dir = os.path.split(mocap_fpath)[0]
                    if not os.path.exists(to_dir): os.makedirs(to_dir)
                    shutil.copyfile(prepro_fpath, mocap_fpath)



    def modify_data_files(self):
        """Make modifications to data files after being extracted from C3D:
        delete unused markers in TRC files, set data capture frequency, etc.
        """
        print('\nin here')
        import pdb
        pdb.set_trace()
        # loop through the subjects
        for subj in self.study.all_subjects:
            print(subj)
            # loop through the loaded or unloaded conditions (sometimes staticredo)
            loaddir = os.path.join(self.preprocess_path, subj)
            loaddirlist = []
            for name in os.listdir(loaddir):
                if os.path.isdir(os.path.join(loaddir, name)):
                    loaddirlist.append(name)
                    print('-->',name)
                    # loop through the conditions
                    conddir = os.path.join(loaddir, name)
                    conddirlist = []
                    for each in os.listdir(conddir):
                        if os.path.isdir(os.path.join(conddir, each)):
                            conddirlist.append(each)
                            print('-->-->',each)
                            # loop through trial numbers
                            trialdir = os.path.join(conddir, each)
                            trialdirlist = []
                            # print(trialdir)        
                            # now for trials
                            for trial in os.listdir(trialdir):
                                if os.path.isdir(os.path.join(trialdir,trial)):
                                    trialdirlist.append(trial)
                                    print('-->-->-->',trial)

                                    mocapdir = os.path.join(trialdir,trial,'expdata','motion_capture.trc')
                                    # now what?
                                    # pdb.set_trace()

                                    helpers.filemodifier(self,mocapdir,subj,name,each,trial)
                                    # pdb.set_trace()
            print('\n')
            pdb.set_trace()



    def add_virtual_markers(self,acq):
        import pdb
        pdb.set_trace()
        # Midpoint between shoulder markers
        # acq.RemovePoint("MidAC")
        RAC = acq.getDependentColumn("R.Shoulder")
        LAC = acq.getDependentColumn("L.Shoulder")
        pdb.set_trace()
        MidAC = btk.btkPoint("MidAC", acq.GetPointFrameNumber())
        MidAC.SetValues( (RAC.GetValues() + LAC.GetValues() )/2.0 )
        acq.AppendPoint(MidAC)
        '''
        # midpoint of clavicle
        Rclav = acq.GetPoint("R.Clavicle")
        Lclav = acq.GetPoint("L.Clavicle")
        MidClav = btk.btkPoint("MidClav", acq.GetPointFrameNumber())
        MidClav.SetValues( (Rclav.GetValues() + Lclav.GetValues() )/2.0 )
        acq.AppendPoint(MidClav)




        # Midpoint between iliac crests
        acq.RemovePoint("MidIC")
        RIC = acq.GetPoint("RIC")
        LIC = acq.GetPoint("LIC")
        MidIC = btk.btkPoint("MidIC", acq.GetPointFrameNumber())
        MidIC.SetValues( ( RIC.GetValues() + LIC.GetValues() )/2.0 )
        acq.AppendPoint(MidIC)

        # Midpoint between trochanter markers
        acq.RemovePoint("MidTC")
        RTC = acq.GetPoint("RTC")
        LTC = acq.GetPoint("LTC")
        MidTC = btk.btkPoint("MidTC", acq.GetPointFrameNumber())
        MidTC.SetValues( (RTC.GetValues() + LTC.GetValues() )/2.0 )
        acq.AppendPoint(MidTC)


        # Midpoint between ASIS markers
        acq.RemovePoint("MidASIS")
        RASIS = acq.GetPoint("RASIS")
        LASIS = acq.GetPoint("LASIS")
        MidASIS = btk.btkPoint("MidASIS", acq.GetPointFrameNumber())
        MidASIS.SetValues( ( RASIS.GetValues() + LASIS.GetValues() )/2.0 )
        acq.AppendPoint(MidASIS)

        # Midpoint between PSIS markers
        acq.RemovePoint("MidPSIS")
        RPSIS = acq.GetPoint("RPSIS")
        LPSIS = acq.GetPoint("LPSIS")
        MidPSIS = btk.btkPoint("MidPSIS", acq.GetPointFrameNumber())
        MidPSIS.SetValues( ( RPSIS.GetValues() + LPSIS.GetValues() )/2.0 )
        acq.AppendPoint(MidPSIS)

        # Mid pelvis virtual marker
        acq.RemovePoint("MidPELV")
        MidPELV = btk.btkPoint("MidPELV", acq.GetPointFrameNumber())
        MidPELV.SetValues( ( MidPSIS.GetValues() + MidASIS.GetValues() )/2.0 )
        acq.AppendPoint(MidPELV)

        # Projected HJC location (set y-position to MidIC's y-position)
        acq.RemovePoint("RHJC_P")
        RHJC_P = btk.btkPoint("RHJC_P", acq.GetPointFrameNumber())
        RHJC = acq.GetPoint("RHJC")
        RHJCp = RHJC.GetValues()
        RHJCp[:,1] = MidIC.GetValues()[:,1]
        RHJC_P.SetValues(RHJCp)
        acq.AppendPoint(RHJC_P)

        # Right knee joint center
        acq.RemovePoint("RKJC")
        RLK = acq.GetPoint("RLK")
        RMK = acq.GetPoint("RMK")
        RKJC = btk.btkPoint("RKJC", acq.GetPointFrameNumber())
        RKJC.SetValues( ( RLK.GetValues() + RMK.GetValues() )/2.0 )
        acq.AppendPoint(RKJC)

        # Right ankle joint center
        acq.RemovePoint("RAJC")
        RLA = acq.GetPoint("RLA")
        RMA = acq.GetPoint("RMA")
        RAJC = btk.btkPoint("RAJC", acq.GetPointFrameNumber())
        RAJC.SetValues( ( RLA.GetValues() + RMA.GetValues() )/2.0 )
        acq.AppendPoint(RAJC)

        # Right ankle joint center projected onto the floor
        acq.RemovePoint("RAJC_P")
        RAJC_P = btk.btkPoint("RAJC_P", acq.GetPointFrameNumber())
        RAJCp = RAJC.GetValues()
        RAJCp[:,1] = 0 # set vertical axis component to zero
        RAJC_P.SetValues(RAJCp)
        acq.AppendPoint(RAJC_P)

        # Right heel marker projected onto the floor
        # acq.RemovePoint("RHL_P")
        # RHL_P = btk.btkPoint("RHL_P", acq.GetPointFrameNumber())
        # RHL = acq.GetPoint("RHL")
        # RHLp = RHL.GetValues()
        # RHLp[:,1] = 0 # set vertical axis component to zero
        # RHL_P.SetValues(RHL)
        # acq.AppendPoint(RHL_P)

        # RLT marker projected onto the floor
        acq.RemovePoint("RLT_P")
        RLT_P = btk.btkPoint("RLT_P", acq.GetPointFrameNumber())
        RLT = acq.GetPoint("RLT")
        RLTp = RLT.GetValues()
        RLTp[:,1] = 0 # set vertical axis component to zero
        RLT_P.SetValues(RLTp)
        acq.AppendPoint(RLT_P)

        # RMT marker projected onto the floor
        acq.RemovePoint("RMT_P")
        RMT_P = btk.btkPoint("RMT_P", acq.GetPointFrameNumber())
        RMT = acq.GetPoint("RMT")
        RMTp = RMT.GetValues()
        RMTp[:,1] = 0 # set vertical axis component to zero
        RMT_P.SetValues(RMTp)
        acq.AppendPoint(RMT_P)

        # Right mid toe marker projected onto the floor
        acq.RemovePoint("RMidT_P")
        RMidT_P = btk.btkPoint("RMidT_P", acq.GetPointFrameNumber())
        RMidT_P.SetValues( (RLTp+RMTp)/2 )
        acq.AppendPoint(RMidT_P)
        '''








        pdb.set_trace()
        ## old stuff
            # # Midpoint between trochanter markers
            # acq.RemovePoint("MidTC")
            # RTC = acq.GetPoint("RTC")
            # LTC = acq.GetPoint("LTC")
            # MidTC = btk.btkPoint("MidTC", acq.GetPointFrameNumber())
            # MidTC.SetValues( (RTC.GetValues() + LTC.GetValues() )/2.0 )
            # acq.AppendPoint(MidTC)

            # # Midpoint between shoulder markers
            # acq.RemovePoint("MidAC")
            # RAC = acq.GetPoint("RAC")
            # LAC = acq.GetPoint("LAC")
            # MidAC = btk.btkPoint("MidAC", acq.GetPointFrameNumber())
            # MidAC.SetValues( (RAC.GetValues() + LAC.GetValues() )/2.0 )
            # acq.AppendPoint(MidAC)

            # # Midpoint between iliac crests
            # acq.RemovePoint("MidIC")
            # RIC = acq.GetPoint("RIC")
            # LIC = acq.GetPoint("LIC")
            # MidIC = btk.btkPoint("MidIC", acq.GetPointFrameNumber())
            # MidIC.SetValues( ( RIC.GetValues() + LIC.GetValues() )/2.0 )
            # acq.AppendPoint(MidIC)

            # # Midpoint between ASIS markers
            # acq.RemovePoint("MidASIS")
            # RASIS = acq.GetPoint("RASIS")
            # LASIS = acq.GetPoint("LASIS")
            # MidASIS = btk.btkPoint("MidASIS", acq.GetPointFrameNumber())
            # MidASIS.SetValues( ( RASIS.GetValues() + LASIS.GetValues() )/2.0 )
            # acq.AppendPoint(MidASIS)

            # # Midpoint between PSIS markers
            # acq.RemovePoint("MidPSIS")
            # RPSIS = acq.GetPoint("RPSIS")
            # LPSIS = acq.GetPoint("LPSIS")
            # MidPSIS = btk.btkPoint("MidPSIS", acq.GetPointFrameNumber())
            # MidPSIS.SetValues( ( RPSIS.GetValues() + LPSIS.GetValues() )/2.0 )
            # acq.AppendPoint(MidPSIS)

            # # Mid pelvis virtual marker
            # acq.RemovePoint("MidPELV")
            # MidPELV = btk.btkPoint("MidPELV", acq.GetPointFrameNumber())
            # MidPELV.SetValues( ( MidPSIS.GetValues() + MidASIS.GetValues() )/2.0 )
            # acq.AppendPoint(MidPELV)

            # # Projected HJC location (set y-position to MidIC's y-position)
            # acq.RemovePoint("RHJC_P")
            # RHJC_P = btk.btkPoint("RHJC_P", acq.GetPointFrameNumber())
            # RHJC = acq.GetPoint("RHJC")
            # RHJCp = RHJC.GetValues()
            # RHJCp[:,1] = MidIC.GetValues()[:,1]
            # RHJC_P.SetValues(RHJCp)
            # acq.AppendPoint(RHJC_P)

            # # Right knee joint center
            # acq.RemovePoint("RKJC")
            # RLK = acq.GetPoint("RLK")
            # RMK = acq.GetPoint("RMK")
            # RKJC = btk.btkPoint("RKJC", acq.GetPointFrameNumber())
            # RKJC.SetValues( ( RLK.GetValues() + RMK.GetValues() )/2.0 )
            # acq.AppendPoint(RKJC)

            # # Right ankle joint center
            # acq.RemovePoint("RAJC")
            # RLA = acq.GetPoint("RLA")
            # RMA = acq.GetPoint("RMA")
            # RAJC = btk.btkPoint("RAJC", acq.GetPointFrameNumber())
            # RAJC.SetValues( ( RLA.GetValues() + RMA.GetValues() )/2.0 )
            # acq.AppendPoint(RAJC)

            # # Right ankle joint center projected onto the floor
            # acq.RemovePoint("RAJC_P")
            # RAJC_P = btk.btkPoint("RAJC_P", acq.GetPointFrameNumber())
            # RAJCp = RAJC.GetValues()
            # RAJCp[:,1] = 0 # set vertical axis component to zero
            # RAJC_P.SetValues(RAJCp)
            # acq.AppendPoint(RAJC_P)

            # # Right heel marker projected onto the floor
            # # acq.RemovePoint("RHL_P")
            # # RHL_P = btk.btkPoint("RHL_P", acq.GetPointFrameNumber())
            # # RHL = acq.GetPoint("RHL")
            # # RHLp = RHL.GetValues()
            # # RHLp[:,1] = 0 # set vertical axis component to zero
            # # RHL_P.SetValues(RHL)
            # # acq.AppendPoint(RHL_P)

            # # RLT marker projected onto the floor
            # acq.RemovePoint("RLT_P")
            # RLT_P = btk.btkPoint("RLT_P", acq.GetPointFrameNumber())
            # RLT = acq.GetPoint("RLT")
            # RLTp = RLT.GetValues()
            # RLTp[:,1] = 0 # set vertical axis component to zero
            # RLT_P.SetValues(RLTp)
            # acq.AppendPoint(RLT_P)

            # # RMT marker projected onto the floor
            # acq.RemovePoint("RMT_P")
            # RMT_P = btk.btkPoint("RMT_P", acq.GetPointFrameNumber())
            # RMT = acq.GetPoint("RMT")
            # RMTp = RMT.GetValues()
            # RMTp[:,1] = 0 # set vertical axis component to zero
            # RMT_P.SetValues(RMTp)
            # acq.AppendPoint(RMT_P)

            # # Right mid toe marker projected onto the floor
            # acq.RemovePoint("RMidT_P")
            # RMidT_P = btk.btkPoint("RMidT_P", acq.GetPointFrameNumber())
            # RMidT_P.SetValues( (RLTp+RMTp)/2 )
            # acq.AppendPoint(RMidT_P)

class TaskCopyMotionCaptureData(osp.TaskCopyMotionCaptureData):
    REGISTRY = []
    def __init__(self, study, Slack=None, Low=None, Med=None, High=None,
        Max=None):
        regex_replacements = list()

        for subj in study.subjects:

            for datastr, condname in [
                    ('static', 'static'),
                    ('stp000', 'slack'),
                    ('stp025', 'low'),
                    ('stp050', 'med'),
                    ('stp075', 'high'),
                    ('stp100', 'max')]:

                # Marker trajectories.
                regex_replacements.append(
                    (
                        os.path.join(subj.name,
                            '%s.trc' % (datastr)).replace('\\','\\\\'),
                        os.path.join('experiments',
                            subj.name, condname, 'expdata',
                            'marker_trajectories.trc').replace('\\','\\\\')
                        ))
                # Ground reaction.
                regex_replacements.append(
                    (
                        os.path.join(subj.name,
                            '%s.mot' % (datastr)).replace(
                                '\\','\\\\'),
                        os.path.join('experiments', subj.name, condname,
                            'expdata','ground_reaction_orig.mot').replace(
                                '\\','\\\\')
                        ))

        super(TaskCopyMotionCaptureData, self).__init__(study,
                regex_replacements)

class TaskFilterGRFs(osp.TrialTask):
    REGISTRY = []
    def __init__(self, trial):
        super(TaskFilterGRFs, self).__init__(trial)
        self.name = trial.id + '_filter_grfs'
        self.add_action(
            [os.path.join(trial.expdata_path, 'ground_reaction_orig.mot')],
            [os.path.join(trial.expdata_path, 'ground_reaction_filt.mot')],
            self.filter_grfs)

    def filter_grfs(self, file_dep, target):
        from scipy.signal import butter
        from scipy.signal import filtfilt

        data = util.storage2numpy(file_dep[0])
        data_unfilt = np.array(data)
        data_filt = np.array(data)
        names = data_unfilt.dtype.names

        # Create a 4th-order, 10 Hz, low-pass filter
        b, a = butter(4, 0.01, btype='low', analog=False) # middle field is the frequency -> less aggressive!!!

        # Filter data and create new storage file
        for name in names:
            if not name == 'time':
                data_filt[name] = data_unfilt[name]
                # data_filt[name] = filtfilt(b, a, data_unfilt[name]) # this filters GRF

            ## we always want to do this, but may not always want to filter things
            # Convert torques from N-mm to N-m
            if 'torque' in name:
                data_filt[name] /= 1000.0 
        
        util.ndarray2storage(data_filt, target[0])

class TaskUpdateGroundReactionColumnLabels(osp.TrialTask):  # This is just taking the GRF file and renaming the columns so that they match what ID tool wants
    REGISTRY = []
    def __init__(self, trial):
        super(TaskUpdateGroundReactionColumnLabels, self).__init__(trial)
        self.name = trial.id + '_update_grf_column_labels'
        self.add_action(
                [os.path.join(trial.expdata_path, 'ground_reaction_filt.mot')],  # this can be ground_reaction_filt for the filtered grfs
                [trial.ground_reaction_fpath],
                self.dispatch)
    def dispatch(self, file_dep, target):
        import re
        data = pp.storage2numpy(file_dep[0])
        new_names = list()
        for name in data.dtype.names:
            if name == 'time':
                new_name = name
            elif name.startswith('1_'):
                new_name = re.sub('1_ground_(.*)_(.*)', 'ground_\\1_l_\\2',
                        name)
            elif name.startswith('2_'):
                new_name = re.sub('2_ground_(.*)_(.*)', 'ground_\\1_r_\\2',
                        name)
            else:
                Exception("Not a valid GRF column header name.")

            new_names.append(new_name)
        data.dtype.names = new_names
        util.ndarray2storage(data, target[0])

class TaskGenericModelSetup(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, default_model_fpath,
                 muscle_name_map=None, marker_map=None):
        super(TaskGenericModelSetup, self).__init__(study)
        self.name = '_'.join([study.name, 'generic_model_setup'])
        self.doc = 'Setup generic model to fit study needs.'
        self.muscle_name_map = muscle_name_map
        self.marker_map = marker_map
        self.default_model_fpath = default_model_fpath
        self.modified_generic_model_fpath = \
            self.study.source_generic_model_fpath

        self.actions += [self.load_model]

        if muscle_name_map:
            self.actions += [self.rename_muscles]

        self.actions += [self.set_clamped_joints]
        if marker_map:
            self.actions += [self.create_markersets]

        self.actions += [self.add_metabolics_probes,
                         self.print_modified_model]
        ##!!! make the whole clamping thing here!!!

    def load_model(self):
        self.model = osm.Model(self.default_model_fpath)

    def rename_muscles(self):
        muscle_set = self.model.updMuscles()

        for k, v in self.muscle_name_map.items():
            # Only right leg muscles
            if muscle_set.contains("%s_r" % k):
                muscle_set.get("%s_r" % k).setName("%s_r" % v)

    def set_clamped_joints(self):
        # this is going to be used to clamp any of the joints in the model
        # made because the knee in the model is not clamped, but should be?
        #!!! need to edit the modified one to clamp the knee
        #!!! need to add the action to the init above
        print ("\nhello need to finish the clamping editing function!!!\n")
        


    def create_markersets(self):

        marker_set = self.model.updMarkerSet()

        # There are three possible use cases for the marker_map dictionary:
        #     1) Deleting an unused marker
        #           'delete_marker': 'MARKER_TO_DELETE_NAME'
        #
        #     2) Adjust an exisiting marker
        #           'adjust_marker': ('CURRENT_NAME','NEW_NAME',
        #                             [x_loc_body,y_loc_body,z_loc_body],
        #                             'body_name',fixedBool)
        #
        #     3) Adding a new marker
        #           'add_marker': ('NAME',[x_loc_body,y_loc_body,z_loc_body],
        #                          'body_name',fixedBool)
        #

        for k, values in self.marker_map.items():
            for v in values:
                if k == 'delete_marker':
                    marker_set.remove(marker_set.get(v))

                elif k == 'adjust_marker':
                    offset = osm.Vec3(v[2][0],v[2][1],v[2][2])
                    body = self.model.getBodySet().get(v[3])

                    marker_set.get(v[0]).setBodyName(v[3])
                    marker_set.get(v[0]).changeBody(body)
                    marker_set.get(v[0]).setOffset(offset)
                    marker_set.get(v[0]).setFixed(v[4])
                    marker_set.get(v[0]).setName(v[1])

                elif k == 'add_marker':
                    marker = osm.Marker()
                    body = self.model.getBodySet().get(v[2])
                    offset = osm.Vec3(v[1][0],v[1][1],v[1][2])

                    marker.setFixed(v[3])
                    marker.setName(v[0])
                    marker.setBodyName(v[2])
                    marker.changeBody(body)
                    marker.setOffset(offset)

                    marker_set.cloneAndAppend(marker)

        marker_set.printToXML('templates/scale/prescale_markerset.xml')






    def print_modified_model(self):
        self.model.printToXML(self.modified_generic_model_fpath)
        
    def add_metabolics_probes(self):
        util.add_metabolics_probes(self.model, twitch_ratio_set='gait2392')

class TaskMRSDeGrooteModPost(osp.TaskMRSDeGrooteModPost):
    REGISTRY = []
    def __init__(self, trial, mrsmod_task, **kwargs):
        super(TaskMRSDeGrooteModPost, self).__init__(trial, mrsmod_task, 
            **kwargs)
        self.cost = mrsmod_task.mrs_setup_task.cost
        self.cycle = mrsmod_task.mrs_setup_task.tricycle

        # "slack" condition solution
        self.slack_results_output_fpath = os.path.join(
            trial.subject.results_exp_path, 'slack', 'mrs', self.cycle.name,
            '' if self.cost == 'Default' else self.cost, 
            '%s_%s_slack_%s_mrs.mat' % (self.study.name, self.subject.name, 
            self.cycle.name))

        self.add_action([self.slack_results_output_fpath,
                         self.results_output_fpath],
                         [os.path.join(self.path, 'metabolic_reductions.pdf')],
                         self.plot_metabolic_reductions)

        self.add_action([self.slack_results_output_fpath,
                         self.results_output_fpath],
                         [os.path.join(self.path, 
                            'muscle_activity_reductions.pdf')],
                         self.plot_muscle_activity_reductions)

    def plot_joint_moment_breakdown(self, file_dep, target):
        # Load mat file fields
        muscle_names = util.hdf2list(file_dep[0], 'MuscleNames', type=str)
        dof_names = util.hdf2list(file_dep[0],'DatStore/DOFNames',type=str)
        num_dofs = len(dof_names)
        num_muscles = len(muscle_names)
        joint_moments_exp = util.hdf2numpy(file_dep[0], 'DatStore/T_exp')
        tendon_forces = util.hdf2numpy(file_dep[0], 'TForce')
        exp_time = util.hdf2numpy(file_dep[0], 'DatStore/time').transpose()[0]
        time = util.hdf2numpy(file_dep[0], 'Time').transpose()[0]
        moment_arms_exp = util.hdf2numpy(file_dep[0], 'DatStore/dM').transpose()

        # Clip large tendon forces at final time
        from warnings import warn
        for imusc in range(len(muscle_names)):
            tendon_force = tendon_forces[:,imusc]
            if (tendon_force[-1] > 10*tendon_force[-2]):
                tendon_force[-1] = tendon_force[-2]
                tendon_forces[:,imusc] = tendon_force
                warn('WARNING: large %s tendon force at final time. '
                    'Clipping...' % muscle_names[imusc])

        # Get device torques
        device_torques = list()
        device_names = list()
        device_colors = list()
        act_torques = util.hdf2pandas(file_dep[0], 
            'DatStore/ExoTorques_Act', labels=dof_names)
        device_torques.append(act_torques)
        device_names.append('active')
        device_colors.append('green')

        # Interpolate to match solution time
        from scipy.interpolate import interp1d
        ma_shape = (len(time), moment_arms_exp.shape[1], 
            moment_arms_exp.shape[2])
        moment_arms = np.empty(ma_shape)
        for i in range(moment_arms_exp.shape[2]):
            func_moment_arms_interp = interp1d(exp_time, 
                moment_arms_exp[:,:,i].squeeze(), axis=0)
            moment_arms[:,:,i] = func_moment_arms_interp(time)

        func_joint_moments_interp = interp1d(exp_time, joint_moments_exp,
            axis=0)
        joint_moments = func_joint_moments_interp(time)

        # Generate plots
        pp.plot_joint_moment_breakdown(time, joint_moments, tendon_forces,
            moment_arms, dof_names, muscle_names, target[0], target[1],
            mass=self.subject.mass, ext_moments=device_torques,
            ext_names=device_names, ext_colors=device_colors)



    def plot_metabolic_reductions(self, file_dep, target):
        # Load mat file fields from original, "no-mod" solution
        muscle_names = util.hdf2list(file_dep[0], 'MuscleNames', type=str)
        num_muscles = len(muscle_names)
        mrs_whole_body_metabolic_rate = util.hdf2pandas(file_dep[0], 
            'DatStore/MetabolicRate/whole_body')
        mrs_muscle_metabolic_rates = util.hdf2pandas(file_dep[0],
            'DatStore/MetabolicRate/individual_muscles', labels=muscle_names)

        # Load mat file fields from modified solution
        muscle_names = util.hdf2list(file_dep[1], 'MuscleNames', type=str)
        num_muscles = len(muscle_names)

        mrsmod_whole_body_metabolic_rate = util.hdf2pandas(file_dep[1], 
            'DatStore/MetabolicRate/whole_body')
        mrsmod_muscle_metabolic_rates = util.hdf2pandas(file_dep[1],
            'DatStore/MetabolicRate/individual_muscles', labels=muscle_names)


        reductions = list()
        reduc_names = list()
        colors = list()

        for musc in muscle_names:
            muscle_reduction = 100.0 * ((mrsmod_muscle_metabolic_rates[musc] -
                mrs_muscle_metabolic_rates[musc]) / 
                mrs_muscle_metabolic_rates[musc])
            reductions.append(muscle_reduction)
            reduc_names.append(musc)
            colors.append('b')

        whole_body_reduction = 100.0 * (mrsmod_whole_body_metabolic_rate - 
            mrs_whole_body_metabolic_rate) / mrs_whole_body_metabolic_rate
        reductions.append(whole_body_reduction[0])
        reduc_names.append('whole_body')
        colors.append('r')

        # Plot metabolic reductions
        fig = pl.figure(figsize=(40,8.5))
        ax = fig.add_subplot(1,1,1)
        pos = np.arange(len(muscle_names)+1)

        reduc_array = np.array(reductions)
        # print reduction_array
        reductions_array= []
        for reduc in reduc_array:
        	reductions_array.append(reduc[0])

        # ax.bar(pos, reductions, align='center', color=colors)
        ax.bar(pos, reductions_array, align='center', color=colors)        
        ax.set_xticks(pos)
        ax.set_xticklabels(reduc_names, fontsize=10)
        ax.set_yticks(np.arange(-100,105,5))
        for label in ax.yaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
        ax.set_ylabel('Percent Change in Metabolic Rate', fontsize=12)
        ax.grid(which='both', axis='both', linestyle='--')
        ax.set_axisbelow(True)

        fig.tight_layout()
        fig.savefig(target[0])
        pl.close(fig)

    def plot_muscle_activity_reductions(self, file_dep, target):
        # Load mat file fields from original, "no-mod" solution

        muscle_names = util.hdf2list(file_dep[0], 'MuscleNames', type=str)
        num_muscles = len(muscle_names)

        mrs_excitations = util.hdf2pandas(file_dep[0], 
            'MExcitation', labels=muscle_names)
        mrs_activations = util.hdf2pandas(file_dep[0],
            'MActivation', labels=muscle_names)

        # Load mat file fields from modified solution
        muscle_names = util.hdf2list(file_dep[1], 'MuscleNames', type=str)
        num_muscles = len(muscle_names)

        mrsmod_excitations = util.hdf2pandas(file_dep[1], 
            'MExcitation', labels=muscle_names)
        mrsmod_activations = util.hdf2pandas(file_dep[1],
            'MActivation', labels=muscle_names)


        exc_reductions = list()
        act_reductions = list()
        reduc_names = list()
        exc_colors = list()
        act_colors = list()

        from matplotlib import colors as mcolors
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

        # Individual muscles
        for musc in muscle_names:
            reduc_names.append(musc)
            
            mrs_exc = mrs_excitations[musc]
            mrsmod_exc = mrsmod_excitations[musc]
            diff_exc = mrsmod_exc - mrs_exc
            reduc_exc = 100.0 * (np.nansum(diff_exc) / np.nansum(mrs_exc))
            exc_reductions.append(reduc_exc)
            exc_colors.append(colors['khaki'])
            
            mrs_act = mrs_activations[musc]
            mrsmod_act = mrsmod_activations[musc]
            diff_act = mrsmod_act - mrs_act
            reduc_act = 100.0 * (np.nansum(diff_act) / np.nansum(mrs_act))
            act_reductions.append(reduc_act)
            act_colors.append(colors['palegreen'])



        # Whole body
        reduc_names.append('whole_body')       
        whole_reduc_exc = sum(exc_reductions)
        exc_reductions.append(whole_reduc_exc)
        exc_colors.append('gold')

        whole_reduc_act = sum(act_reductions)
        act_reductions.append(whole_reduc_act)
        act_colors.append('seagreen')
        

        # Plot activity reductions
        fig = pl.figure(figsize=(40,8.5))
        ax = fig.add_subplot(1,1,1)
        pos = np.arange(len(muscle_names)+1)
        width = 0.4


        bar1 = ax.bar(pos, exc_reductions, width, color=exc_colors)
        bar2 = ax.bar(pos + width, act_reductions, width, color=act_colors)
        ax.set_xticks(pos + width / 2)
        ax.set_xticklabels(reduc_names, fontsize=10)
        # ax.set_yticks(np.arange(-100,105,5))
        # for label in ax.yaxis.get_ticklabels()[1::2]:
        #     label.set_visible(False)
        ax.set_ylabel('Percent Change in Muscle Activity', fontsize=12)
        ax.grid(which='both', axis='both', linestyle='--')
        ax.set_axisbelow(True)
        ax.legend([bar1, bar2], ['excitations', 'activations'])
        fig.tight_layout()
        fig.savefig(target[0])
        pl.close(fig)

def construct_multiindex_tuples(study, subjects, conditions, cycle_nums,
    muscle_level=False):
    ''' Construct multiindex tuples and list of cycles for DataFrame indexing.
    '''
    
    multiindex_tuples = list()
    cycles = list()

    for subject in study.subjects:
        if not subject.num in subjects: continue
        for cond_name in conditions:
            cond = subject.get_condition(cond_name)
            if not cond: continue
            # We know there is only one overground trial, but perhaps it
            # has not yet been added for this subject.
            assert len(cond.trials) <= 1
            if len(cond.trials) == 1:
                trial = cond.trials[0]
                for cycle in trial.cycles:
                    if cycle.name not in cycle_nums: 
                        continue
                    else:
                        cycles.append(cycle)
                        if not muscle_level:
                            multiindex_tuples.append((
                                cycle.subject.name,
                                # This must be the full ID, not just the cycle
                                # name, because 'cycle01' from subject 1 has
                                # nothing to do with 'cycle01' from subject 2
                                # (whereas the 'walk2' condition for subject 1 # isrelated to 'walk2' for subject 2).
                                cycle.name))
                        if muscle_level:
                            for mname in study.muscle_names:
                                multiindex_tuples.append((
                                    cycle.subject.name,
                                    cycle.name,
                                    mname))

    return multiindex_tuples, cycles

class TaskAggregateMetabolicRateSlackKinematics(osp.StudyTask):
    """Aggregate metabolic rate without and with mods across all subjects and
    gait cycles for each condition provided."""
    REGISTRY = []
    def __init__(self, study, mods=None, subjects=None, cycle_nums=None,
            assist_conds=['low','med','high','max']):
        super(TaskAggregateMetabolicRateSlackKinematics, self).__init__(study)
        
        print('\n!!!!!!!\nTaskAggregateMetabolicRateSlackKinematics unedited whole_fpath - edited assisted for both cost functions\n!!!!\n')
        
        self.name = 'aggregate_metabolic_rate_slack_kinematics'
        self.whole_fpath = os.path.join('metabolic_rate', 
            'whole_body_metabolic_rates_slack_kinematics.csv')
        self.muscs_fpath = os.path.join('metabolic_rate', 
            'muscle_metabolic_rates_slack_kinematics.csv')   
        self.doc = 'Aggregate metabolic rate.'
        self.study = study

        if cycle_nums:
            self.cycle_nums = cycle_nums
        else:
            self.cycle_nums = ['cycle%02i'% i for i in range(1, 11)]
    
        if subjects == None:
            subjects = [s.num for s in study.subjects]

        # Get multiindex tuples for DataFrame indexing for both whole body,
        # and muscle level metabolic rate. Also get cycles list.
        self.multiindex_tuples, cycles = construct_multiindex_tuples(study, 
            subjects, ['slack'], cycle_nums, muscle_level=False)

        self.mod_for_file_dep = list()
        deps = list()


        ##!!! -> do we need two cases of this for the default and the Met trials??
        # lets give it a shot!
        print ("\n_______________________________________________________________________________")
        print ("Reminder to check which files are being looked at: Met or Default funcitons.")
        print ("Look in tasks.py in the TaskAggregateMetabolicRateSlackKinematics task.")
        print ("_______________________________________________________________________________\n")
        # for cycle in cycles:
        #     if cycle.name in self.cycle_nums:
        #         self.mod_for_file_dep.append('experiment')
                
        # print self.mod_for_file_dep
        
        ##!!! TODO: not sure if any of this division is necessary, or useful
        # for using the 'Default' cost function
        # """
        # Prepare for processing simulations of experiments.
        for cycle in cycles:
            if cycle.name in self.cycle_nums:
                self.mod_for_file_dep.append('experiment')
                deps.append(os.path.join(
                        self.study.config['results_path'], 'experiments', 
                        cycle.subject.name, 'slack', 'mrs', cycle.name,
                        '%s_%s_slack_%s_mrs.mat' % (
                            study.name, cycle.subject.name, cycle.name))
                        )

        # Prepare for processing simulations of mods.
        for cond in assist_conds:
            for cycle in cycles:
                if cycle.name in self.cycle_nums:
                    self.mod_for_file_dep.append('mrsmod_%s' % cond)
                    deps.append(os.path.join(
                            self.study.config['results_path'],
                            'mrsmod_%s' % cond, cycle.subject.name, 'slack', # this use to be 'slack' 
                            'mrs',cycle.name, '%s_%s_slack_%s_mrs.mat' % # this use to have slack in it third
                                (study.name, cycle.subject.name, cycle.name))
                        )
        # """

        # For using the 'Met' cost function        
        """
        # Prepare for processing simulations of experiments.
        for cycle in cycles:
            if cycle.name in self.cycle_nums:
                self.mod_for_file_dep.append('experiment')
                deps.append(os.path.join(
                        self.study.config['results_path'], 'experiments', 
                        cycle.subject.name, 'slack', 'mrs', cycle.name,
                        'Met', '%s_%s_slack_%s_mrs.mat' % (
                            study.name, cycle.subject.name, cycle.name))
                        )

        # Prepare for processing simulations of mods.
        for cond in assist_conds:
            for cycle in cycles:
                if cycle.name in self.cycle_nums:
                    self.mod_for_file_dep.append('mrsmod_%s' % cond)
                    deps.append(os.path.join(
                            self.study.config['results_path'],
                            'mrsmod_%s' % cond, cycle.subject.name, 'slack', 
                            'mrs',cycle.name, 'Met', '%s_%s_slack_%s_mrs.mat' %
                                (study.name, cycle.subject.name, cycle.name))
                        )
        """
        
        self.add_action(deps,
                [os.path.join(study.config['analysis_path'], 
                self.whole_fpath)],
                self.aggregate_metabolic_rate)

    def aggregate_metabolic_rate(self, file_dep, target):
        import numpy as np
        from collections import OrderedDict
        metabolic_rate = OrderedDict()

        for ifile, fpath in enumerate(file_dep):
            df = util.hdf2pandas(fpath, 'DatStore/MetabolicRate/whole_body')
            this_mod = self.mod_for_file_dep[ifile]
            if not this_mod in metabolic_rate:
                metabolic_rate[this_mod] = list()
            metabolic_rate[this_mod].append(df[0][0])


        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced-hierarchical
        index = pd.MultiIndex.from_tuples(self.multiindex_tuples,
                names=['subject','cycle'])
        df = pd.DataFrame(metabolic_rate, index=index)

        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# columns contain whole body metabolic rate normalized by'
                    ' subject mass (W/kg)\n')
            df.to_csv(f)

    def aggregate_metabolic_rate_muscles(self, file_dep, target):
        import numpy as np
        from collections import OrderedDict
        metabolic_rate = OrderedDict()
        for ifile, fpath in enumerate(file_dep):
            df = util.hdf2pandas(fpath, 
                'DatStore/MetabolicRate/individual_muscles',
                labels=self.study.muscle_names)
            this_mod = self.mod_for_file_dep[ifile]
            if not this_mod in metabolic_rate:
                metabolic_rate[this_mod] = list()
            for muscle in self.study.muscle_names:
                metabolic_rate[this_mod].append(df[muscle][0])
       
        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced-hierarchical
        index = pd.MultiIndex.from_tuples(self.multiindex_tuples_musc,
                names=['subject', 'condition', 'cycle', 'muscle'])

        df = pd.DataFrame(metabolic_rate, index=index)

        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# columns contain muscle metabolic rates normalized by '
                    'subject mass (W/kg)\n')
            df.to_csv(f)



class TaskAggregatePlotStrideLengthsAssistedKinematics(osp.StudyTask):
    """Aggregate stride lengths for all subjects and
    gait cycles for each condition provided."""
    REGISTRY = []
    def __init__(self, study, mods=None, subjects=None, cycle_nums=None,
            assist_conds=['low','med','high','max']):
            # assist_conds=['max']):
        super(TaskAggregatePlotStrideLengthsAssistedKinematics, self).__init__(
            study)
        self.name = 'aggregate_plot_stride_lengths_assisted_kinematics'
        results_fpath = 'stride_lengths_assisted_kinematics.csv'
        plot_fpath = 'stride_lengths_assisted_kinematics.png'
        self.doc = 'Aggregate stride lengths.'
        self.study = study

        if cycle_nums:
            self.cycle_nums = cycle_nums
        else:
            self.cycle_nums = ['cycle%02i'% i for i in range(1, 11)]

        if subjects == None:
            subjects = [s.num for s in study.subjects]

        # Get multiindex tuples for DataFrame indexing for both whole body,
        # and muscle level metabolic rate. Also get cycles list.
        self.multiindex_tuples, cycles = construct_multiindex_tuples(study, 
            subjects, ['slack'], cycle_nums, muscle_level=False)

        # print self.multiindex_tuples
        # print cycles
        done_subjs = list()
        try:
            self.mod_for_file_dep = list()
            deps = list()
            # need to go through the subjects to get the lengths
            for cycle in cycles:
                if cycle.name in self.cycle_nums:
                    if cycle.subject.name not in done_subjs:
                        done_subjs.append(cycle.subject.name)
                        # print "appended one"
                        self.mod_for_file_dep.append('experiment')
                        deps.append(os.path.join(
                                self.study.config['results_path'], 'experiments', 
                                cycle.subject.name, 'slack', 'ik',
                                'avg_stride_lengths.csv'))

            del done_subjs
            done_subjs = list()
            # Prepare for processing simulations of mods.
            for cond in assist_conds:
                for cycle in cycles:
                    if cycle.name in self.cycle_nums:
                        if cycle.subject.name not in done_subjs:
                            done_subjs.append(cycle.subject.name)
                            self.mod_for_file_dep.append('mrsmod_%s' % cond)
                            deps.append(os.path.join(
                                    self.study.config['results_path'],
                                    'experiments', cycle.subject.name, cond, 
                                    'ik', 'avg_stride_lengths.csv'))
                del done_subjs
                done_subjs = list()
            # print deps
            self.add_action(deps,
                [os.path.join(study.config['analysis_path'], 
                results_fpath), os.path.join(study.config['analysis_path'],
                plot_fpath)],
                self.aggregate_stride_lengths)

        except:
            print('did not work for stride lengths')


        
    def aggregate_stride_lengths(self, file_dep, target):
        import numpy as np
        from collections import OrderedDict
        
        strides = OrderedDict()
        df = pd.DataFrame()
        for ifile, fpath in enumerate(file_dep):
            df_temp = pd.read_csv(fpath, index_col=[0,1], skiprows=1)
            df = df.append(df_temp)

        # print target
        with file(target[0], 'w') as f:
            f.write(' columns contain step lengths normalized by height\n')
            df.to_csv(f)

        ## TODO now for the averaging and plotting parts
        print ('\nplotting stuff')

        sns.set(style="whitegrid")
        # df.reset_index(inplace=True)
        # axtry = sns.boxplot(x=df["step_length"], order=['slack', 'max'])

        dfnew = df.reset_index(level='condition')
        dfnewnew = dfnew.reset_index(level='subject')
        print (dfnewnew.keys())
        axnew = sns.boxplot(x='condition',y='step_length',data=dfnewnew)
        fignew = axnew.get_figure()
        fignew.savefig(target[1])
        pl.close(fignew)



class TaskAggregateMetabolicRateAssistedKinematics(osp.StudyTask):
    """Aggregate metabolic rate without and with mods across all subjects and
    gait cycles for each condition provided."""
    REGISTRY = []
    def __init__(self, study, mods=None, subjects=None, cycle_nums=None,
            assist_conds=['max']):
        # assist_conds=['low','med','high','max']
        super(TaskAggregateMetabolicRateAssistedKinematics, self).__init__(
            study)
        self.name = 'aggregate_metabolic_rate_assisted_kinematics'
        
        # this is the path where the files will get saved
        # self.whole_fpath = os.path.join('metabolic_rate', 
            # 'whole_body_metabolic_rates_assisted_kinematics.csv')
        # attempt to make this into a tuple, one for default, one for Met
        default_fpath = 'metabolic_rate/default_whole_body_metabolic_rates_assisted_kinematics.csv'
        Met_fpath = 'metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        default_noassist_fpath = 'metabolic_rate/default_noassist_whole_body_metabolic_rates_assisted_kinematics.csv'
        Met_noassist_fpath = 'metabolic_rate/Met_noassist_whole_body_metabolic_rates_assisted_kinematics.csv'
        both_fpath = (default_fpath, Met_fpath, default_noassist_fpath, Met_noassist_fpath)
        self.whole_fpath = both_fpath
        ## TODO make sure that all the whole_fpath places refrerence a single entry

        self.muscs_fpath = os.path.join('metabolic_rate', 
            'muscle_metabolic_rates_assisted_kinematics.csv')   
        self.doc = 'Aggregate metabolic rate.'
        self.study = study


        if cycle_nums:
            self.cycle_nums = cycle_nums
        else:
            self.cycle_nums = ['cycle%02i'% i for i in range(1, 11)]
    
        if subjects == None:
            subjects = [s.num for s in study.subjects]

        # Get multiindex tuples for DataFrame indexing for both whole body,
        # and muscle level metabolic rate. Also get cycles list.
        self.multiindex_tuples, cycles = construct_multiindex_tuples(study, 
            subjects, ['slack'], cycle_nums, muscle_level=False)


        # Prepare for processing simulations of experiments.
        ##TODO alter to do the noassist cases as well

        print ("here to alter what all conditions etc get aggregated: line 1122/1511 ->")
        if study.costFunction == 'Met':
            # try the metabolic cost function with assist
            try:
                self.mod_for_file_dep = list()
                deps = list()
                # for the met function
                # Prepare for processing simulations of experiments.
                for cycle in cycles:
                    if cycle.name in self.cycle_nums:
                        self.mod_for_file_dep.append('experiment')
                        deps.append(os.path.join(
                                self.study.config['results_path'], 'experiments', 
                                cycle.subject.name, 'slack', 'mrs', cycle.name,
                                'Met', '%s_%s_slack_%s_mrs.mat' % (
                                    study.name, cycle.subject.name, cycle.name))
                                )                
                # Prepare for processing simulations of mods.
                for cond in assist_conds:
                    for cycle in cycles:
                        if cycle.name in self.cycle_nums:
                            self.mod_for_file_dep.append('mrsmod_%s' % cond)
                            deps.append(os.path.join(
                                    self.study.config['results_path'],
                                    'mrsmod_%s' % cond, cycle.subject.name, cond, 
                                    'mrs',cycle.name, 'Met', '%s_%s_%s_%s_mrs.mat' %
                                        (study.name, cycle.subject.name, cond, 
                                            cycle.name))
                                )

                self.add_action(deps,
                    [os.path.join(study.config['analysis_path'], 
                    self.whole_fpath[1])],
                    self.aggregate_metabolic_rate)
                # self.add_action(deps,
                #     [os.path.join(study.config['analysis_path'], 
                #     self.whole_fpath[1])],
                #     self.aggregate_metabolic_rate_muscles)

            except:
                print('did not work for the met function')

        else:    
            # for the default cost function with assist
            try:    
                self.mod_for_file_dep = list()
                deps = list()
                for cycle in cycles:
                    if cycle.name in self.cycle_nums:
                        self.mod_for_file_dep.append('experiment')
                        deps.append(os.path.join(
                                self.study.config['results_path'], 'experiments', 
                                cycle.subject.name, 'slack', 'mrs', cycle.name,
                                '%s_%s_slack_%s_mrs.mat' % (
                                    study.name, cycle.subject.name, cycle.name))
                                )
                # Prepare for processing simulations of mods.
                for cond in assist_conds:
                    for cycle in cycles:
                        if cycle.name in self.cycle_nums:
                            self.mod_for_file_dep.append('mrsmod_%s' % cond)
                            deps.append(os.path.join(
                                    self.study.config['results_path'],
                                    'mrsmod_%s' % cond, cycle.subject.name, cond, 
                                    'mrs',cycle.name, '%s_%s_%s_%s_mrs.mat' %
                                        (study.name, cycle.subject.name, cond, 
                                            cycle.name))
                                )

                self.add_action(deps,
                    [os.path.join(study.config['analysis_path'], 
                    self.whole_fpath[0])],
                    self.aggregate_metabolic_rate)
                # self.add_action(deps,
                #     [os.path.join(study.config['analysis_path'], 
                #     self.whole_fpath[1])],
                #     self.aggregate_metabolic_rate_muscles)

            except:
                print('did not work for the default')

        
        # ### adding for third loop all noassist on default
        # try:    
        #     self.mod_for_file_dep = list()
        #     deps = list()
        #     for cycle in cycles:
        #         if cycle.name in self.cycle_nums:
        #             self.mod_for_file_dep.append('experiment')
        #             deps.append(os.path.join(
        #                     self.study.config['results_path'], 'experiments', 
        #                     cycle.subject.name, 'slack', 'mrs', cycle.name,
        #                     '%s_%s_slack_%s_mrs.mat' % (
        #                         study.name, cycle.subject.name, cycle.name))
        #                     )
        #     # Prepare for processing simulations of mods.
        #     for cond in assist_conds:
        #         for cycle in cycles:
        #             if cycle.name in self.cycle_nums:
        #                 self.mod_for_file_dep.append('mrsmod_%s' % cond)
        #                 deps.append(os.path.join(
        #                         self.study.config['results_path'], 'experiments',
        #                         cycle.subject.name, cond, 'mrs', cycle.name, 
        #                         '%s_%s_%s_%s_mrs.mat' %
        #                             (study.name, cycle.subject.name, cond,
        #                                 cycle.name))
        #                     )

        #     self.add_action(deps,
        #         [os.path.join(study.config['analysis_path'], 
        #         self.whole_fpath[2])],
        #         self.aggregate_metabolic_rate)
        #     # self.add_action(deps,
        #     #     [os.path.join(study.config['analysis_path'], 
        #     #     self.whole_fpath[1])],
        #     #     self.aggregate_metabolic_rate_muscles)

        # except:
        #     print('did not work for the noassist default')


        # ### adding forth loop for all noassist on the Met
        # try:
        #     self.mod_for_file_dep = list()
        #     deps = list()
        #     # for the met function
        #     # Prepare for processing simulations of experiments.
        #     for cycle in cycles:
        #         if cycle.name in self.cycle_nums:
        #             self.mod_for_file_dep.append('experiment')
        #             deps.append(os.path.join(
        #                     self.study.config['results_path'], 'experiments', 
        #                     cycle.subject.name, 'slack', 'mrs', cycle.name,
        #                     'Met', '%s_%s_slack_%s_mrs.mat' % (
        #                         study.name, cycle.subject.name, cycle.name))
        #                     )
        #     for cond in assist_conds:
        #         for cycle in cycles:
        #             if cycle.name in self.cycle_nums:
        #                 self.mod_for_file_dep.append('mrsmod_%s' % cond)
        #                 deps.append(os.path.join(
        #                         self.study.config['results_path'], 'experiments',
        #                         cycle.subject.name, '%s' % cond, 'mrs', cycle.name,
        #                         'Met', '%s_%s_%s_%s_mrs.mat' % 
        #                         (study.name, cycle.subject.name, cond, cycle.name))
        #                     )
                        
        #     self.add_action(deps,
        #         [os.path.join(study.config['analysis_path'], 
        #         self.whole_fpath[3])],
        #         self.aggregate_metabolic_rate)        
        # except:
        #     print('didnt work on noassist met')

        # ## right now this gets called in each try loop
        # # self.add_action(deps,
        # #         [os.path.join(study.config['analysis_path'], 
        # #         self.whole_fpath[1])],
        # #         self.aggregate_metabolic_rate)


        
    def aggregate_metabolic_rate(self, file_dep, target):
        import numpy as np
        from collections import OrderedDict
        metabolic_rate = OrderedDict()
        

        for ifile, fpath in enumerate(file_dep):
            # print fpath
            df = util.hdf2pandas(fpath, 'DatStore/MetabolicRate/whole_body')
            this_mod = self.mod_for_file_dep[ifile]
            if not this_mod in metabolic_rate:
                metabolic_rate[this_mod] = list()
            metabolic_rate[this_mod].append(df[0][0])
       
        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced-hierarchical
        index = pd.MultiIndex.from_tuples(self.multiindex_tuples,
                names=['subject', 'cycle'])
        df = pd.DataFrame(metabolic_rate, index=index)

        # print "here to figure out what is going on"
        # print df
        # print df['subject024']['cycle01']
        # maybe add in a conditional statement here
        
        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# columns contain whole body metabolic rate normalized by'
                    ' subject mass (W/kg)\n')
            df.to_csv(f)


    def aggregate_metabolic_rate_muscles(self, file_dep, target):
        import numpy as np
        from collections import OrderedDict


        metabolic_rate = OrderedDict()
        for ifile, fpath in enumerate(file_dep):
            df = util.hdf2pandas(fpath, 
                'DatStore/MetabolicRate/individual_muscles',
                labels=self.study.muscle_names)
            this_mod = self.mod_for_file_dep[ifile]
            if not this_mod in metabolic_rate:
                metabolic_rate[this_mod] = list()
            for muscle in self.study.muscle_names:
                metabolic_rate[this_mod].append(df[muscle][0])
       
        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced-hierarchical
        index = pd.MultiIndex.from_tuples(self.multiindex_tuples_musc,
                names=['subject', 'cycle', 'muscle'])

        df = pd.DataFrame(metabolic_rate, index=index)

        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# columns contain muscle metabolic rates normalized by '
                    'subject mass (W/kg)\n')
            df.to_csv(f)


class TaskPlotMetabolicReductions(osp.StudyTask):
    """Plot reductions from simulations and compare to Quinlivan et al. 2017
    results.

    This will plot the reduction in metabolic cost based on the slack 
    kinematics with and without the assistance profiles. 
    """

    ### TODO edit the plotting function to handle both of the cost functions and
    #   some new plots as well
    REGISTRY = []
    def __init__(self, study):
        super(TaskPlotMetabolicReductions, self).__init__(study)
        self.name = 'plot_metabolic_reductions'

        # potentially going to have many different things in here doubled
        default_analysis_path_assist = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'default_whole_body_metabolic_rates_assisted_kinematics.csv')
        Met_analysis_path_assist = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'Met_whole_body_metabolic_rates_assisted_kinematics.csv')
        default_analysis_path_slack = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'default_whole_body_metabolic_rates_slack_kinematics.csv')
        Met_analysis_path_slack = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'Met_whole_body_metabolic_rates_slack_kinematics.csv')
        ### TODO add in the unassisted cases and plots for them throughout !!!
        default_analysis_path_noassist_assist = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'default_noassist_whole_body_metabolic_rates_assisted_kinematics.csv')
        Met_analysis_path_noassist_assist = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'Met_noassist_whole_body_metabolic_rates_assisted_kinematics.csv')
        ## adding in another slack kin one for the noassist cases - not used yet (1/22/19)
        default_noassist_analysis_path_slack = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'default_noassist_whole_body_metabolic_rates_slack_kinematics.csv')
        Met_noassist_analysis_path_slack = os.path.join(study.config['analysis_path'], 'metabolic_rate',
            'Met_noassist_whole_body_metabolic_rates_slack_kinematics.csv')

        # attempt to just tuple everything (default = 0, Met = 1)
        
        ### TODO finish the slack kinematics pipeline
        # self.slack_kin_results_fpath = os.path.join(
            # study.config['analysis_path'], 'metabolic_rate', 
            # 'whole_body_metabolic_rates_slack_kinematics.csv')


        # assisted kinematics pipe
        # self.assisted_kin_results_fpath = os.path.join(
            # study.config['analysis_path'], 'metabolic_rate', 
            # 'whole_body_metabolic_rates_assisted_kinematics.csv')
        self.assisted_kin_results_fpath = (default_analysis_path_assist, Met_analysis_path_assist)
        self.slack_kin_results_fpath = (default_analysis_path_slack, Met_analysis_path_slack)
        self.unassisted_assisted_kin_results_fpath = (default_analysis_path_noassist_assist, Met_analysis_path_noassist_assist)
        self.slack_kin_noassist_result_fpath = (default_noassist_analysis_path_slack, Met_noassist_analysis_path_slack)
        


        # Define the plots:
        ## individual subject plots divided by cost function
        self.default_metabolic_plot_individual_png = os.path.join(study.config['analysis_path'], 
            'metabolic_rate', 'default_metabolic_reductions_individual.png')
        self.default_metabolic_plot_individual_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_metabolic_reductions_individual.pdf')
        self.Met_metabolic_plot_individual_png = os.path.join(study.config['analysis_path'], 
            'metabolic_rate', 'Met_metabolic_reductions_individual.png')
        self.Met_metabolic_plot_individual_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_metabolic_reductions_individual.pdf')
        ### NEW
        self.default_noassist_metabolic_plot_individual_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_individual.png')
        self.default_noassist_metabolic_plot_individual_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_individual.pdf')
        self.Met_noassist_metabolic_plot_individual_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_individual.png')
        self.Met_noassist_metabolic_plot_individual_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_individual.pdf')

        # figure with only the population averages
        self.default_metabolic_plot_average_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_metabolic_reductions_average.png')
        self.default_metabolic_plot_average_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_metabolic_reductions_average.pdf')
        self.Met_metabolic_plot_average_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_metabolic_reductions_average.png')
        self.Met_metabolic_plot_average_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_metabolic_reductions_average.pdf')
        ### new
        self.default_noassist_metabolic_plot_average_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_average.png')
        self.default_noassist_metabolic_plot_average_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_average.pdf')
        self.Met_noassist_metabolic_plot_average_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_average.png')
        self.Met_noassist_metabolic_plot_average_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_average.pdf')

        # combine the population averages on a single plot
        self.both_metabolic_plot_average_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'both_metabolic_reductions_average.png')
        self.both_metabolic_plot_average_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'both_metabolic_reductions_average.pdf')
        ### going to try and splice the third one into the above as well


        # figures with the individuals and the population average all together
        self.default_metabolic_plot_combined_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_metabolic_reductions_combined.png')
        self.default_metabolic_plot_combined_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_metabolic_reductions_combined.pdf')
        self.Met_metabolic_plot_combined_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_metabolic_reductions_combined.png')
        self.Met_metabolic_plot_combined_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_metabolic_reductions_combined.pdf')
        ### adapted kinematics without assistance
        self.default_noassist_metabolic_plot_combined_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_combined.png')
        self.default_noassist_metabolic_plot_combined_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'default_noassist_metabolic_reductions_combined.pdf')
        self.Met_noassist_metabolic_plot_combined_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_combined.png')
        self.Met_noassist_metabolic_plot_combined_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_noassist_metabolic_reductions_combined.pdf')
        # raw figures - MET
        self.raw_metabolic_plot_combined_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_metabolic_plot_combined.png')
        self.raw_metabolic_plot_combined_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_metabolic_plot_combined.pdf')
        self.raw_metabolic_plot_avg_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_metabolic_plot_avg.png')
        self.raw_metabolic_plot_avg_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_metabolic_plot_avg.pdf')

        # add in this one for the making of figs that compare the different models and 
        # all the cost functions and things
        self.Met_reduction_models_avg_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_reduction_models_avg.png')
        self.Met_reduction_models_avg_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'Met_reduction_models_avg.pdf')
        self.raw_reduction_models_avg_png = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_reduction_models_avg.png')
        self.raw_reduction_models_avg_pdf = os.path.join(study.config['analysis_path'],
            'metabolic_rate', 'raw_reduction_models_avg.pdf')



        # raw figures - DEFAULT - act^2
        # self.raw_actsqr_plot_combined_png = os.path.join(study.config['analysis_path'],
        #     'metabolic_rate', 'raw_actsqr_plot_combined.png')
        # self.raw_actsqr_plot_combined_pdf = os.path.join(study.config['analysis_path'],
        #     'metabolic_rate', 'raw_actsqr_plot_combined.pdf')
        # self.raw_actsqr_plot_avg_png = os.path.join(study.config['analysis_path'],
        #     'metabolic_rate', 'raw_actsqr_plot_avg.png')
        # self.raw_actsqr_plot_avg_pdf = os.path.join(study.config['analysis_path'],
        #     'metabolic_rate', 'raw_actsqr_plot_avg.pdf')

        # so right now I have it so that all the cost functions and all the subjects will be plotted
        # on a single figure, but in the future could create a figure for all of the comparisons 
        # of functions vs experiment individually on their own graph and then do avg all together


        ## need to figure out how to add things in here
        # figure this out
        # self.add_action([self.slack_kin_results_fpath,
        #                  self.assisted_kin_results_fpath],
        #                 [self.metabolic_plot_png, 
        #                  self.metabolic_plot_pdf],
        #                  self.plot_metabolics)
        # self.add_action([self.slack_kin_results_fpath,
        #                  self.assisted_kin_results_fpath],
        #                 [self.metabolic_plot_individual_png, self.metabolic_plot_individual_pdf,
        #                 self.metabolic_plot_average_png, self.metabolic_plot_average_pdf,
        #                 self.metabolic_plot_combined_png, self.metabolic_plot_combined_pdf],
        #                  self.plot_metabolics)

        # the stuff in the slack_kin_results_fpath is only for the slack kin sims, so not used right now (1/22/19)


        # print('\nhere we go')
        # print(self.slack_kin_results_fpath)
        # print(self.assisted_kin_results_fpath)
        # print('\n')
        # print(self.unassisted_assisted_kin_results_fpath)
        # print('end\n')
        # best way could be to just do them one at a time by function
        # self.add_action([self.slack_kin_results_fpath[0],
        #                  self.assisted_kin_results_fpath[0]],
        #                 [self.default_metabolic_plot_individual_png, self.default_metabolic_plot_individual_pdf,
        #                 self.default_metabolic_plot_average_png, self.default_metabolic_plot_average_pdf,
        #                 self.default_metabolic_plot_combined_png, self.default_metabolic_plot_combined_pdf,
        #                 self.both_metabolic_plot_average_png, self.both_metabolic_plot_average_pdf],
        #                  self.plot_metabolics)


        # self.add_action([self.slack_kin_results_fpath[1],
        #                  self.assisted_kin_results_fpath[1]],
        #                 [self.Met_metabolic_plot_individual_png, self.Met_metabolic_plot_individual_pdf,
        #                 self.Met_metabolic_plot_average_png, self.Met_metabolic_plot_average_pdf,
        #                 self.Met_metabolic_plot_combined_png, self.Met_metabolic_plot_combined_pdf],
        #                  self.plot_metabolics)

        ## NEW attempt to create the actions for the noassist case
        # self.add_action([self.slack_kin_noassist_result_fpath[0],
        #                  self.unassisted_assisted_kin_results_fpath[0]],
        #                  [self.default_noassist_metabolic_plot_individual_png, self.default_noassist_metabolic_plot_individual_pdf,
        #                  self.default_noassist_metabolic_plot_average_png, self.default_noassist_metabolic_plot_average_pdf,
        #                  self.default_noassist_metabolic_plot_combined_png, self.default_noassist_metabolic_plot_combined_pdf],
        #                  self.plot_metabolics)

        # self.add_action([self.slack_kin_noassist_result_fpath[1],
        #                  self.unassisted_assisted_kin_results_fpath[1]],
        #                  [self.Met_noassist_metabolic_plot_individual_png, self.Met_noassist_metabolic_plot_individual_pdf,
        #                  self.Met_noassist_metabolic_plot_average_png, self.Met_noassist_metabolic_plot_average_pdf,
        #                  self.Met_noassist_metabolic_plot_combined_png, self.Met_noassist_metabolic_plot_combined_pdf],
        #                  self.plot_metabolics)

        ### NEW - make one for plotting the raw data one only
        # self.add_action([self.slack_kin_results_fpath[1],
        #                  self.assisted_kin_results_fpath[1]],
        #                 [self.raw_metabolic_plot_avg_png, self.raw_metabolic_plot_avg_pdf,
        #                 self.raw_metabolic_plot_combined_png, self.raw_metabolic_plot_combined_pdf],
        #                  self.plot_metabolics)


        self.add_action([self.slack_kin_results_fpath[1],
                        self.assisted_kin_results_fpath[1]],
                        [self.Met_reduction_models_avg_png, self.Met_reduction_models_avg_pdf,
                        self.raw_reduction_models_avg_png, self.raw_reduction_models_avg_pdf],
                        self.plot_metabolics)


    # TODO need to get rid of subject 126 from the raw experimental averages just for consistency 
    def plot_metabolics(self, file_dep, target):

        def get_metabolic_reductions(csv_file, mods, subject=None):
            # Process metabolic rate
            # -----------------------
            # The first three "columns" form a MultiIndex.
            # Skip the first line, which has comments.
            df = pd.read_csv(csv_file, index_col=[0, 1], skiprows=1)

            # Using walk2 condition to start. TODO: update if adding 
            # more conditions.
            # "xs" stands for "cross section"
            # df_walk2 = df.xs('walk2', level='condition')

            # Subtract the no-assist cost from all other columns.
            df_relchange = df.subtract(df['experiment'],
                    axis='index').divide(df['experiment'], axis='index')

            # Delete the 'experiments' column, whicih is no longer needed.
            df_relchange.drop('experiment', axis='columns', inplace=True)

            # if a subject is passed, it will go and get that specifically
            if subject:
                df_subj = df_relchange.xs(subject, level='subject')
            # otherwise, it takes the whole population
            else:
                df_subj = df_relchange.groupby(level='subject').mean()
            
            met_relchange_pcent_mean = df_subj.mean()[mods] * 100
            met_relchange_pcent_std = df_subj.std()[mods] * 100

            return met_relchange_pcent_mean, met_relchange_pcent_std
        
        def get_raw_metabolics(csv_file, mods, subject=None):
            # Process metabolic rate
            # -----------------------
            # The first three "columns" form a MultiIndex.
            # Skip the first line, which has comments.
            df = pd.read_csv(csv_file, index_col=[0, 1], skiprows=1)

            # # Subtract the no-assist cost from all other columns.
            # df_relchange = df.subtract(df['experiment'],
            #         axis='index').divide(df['experiment'], axis='index')

            # # Delete the 'experiments' column, whicih is no longer needed.
            # df_relchange.drop('experiment', axis='columns', inplace=True)

            # if a subject is passed, it will go and get that specifically
            if subject:
                df_subj = df.xs(subject, level='subject')
            # otherwise, it takes the whole population
            else:
                df_subj = df.groupby(level='subject').mean()
            
            met_relchange_pcent_mean = df_subj.mean()[mods]
            met_relchange_pcent_std = df_subj.std()[mods]

            return met_relchange_pcent_mean, met_relchange_pcent_std
        

        print ("\n#####################################################################")
        print ("\nokay here we are, time to get cracking")


        #############################################################################
        # defining figures
        #############################################################################
        # Plot changes in metabolic rate
        # figure for individual subjects
        fig_ind = pl.figure(figsize=(6, 6))
        ax_ind = fig_ind.add_subplot(1, 1, 1)
        # figure for population avg
        fig_avg = pl.figure(figsize=(6, 6))
        ax_avg = fig_avg.add_subplot(1, 1, 1)
        # figure for combined individual and population average overlaid
        fig_comb = pl.figure(figsize=(6, 6))
        ax_comb = fig_comb.add_subplot(1, 1, 1)
        # figure for the combined averages
        both_fig_avg = pl.figure(figsize=(6, 4))
        both_ax_avg = both_fig_avg.add_subplot(1,1,1)
        # plotting the raw metabolics values - averages
        rawboth_fig_avg = pl.figure(figsize=(6,4))
        rawboth_ax_avg = rawboth_fig_avg.add_subplot(1,1,1)
        # plotting the raw metabolics values - subjects averages
        rawboth_fig_comb = pl.figure(figsize=(6,4))
        rawboth_ax_comb = rawboth_fig_comb.add_subplot(1,1,1)

        # figures with model comparisons
        # % reductions
        models_fig_avg = pl.figure(figsize=(8,5))
        models_ax_avg = models_fig_avg.add_subplot(1,1,1)
        # raw metabolic values
        rawmodels_fig_avg = pl.figure(figsize=(8,5))
        rawmodels_ax_avg = rawmodels_fig_avg.add_subplot(1,1,1)


        # axis line characteristics
        ax_ind.axhline(linewidth=2, color='k')
        ax_avg.axhline(linewidth=2, color='k')
        ax_comb.axhline(linewidth=2, color='k')
        both_ax_avg.axhline(linewidth=2, color='k')
        rawboth_ax_avg.axhline(linewidth=2, color='k')
        rawboth_ax_comb.axhline(linewidth=2, color='k')
        models_ax_avg.axhline(linewidth=2, color='k')
        rawmodels_ax_avg.axhline(linewidth=2, color='k')


        ## !! going to change the mods to handle only the max and slack cases for the various model analyses
        # this is defining the x-axis for the plot as the different conditions
        mods = ['mrsmod_low', 'mrsmod_med', 'mrsmod_high', 'mrsmod_max']
        # mods = ['mrsmod_max']
        x_pos = np.arange(len(mods))

        # make seperate for he raw plot
        raw_mods = ['experiment', 'mrsmod_low', 'mrsmod_med', 'mrsmod_high', 'mrsmod_max']
        # raw_mods = ['experiment', 'mrsmod_max']
        raw_x_pos = np.arange(len(raw_mods))
        # note do not plot the % reductions for the raw values 

        # these are the experimental values
        experiment_mean = [-3.59, -6.45, -14.79, -22.83]
        # colors to plot the bars
        colors = [(0.662, 0.047, 0.913), (0.094, 0.756, 0.905), 
                  (0.094, 0.905, 0.266), (0.964, 0.607, 0.074)]
        # plotting the bars from the experiment on each figure
        ax_ind.bar(x_pos, experiment_mean, color=colors, alpha=0.5, label='experimental')
        ax_avg.bar(x_pos, experiment_mean, color=colors, alpha=0.5, label='experimental')
        ax_comb.bar(x_pos, experiment_mean, color=colors, alpha=0.5, label='experimental')
        both_ax_avg.bar(x_pos, experiment_mean, color=colors, alpha=0.5, label='experimental')
        # models_ax_avg.bar(x_pos, experiment_mean, color=colors, alpha=0.5, label='experimental')


        #############################################################################
        # individual % change plot
        #############################################################################
        # get the data for the individual subjects - based on which file-dep was input
        met_mean_assisted24, met_std_assisted24 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject024')
        met_mean_assisted77, met_std_assisted77 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject077')
        met_mean_assisted88, met_std_assisted88 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject088')
        met_mean_assisted112, met_std_assisted112 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject112')
        met_mean_assisted127, met_std_assisted127 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject127')
        met_mean_assisted128, met_std_assisted128 = get_metabolic_reductions(
            file_dep[1], mods, subject='subject128')

        # plot the data for the individual plots
        ax_ind.plot(x_pos, met_mean_assisted24, 
            color='red', marker='^', linestyle='--')
        ax_ind.plot(x_pos, met_mean_assisted77, 
            color='orange', marker='^', linestyle='--')
        ax_ind.plot(x_pos, met_mean_assisted88, 
            color='yellow', marker='^', linestyle='--')
        ax_ind.plot(x_pos, met_mean_assisted112, 
            color='green', marker='^', linestyle='--')
        ax_ind.plot(x_pos, met_mean_assisted127, 
            color='blue', marker='^', linestyle='--')
        ax_ind.plot(x_pos, met_mean_assisted128, 
            color='purple', marker='^', linestyle='--')        


        #############################################################################
        # defining all extra file paths
        #############################################################################
        default_assist = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/analysis/metabolic_rate/default_whole_body_metabolic_rates_assisted_kinematics.csv'
        Met_assist = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        
        default_noassist = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/analysis/metabolic_rate/default_noassist_whole_body_metabolic_rates_assisted_kinematics.csv'
        Met_noassist = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/analysis/metabolic_rate/Met_noassist_whole_body_metabolic_rates_assisted_kinematics.csv'
        
        assist_general_kin = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/results/revamp_muscle_metabolic_rates_Quinlivan2017.csv'
        raw_met = 'C:/Users/JP/code/repos/Stanford/delplab/data/AnkleHipExosuit/Metabolic_results_forplot.csv'
        raw_met_nostanding = 'C:/Users/JP/code/repos/Stanford/delplab/data/AnkleHipExosuit/Metabolic_results_forplot_withoutstanding.csv'

        # the multiple models stuff

        ## simple model first iteration
        # a^2
        simplemodel_a2 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_03182019/analysis/metabolic_rate/default_whole_body_metabolic_rates_assisted_kinematics.csv'
        # met
        simplemodel_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_03182019/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## simple model added weight
        # a^2
        simple_addedweight_a2 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_04092019_addedmass/analysis/metabolic_rate/default_whole_body_metabolic_rates_assisted_kinematics.csv'
        # met
        simple_addedweight_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_04092019_addedmass/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## full model
        # met
        full_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_05222019_fullmuscle/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## simple model calibrated
        # met
        simple_cal_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_06202019_simple_calibrated/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## simple model calibrated added weight
        # met 
        simple_cal_addedweight_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_06242019_simple_calibrated_addedmass/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## full muscle calibrated
        # met
        full_cal_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_07032019_fullmuscle_calibrated/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## full muscle calibrated 5 dof
        # met
        full_cal_5dof_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_07102019_fullmuscle_calibrated_5dof/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## full model calibrated with new weighting fixed 5 dof
        # met
        full_cal_new_5dof_met = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_08122019_newweighted_fullmusclemodel_calibrated_5dof/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## simple model cal new weighted minetti implicit
        # met
        simple_cal_new_minetti = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_09032019_simplemodelcalbrated_newweight_minettiimplicit/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        ## full model cal new weighted minetti implicit 5 dof
        # met
        full_cal_new_minetti_5dof = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_09112019_fullmuscle_calibrated_newweight_minettiimplicit_5dof/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'


        ## uchida models
        full_cal_5dof_uchida_25 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_09282019_fullmuscle_calibrated_5dof_uchida/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'
        full_cal_5dof_uchida_30 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit_09282019_fullmuscle_calibrated_5dof_uchida_sigmapt3/analysis/metabolic_rate/Met_whole_body_metabolic_rates_assisted_kinematics.csv'



        #############################################################################
        # gathereing raw met values
        #############################################################################
        # getting the data for the raw figure - metabolic cost function in file_dep[1] 
        raw_met_mean_assisted24, raw_met_std_assisted24 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject024')
        raw_met_mean_assisted77, raw_met_std_assisted77 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject077')
        raw_met_mean_assisted88, raw_met_std_assisted88 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject088')
        raw_met_mean_assisted112, raw_met_std_assisted112 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject112')
        raw_met_mean_assisted127, raw_met_std_assisted127 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject127')
        raw_met_mean_assisted128, raw_met_std_assisted128 = get_raw_metabolics(
            file_dep[1], raw_mods, subject='subject128')
        # this should maybe change to Met_assist

        # # activations squared cost functions specifically
        # raw_default_mean_assisted24, raw_default_std_assisted24 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject024')
        # raw_default_mean_assisted77, raw_default_std_assisted77 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject077')
        # raw_default_mean_assisted88, raw_default_std_assisted88 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject088')
        # raw_default_mean_assisted112, raw_default_std_assisted112 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject112')
        # raw_default_mean_assisted127, raw_default_std_assisted127 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject127')
        # raw_default_mean_assisted128, raw_default_std_assisted128 = get_raw_metabolics(
        #     default_assist, raw_mods, subject='subject128')

        # # activations squared cost functions specifically - nonassisted
        # raw_default_mean_noassisted24, raw_default_std_noassisted24 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject024')
        # raw_default_mean_noassisted77, raw_default_std_noassisted77 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject077')
        # raw_default_mean_noassisted88, raw_default_std_noassisted88 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject088')
        # raw_default_mean_noassisted112, raw_default_std_noassisted112 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject112')
        # raw_default_mean_noassisted127, raw_default_std_noassisted127 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject127')
        # raw_default_mean_noassisted128, raw_default_std_noassisted128 = get_raw_metabolics(
        #     default_noassist, raw_mods, subject='subject128')

        #  # metabolic cost functions specifically - nonassisted
        #  # activations squared cost functions specifically - nonassisted
        # raw_met_mean_noassisted24, raw_met_std_noassisted24 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject024')
        # raw_met_mean_noassisted77, raw_met_std_noassisted77 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject077')
        # raw_met_mean_noassisted88, raw_met_std_noassisted88 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject088')
        # raw_met_mean_noassisted112, raw_met_std_noassisted112 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject112')
        # raw_met_mean_noassisted127, raw_met_std_noassisted127 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject127')
        # raw_met_mean_noassisted128, raw_met_std_noassisted128 = get_raw_metabolics(
        #     Met_noassist, raw_mods, subject='subject128')

        #############################################################################
        # population averages - % change
        #############################################################################
        met_mean_assisted, met_std_assisted = get_metabolic_reductions(
                file_dep[1], mods)

        # default_met_mean_assisted, default_met_std_assisted = get_metabolic_reductions(
        #         default_assist, mods)
        Met_met_mean_assisted, Met_met_std_assisted = get_metabolic_reductions(
                Met_assist, mods)
        # default_noassist_met_mean_assisted, default_noassist_met_std_assisted = get_metabolic_reductions(
        #         default_noassist, mods)
        # Met_noassist_met_mean_assisted, Met_noassist_met_std_assisted = get_metabolic_reductions(
        #         Met_noassist, mods)
        
        # as a note some of these are redundant I think based on what is passed in file_dep        

        ## here I am going to get the population avg for all the models and functions
        simplemodel_a2_mean_assisted, simplemodel_a2_std_assisted = get_metabolic_reductions(
                simplemodel_a2, mods)
        simplemodel_met_mean_assisted, simplemodel_met_std_assisted = get_metabolic_reductions(
                simplemodel_met, mods)
        simple_addedweight_a2_mean_assisted, simple_addedweight_a2_std_assisted = get_metabolic_reductions(
                simple_addedweight_a2, mods)
        simple_addedweight_met_mean_assisted, simple_addedweight_met_std_assisted = get_metabolic_reductions(
                simple_addedweight_met, mods)
        full_met_mean_assisted, full_met_std_assisted = get_metabolic_reductions(
                full_met, mods)
        simple_cal_met_mean_assisted, simple_cal_met_std_assisted = get_metabolic_reductions(
                simple_cal_met, mods)
        simple_cal_addedweight_met_mean_assisted, simple_cal_addedweight_met_std_assisted = get_metabolic_reductions(
                simple_cal_addedweight_met, mods)
        full_cal_met_mean_assisted, full_cal_met_std_assisted = get_metabolic_reductions(
                full_cal_met, mods)
        full_cal_5dof_met_mean_assisted, full_cal_5dof_met_std_assisted = get_metabolic_reductions(
                full_cal_5dof_met, mods)
        full_cal_new_5dof_met_mean_assisted, full_cal_new_5dof_met_std_assisted = get_metabolic_reductions(
                full_cal_new_5dof_met, mods)
        simple_cal_new_minetti_mean_assisted, simple_cal_new_minetti_std_assisted = get_metabolic_reductions(
                simple_cal_new_minetti, mods)
        full_cal_new_minetti_5dof_mean_assisted, full_cal_new_minetti_5dof_std_assisted = get_metabolic_reductions(
                full_cal_new_minetti_5dof, mods)


        full_cal_5dof_uchida_25_mean_assisted, full_cal_5dof_uchida_25_std_assisted = get_metabolic_reductions(
                full_cal_5dof_uchida_25, mods)
        full_cal_5dof_uchida_30_mean_assisted, full_cal_5dof_uchida_30_std_assisted = get_metabolic_reductions(
                full_cal_5dof_uchida_30, mods)


        # plot the population average figure
        ax_avg.errorbar(x_pos, met_mean_assisted, 
            yerr=met_std_assisted,
            color='blue', marker='o', linestyle='--',
            label='assisted kinematics')


        # include numbers for the general walking kinematics (from nick and chris)
        # general_kinematics_mean_assisted = [-4.9031786, -7.313282274, -9.615899673, -11.65958781]
        df_gen = pd.read_csv(assist_general_kin, skiprows=1) # index_col=[0, 1], skiprows=1)
        df_gen2 = df_gen.copy(deep=True)
        df_gen2.drop(['condition','cycle','muscle','experiment','Quinlivan2017_opt_peak','Quinlivan2017_force_level_1',
                                'Quinlivan2017_force_level_2','Quinlivan2017_force_level_3','Quinlivan2017_force_level_4','Quinlivan2017_force_level_5',
                                'Quinlivan2017_force_level_6','Quinlivan2017_force_level_7','Quinlivan2017_force_level_8','Quinlivan2017_force_level_9',
                                'Quinlivan2017_force_level_10','experiment - cycle average','experiment - total average','low - cycleavg','low - totalavg',
                                'med - cycleavg','med - totalavg','high - cycleavg','high - totalavg','max - cycleavg','max - totalavg'],
                                axis='columns', inplace=True)
        df_gen_raw = df_gen2.copy(deep=True)
        df_gen2['low - musclesum'] = (df_gen2['low - musclesum'] - df_gen2['experiment - muscle sum']) / df_gen2['experiment - muscle sum']
        df_gen2['med - musclesum'] = (df_gen2['med - musclesum'] - df_gen2['experiment - muscle sum']) / df_gen2['experiment - muscle sum']        
        df_gen2['high - musclesum'] = (df_gen2['high - musclesum'] - df_gen2['experiment - muscle sum']) / df_gen2['experiment - muscle sum']
        df_gen2['max - musclesum'] = (df_gen2['max - musclesum'] - df_gen2['experiment - muscle sum']) / df_gen2['experiment - muscle sum']
        df_gen2.drop('experiment - muscle sum', axis='columns', inplace=True)
        # print(df_gen2)
        # df_gen3 = df_gen2.subtract(df_gen2['experiment - muscle sum'],
        #             axis='index').divide(df_gen2['experiment - muscle sum'], axis='index')
        # print(df_gen3)
        # df_gen2.drop('experiment - muscle sum', axis='columns', inplace=True)
        # print(df_gen2)
        # if subject:
        #     df_subj = df_relchange.xs(subject, level='subject')
        #     # otherwise, it takes the whole population
        # else:
        df_gen3 = df_gen2.groupby('subject').mean()
        general_kin_met_pcnt_mean = df_gen3.mean(axis='rows') * 100 #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100
        general_kin_met_pcnt_std = df_gen3.std(axis='rows') * 100   #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100


        ## manipulate the raw metabolic dataframe from Nick's study
        df_gen_raw2 = df_gen_raw.groupby('subject').mean()
        raw_genkin_mean = df_gen_raw2.mean(axis='rows')
        raw_genkin_std = df_gen_raw2.std(axis='rows')


        # add to the population average plot
        ax_avg.errorbar(x_pos, general_kin_met_pcnt_mean,
            yerr=general_kin_met_pcnt_std, color='black', marker='o',
            linestyle='--', label='Bianco et al. 2018')
        ax_avg.legend(fontsize=14)

        # plot the both function avgs
        # both_ax_avg.errorbar(x_pos, default_met_mean_assisted,
        #     yerr=default_met_std_assisted,
        #     color='red', marker='o', linestyle='--', label='$a^2$')
        # both_ax_avg.errorbar(x_pos, default_noassist_met_mean_assisted,
        #     yerr=default_noassist_met_std_assisted,
        #     color='tomato', marker='o', linestyle='--', label='$a^2$ no assistance')
        both_ax_avg.errorbar(x_pos, Met_met_mean_assisted,
            yerr=Met_met_std_assisted,
            color='blue', marker='o', linestyle='--', label='metabolic cost')
        # both_ax_avg.errorbar(x_pos, Met_noassist_met_mean_assisted,
        #     yerr=Met_noassist_met_std_assisted,
        #     color='deepskyblue', marker='o', linestyle='--', label='metabolic cost no assistance')
        both_ax_avg.errorbar(x_pos, general_kin_met_pcnt_mean,
            yerr=general_kin_met_pcnt_std,
            color='black', marker='^', linestyle='--', label='Bianco et al. 2018')
        both_ax_avg.errorbar
        both_ax_avg.legend(fontsize=14)




        ## now to make the figure with all the models
        qualspos = x_pos[:-1]
        pctplot = ['Experimental','Umberger 2010','Uchida 2016']
        pctmean = [-22.83, full_cal_5dof_met_mean_assisted['mrsmod_max'], full_cal_5dof_uchida_25_mean_assisted['mrsmod_max']]
        pctstd = [3.17, full_cal_5dof_met_std_assisted['mrsmod_max'], full_cal_5dof_uchida_25_std_assisted['mrsmod_max']]

        # import time
        # time.sleep(5)
        models_ax_avg.bar(qualspos, pctmean, yerr=pctstd, align='center', color=['darkslategrey','forestgreen','red'], ecolor='black', capsize=10, alpha=0.8)
        # models_ax_avg.set_ylabel('Percent Metabolic Reduction [%]')
        # models_ax_avg.set_xticks(qualspos)
        # models_ax_avg.set_xticklabels(pctplot)
        # models_ax_avg.grid(True)



        # models_ax_avg.plot(x_pos, general_kin_met_pcnt_mean,
        #     color ='black', marker='o', linestyle='-', label='Bianco 2017')
        # models_ax_avg.plot(x_pos, simplemodel_a2_mean_assisted,
        #     color='red', marker='o', linestyle=':', label='simple a2')
        # models_ax_avg.plot(x_pos, simple_addedweight_a2_mean_assisted,
        #     color='pink', marker='o', linestyle=':', label='simple add weight a2')
        # models_ax_avg.plot(x_pos, simplemodel_met_mean_assisted,
        #     color='mediumblue', marker='o', linestyle=':', label='simple met')
        # models_ax_avg.plot(x_pos, simple_addedweight_met_mean_assisted,
        #     color='cornflowerblue', marker='o', linestyle=':', label='simple add weight met')
        # models_ax_avg.plot(x_pos, simple_cal_met_mean_assisted,
        #     color='forestgreen', marker='o', linestyle=':', label='simple cal. met')
        # models_ax_avg.plot(x_pos, simple_cal_addedweight_met_mean_assisted,
        #     color='limegreen', marker='o', linestyle=':', label='simple cal. add weight met')
        # models_ax_avg.plot(x_pos, full_met_mean_assisted,
        #     color='violet', marker='^', linestyle=':', label='full met 3dof')
        # models_ax_avg.plot(x_pos, full_cal_met_mean_assisted,
        #     color='blueviolet', marker='^', linestyle=':', label='full cal. met 3dof')

        # models_ax_avg.plot(x_pos, full_cal_5dof_met_mean_assisted,
        #     color='forestgreen', marker='o', linestyle=':', label='Umberger 2010')
        
        # models_ax_avg.plot(x_pos, full_cal_new_5dof_met_mean_assisted,
        #     color='magenta', marker='^', linestyle=':', label='full cal. fixed 5dof')
        # models_ax_avg.plot(x_pos, simple_cal_new_minetti_mean_assisted,
        #     color='lightseagreen', marker='o', linestyle=':', label='simple cal. fixed minetti')
        # models_ax_avg.plot(x_pos, full_cal_new_minetti_5dof_mean_assisted,
        #     color='teal', marker='^', linestyle=':', label='full cal. fixed minetti 5dof')

        # models_ax_avg.plot(x_pos, full_cal_5dof_uchida_25_mean_assisted,
        #     color='red', marker='o', linestyle=':', label='Uchida 2016')
        # models_ax_avg.plot(x_pos, full_cal_5dof_uchida_30_mean_assisted,
        #     color='mediumblue', marker='^', linestyle=':', label='Uchida 2016, sigma=.30 MPa')
                
        # models_ax_avg.legend(fontsize=14)
        # models_ax_avg.legend(loc='upper left', bbox_to_anchor= (1, 0.9), ncol=1, 
        #     borderaxespad=0, frameon=False, fontsize=12)


        #############################################################################
        # population averages - raw cost
        #############################################################################
        # simulation - (called in the metabolic cost function as file_dep[1])
        raw_met_mean_assisted, raw_met_std_assisted = get_raw_metabolics(
                Met_assist, raw_mods)


        # gathering and handling the experimental indirect calorimetry values -> here we are
        df_raw = pd.read_csv(raw_met, index_col = [0], skiprows=1)
        # note we drop BW, height, standing Met rate
        df_raw.drop(['BW[kg]','Height[m]','Standing'], axis='columns', inplace=True)
        df_raw_slack = df_raw.copy(deep=True)
        df_raw_slack.drop(['LOW','MED','HIGH','MAX'], axis='columns', inplace=True)
        df_raw_slack['slackavg'] = df_raw_slack.mean(axis=1)
        df_raw_slack.drop(['Slack #1', 'Slack #2'], axis='columns', inplace=True)
        df_raw_slack['low'] = df_raw['LOW']
        df_raw_slack['med'] = df_raw['MED']
        df_raw_slack['high'] = df_raw['HIGH']
        df_raw_slack['max'] = df_raw['MAX']
        # now have a dataframe with the raw values from experimental
        # each entry is the average cost for one subject, below gives population average
        raw_met_mean = df_raw_slack.mean(axis='rows') #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100
        raw_met_std = df_raw_slack.std(axis='rows')   #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100

        # get the indirect calorimetry with standing subtracted
        df_raw_nostand = pd.read_csv(raw_met_nostanding, index_col = [0], skiprows=1)
        # note we drop BW, height, standing Met rate
        df_raw_nostand.drop(['BW[kg]','Height[m]','Standing'], axis='columns', inplace=True)
        df_raw_nostand_slack = df_raw_nostand.copy(deep=True)
        df_raw_nostand_slack.drop(['LOW','MED','HIGH','MAX'], axis='columns', inplace=True)
        df_raw_nostand_slack['slackavg'] = df_raw_nostand_slack.mean(axis=1)
        df_raw_nostand_slack.drop(['Slack #1', 'Slack #2'], axis='columns', inplace=True)
        df_raw_nostand_slack['low'] = df_raw_nostand['LOW']
        df_raw_nostand_slack['med'] = df_raw_nostand['MED']
        df_raw_nostand_slack['high'] = df_raw_nostand['HIGH']
        df_raw_nostand_slack['max'] = df_raw_nostand['MAX']
        # now have a dataframe with the raw values from experimental
        # each entry is the average cost for one subject, below gives population average
        raw_nostand_met_mean = df_raw_nostand_slack.mean(axis='rows') #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100
        raw_nostand_met_std = df_raw_nostand_slack.std(axis='rows')   #['low - musclesum','med - musclesum','high - musclesum','max - musclesum'] * 100


        # gets the population averages for the activations squared cost function
        # raw_default_mean_assisted, raw_default_std_assisted = get_raw_metabolics(
        #         default_assist, raw_mods)

        # gets the population averages for the activations squared cost function - no assistance
        # raw_default_mean_noassisted, raw_default_std_noassisted = get_raw_metabolics(
        #         default_noassist, raw_mods)

        # gets the population averages for met cost function - no assistance
        # raw_met_mean_noassisted, raw_met_std_noassisted = get_raw_metabolics(
        #         Met_noassist, raw_mods)


        # 2 * for all the sims, since we only have one leg
        # actual plotting for the raw cases 
        rawboth_ax_avg.errorbar(raw_x_pos, raw_met_mean,
            yerr=raw_met_std, color='black', marker='o',
            linestyle='--', label='indirect calorimetry')
        rawboth_ax_avg.errorbar(raw_x_pos, 2*raw_met_mean_assisted,
            yerr=raw_met_std_assisted, color='blue', marker='o',
            linestyle='--', label='metabolic cost sim')

        rawboth_ax_avg.errorbar(raw_x_pos, raw_nostand_met_mean,
            yerr=raw_nostand_met_std, color='grey', marker='o',
            linestyle='--', label='indirect calorimetry minus standing')


        # rawboth_ax_avg.errorbar(raw_x_pos, 2*raw_met_mean_noassisted,
        #     yerr=raw_default_std_assisted, color='deepskyblue', marker='o',
        #     linestyle='--', label='metabolic cost sim no assist')
        # rawboth_ax_avg.errorbar(raw_x_pos, 2*raw_default_mean_assisted,
        #     yerr=raw_default_std_assisted, color='red', marker='o',
        #     linestyle='--', label='activations squared cost sim')
        # rawboth_ax_avg.errorbar(raw_x_pos, 2*raw_default_mean_noassisted,
        #     yerr=raw_default_std_noassisted, color='tomato', marker='o',
        #     linestyle='--', label='activations squared cost sim no assist')

        rawboth_ax_avg.legend(fontsize=14)


        ## here for the full models picture
        raw_simplemodel_a2_mean_assisted, raw_simplemodel_a2_std_assisted = get_raw_metabolics(
                simplemodel_a2, raw_mods)

        raw_simplemodel_met_mean_assisted, raw_simplemodel_met_std_assisted = get_raw_metabolics(
                simplemodel_met, raw_mods)

        raw_simple_addedweight_a2_mean_assisted, raw_simple_addedweight_a2_std_assisted = get_raw_metabolics(
                simple_addedweight_a2, raw_mods)

        raw_simple_addedweight_met_mean_assisted, raw_simple_addedweight_met_std_assisted = get_raw_metabolics(
                simple_addedweight_met, raw_mods)

        raw_full_met_mean_assisted, raw_full_met_std_assisted = get_raw_metabolics(
                full_met, raw_mods)
        raw_full_met_mean_assisted_2 = raw_full_met_mean_assisted.dropna()

        raw_simple_cal_met_mean_assisted, raw_simple_cal_met_std_assisted = get_raw_metabolics(
                simple_cal_met, raw_mods)

        raw_simple_cal_addedweight_met_mean_assisted, raw_simple_cal_addedweight_met_std_assisted = get_raw_metabolics(
                simple_cal_addedweight_met, raw_mods)
        raw_simple_cal_addedweight_met_mean_assisted_2 = raw_simple_cal_addedweight_met_mean_assisted.dropna()

        raw_full_cal_met_mean_assisted, raw_full_cal_met_std_assisted = get_raw_metabolics(
                full_cal_met, raw_mods)
        raw_full_cal_met_mean_assisted_2 = raw_full_cal_met_mean_assisted.dropna()


        raw_full_cal_5dof_met_mean_assisted, raw_full_cal_5dof_met_std_assisted = get_raw_metabolics(
                full_cal_5dof_met, raw_mods)
        raw_full_cal_5dof_met_mean_assisted_2 = raw_full_cal_5dof_met_mean_assisted.dropna()
        raw_full_cal_5dof_met_std_assisted_2 = raw_full_cal_5dof_met_std_assisted.dropna()


        raw_full_cal_new_5dof_met_mean_assisted, raw_full_cal_new_5dof_met_std_assisted = get_raw_metabolics(
                full_cal_new_5dof_met, raw_mods)
        raw_full_cal_new_5dof_met_mean_assisted_2 = raw_full_cal_new_5dof_met_mean_assisted.dropna()

        raw_simple_cal_new_minetti_mean_assisted, raw_simple_cal_new_minetti_std_assisted = get_raw_metabolics(
                simple_cal_new_minetti, raw_mods)
        raw_simple_cal_new_minetti_mean_assisted_2 = raw_simple_cal_new_minetti_mean_assisted.dropna()

        raw_full_cal_new_minetti_5dof_mean_assisted, raw_full_cal_new_minetti_5dof_std_assisted = get_raw_metabolics(
                full_cal_new_minetti_5dof, raw_mods)
        raw_full_cal_new_minetti_5dof_mean_assisted_2 = raw_full_cal_new_minetti_5dof_mean_assisted.dropna()

        ## get the uchida models
        raw_full_cal_5dof_uchida_25_mean_assisted, raw_full_cal_5dof_uchida_25_std_assisted = get_raw_metabolics(
                full_cal_5dof_uchida_25, raw_mods)
        raw_full_cal_5dof_uchida_25_mean_assisted_2 = raw_full_cal_5dof_uchida_25_mean_assisted.dropna()
        raw_full_cal_5dof_uchida_25_std_assisted_2 = raw_full_cal_5dof_uchida_25_std_assisted.dropna()


        raw_full_cal_5dof_uchida_30_mean_assisted, raw_full_cal_5dof_uchida_30_std_assisted = get_raw_metabolics(
                full_cal_5dof_uchida_30, raw_mods)
        raw_full_cal_5dof_uchida_30_mean_assisted_2 = raw_full_cal_5dof_uchida_30_mean_assisted.dropna()
        raw_full_cal_5dof_uchida_30_std_assisted_2 = raw_full_cal_5dof_uchida_30_std_assisted.dropna()

        # now make the raw figure
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_genkin_mean,
        #     color='black', marker='o', linestyle='-', label='Bianco 2017')
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_simplemodel_a2_mean_assisted,
        #     color='red', marker='o', linestyle=':', label='simple a2')
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_simple_addedweight_a2_mean_assisted,
        #     color='pink', marker='o', linestyle=':', label='simple add weight a2')
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_simplemodel_met_mean_assisted,
        #     color='mediumblue', marker='o', linestyle=':', label='simple met')
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_simple_addedweight_met_mean_assisted,
        #     color='cornflowerblue', marker='o', linestyle=':', label='simple add weight met')
        # rawmodels_ax_avg.plot(raw_x_pos, 2*raw_simple_cal_met_mean_assisted,
        #     color='forestgreen', marker='o', linestyle=':', label='simple cal. met')
        
        # rawmodels_ax_avg.plot([0,4], 2*raw_simple_cal_addedweight_met_mean_assisted_2,
        #     color='limegreen', marker='o', linestyle=':', label='simple cal. add weight met')
        # rawmodels_ax_avg.plot([0,4], 2*raw_full_met_mean_assisted_2,
        #     color='violet', marker='^', linestyle=':', label='full met 3dof')
        # rawmodels_ax_avg.plot([0,4], 2*raw_full_cal_met_mean_assisted_2,
        #     color='blueviolet', marker='^', linestyle=':', label='full cal. met 3dof')
        # rawmodels_ax_avg.errorbar([0,4], 2*raw_full_cal_5dof_met_mean_assisted_2, raw_full_cal_5dof_met_std_assisted_2, 
        #     color='forestgreen', marker='o', linestyle=':', label='Umberger 2010')
        # rawmodels_ax_avg.plot([0,4], 2*raw_full_cal_new_5dof_met_mean_assisted_2,
        #     color='magenta', marker='^', linestyle=':', label='full cal. fixed 5dof')
        # rawmodels_ax_avg.plot([0,4], 2*raw_simple_cal_new_minetti_mean_assisted_2,
        #     color='lightseagreen', marker='o', linestyle=':', label='simple cal. fixed minetti')
        # rawmodels_ax_avg.plot([0,4], 2*raw_full_cal_new_minetti_5dof_mean_assisted_2,
        #     color='teal', marker='^', linestyle=':', label='full cal. fixed minetti 5dof')

        # rawmodels_ax_avg.errorbar([0,4], 2*raw_full_cal_5dof_uchida_25_mean_assisted_2, raw_full_cal_5dof_uchida_25_std_assisted_2,
        #     color='red', marker='o', linestyle=':', label='Uchida 2016, sigma=.25 MPa')




        # rawmodels_ax_avg.errorbar(raw_x_pos, raw_met_mean,
        #     yerr=raw_met_std, color='black', marker='o',
        #     linestyle='', label='indirect calorimetry')
        # rawmodels_ax_avg.errorbar(raw_x_pos, raw_met_mean,
        #     color='dimgray', marker='.',
        #     linestyle='', label='indirect calorimetry')
        # rawmodels_ax_avg.fill_between(raw_x_pos, raw_met_mean-raw_met_std, raw_met_mean+raw_met_std,
        #     alpha=0.3, edgecolor='dimgray', facecolor='dimgray')

        # rawmodels_ax_avg.errorbar(raw_x_pos, raw_nostand_met_mean,
        #     yerr=raw_nostand_met_std, color='grey', marker='o',
        #     linestyle='', label='indirect calorimetry minus standing')
        rawmodels_ax_avg.errorbar([],[], linestyle='', label='Mean +/- SD')
        rawmodels_ax_avg.errorbar(raw_x_pos, raw_nostand_met_mean,
            color='darkslategrey', marker='o',
            linestyle='', label='Experimental')
        rawmodels_ax_avg.fill_between(raw_x_pos, raw_nostand_met_mean-raw_nostand_met_std, raw_nostand_met_mean+raw_nostand_met_std,
            alpha=0.3, edgecolor='darkslategrey', facecolor='darkslategrey')

        # rawmodels_ax_avg.errorbar([0,4], 2*raw_full_cal_5dof_uchida_30_mean_assisted_2, raw_full_cal_5dof_uchida_30_std_assisted_2,
        #     color='mediumblue', marker='o', linestyle=':', label='Simulation')

        rawmodels_ax_avg.legend(loc='upper left', bbox_to_anchor= (1, 0.9), ncol=1, 
            borderaxespad=0, frameon=False, fontsize=12)
       
        
        #############################################################################
        # combined figures
        #############################################################################
        # now plot the combined figure - individual and average data used from above
        
        # % change figure
        ax_comb.plot(x_pos, met_mean_assisted24, 
            color='red', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted77, 
            color='orange', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted88, 
            color='yellow', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted112, 
            color='green', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted127, 
            color='blue', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted128, 
            color='purple', marker='^', linestyle=':')
        ax_comb.plot(x_pos, met_mean_assisted,
            color='black', marker='o', label='all subject average')
        ax_comb.legend(fontsize=14)


        # plot the raw combined figure - sims met cost
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted24,
            color='red', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted77,
            color='orange', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted88,
            color='yellow', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted112,
            color='green', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted127,
            color='blue', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted128,
            color='purple', marker='o', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_assisted,
            color='blue', marker='o', linestyle='--')

        # # plot the raw activations squared
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted24,
        #     color='red', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted77,
        #     color='orange', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted88,
        #     color='yellow', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted112,
        #     color='green', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted127,
        #     color='blue', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted128,
        #     color='purple', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_assisted,
        #     color='red', marker='o', linestyle='--')

        # raw combined experimental data calorimetry
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW024',:],
            color='red', marker='.', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW077',:],
            color='orange', marker='.', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW088',:],
            color='yellow', marker='.', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW112',:],
            color='green', marker='.', linestyle='-')
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW127',:],
            color='blue', marker='.', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, df_raw_slack.loc['WW128',:],
            color='purple', marker='.', linestyle=':')
        rawboth_ax_comb.plot(raw_x_pos, raw_met_mean,
            color='black', marker='o', linestyle='--')

        # note will maybe want to color code this better later to match the population averge colors

        # # todo add in the nonassist
        #  # plot the raw activations squared - no assist
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted24,
        #     color='red', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted77,
        #     color='orange', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted88,
        #     color='yellow', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted112,
        #     color='green', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted127,
        #     color='blue', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted128,
        #     color='purple', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_default_mean_noassisted,
        #     color='red', marker='o', linestyle='--')

        #  # plot the raw met cost - no assist
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted24,
        #     color='red', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted77,
        #     color='orange', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted88,
        #     color='yellow', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted112,
        #     color='green', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted127,
        #     color='blue', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted128,
        #     color='purple', marker='+', linestyle=':')
        # rawboth_ax_comb.plot(raw_x_pos, 2*raw_met_mean_noassisted,
        #     color='red', marker='o', linestyle='--')

        #############################################################################        

        # by_subject = True
        # if by_subject:
            # met_mean_slack77, met_std_slack77 = get_metabolic_reductions(
            #     file_d
            # ep[0], mods, subject='subject077')
            # met_mean_assisted77, met_std_assisted77 = get_metabolic_reductions(
                # file_dep[1], mods, subject='subject077')
            # met_mean_slack88, met_std_slack88 = get_metabolic_reductions(
            #     file_dep[0], mods, subject='subject088')
            # met_mean_assisted88, met_std_assisted88 = get_metabolic_reductions(
            #     file_dep[1], mods, subject='subject088')
            # met_mean_slack128, met_std_slack128 = get_metabolic_reductions(
                # file_dep[0], mods, subject='subject128')
            # met_mean_assisted128, met_std_assisted128 = \
                # get_metabolic_reductions(file_dep[1], 
                    # mods, subject='subject128')

            # ax.plot(x_pos, met_mean_slack77, 
            #     color='blue', marker='o', linestyle='-')
            # ax.plot(x_pos, met_mean_assisted77, 
            #     color='blue', marker='^', linestyle='--')
            # ax.plot(x_pos, met_mean_slack88, 
            #     color='green', marker='o', linestyle='-')
            # ax.plot(x_pos, met_mean_assisted88, 
            #     color='green', marker='^', linestyle='--')
            # ax.plot(x_pos, met_mean_slack128, 
            #     color='red', marker='o', linestyle='-')
            # ax.plot(x_pos, met_mean_assisted128, 
            #     color='red', marker='^', linestyle='--')

        # else:
            # met_mean_slack, met_std_slack = get_metabolic_reductions(
            #     file_dep[0], mods)
            # met_mean_assisted, met_std_assisted = get_metabolic_reductions(
                # file_dep[1], mods)

            # ax.errorbar(x_pos, met_mean_slack, 
            #     yerr=met_std_slack,
            #     color='blue', marker='o', linestyle='--')
            # ax.errorbar(x_pos, met_mean_assisted, 
                # yerr=met_std_assisted,
                # color='red', marker='o', linestyle='--')

        # for i, v in enumerate(met_relchange_pcent_mean_sort):
        #     color = 'green' if (v < 0) else 'red'
        #     shift = -6 if (v < 0) else 3
        #     ax.text(v + shift, i, '%.2f' % v, va='center',
        #         color=color, fontweight='bold')
        # textstr = 'Subjects included: \n'
        # for subj in list(df_by_subjs.index):
        #     textstr += '  ' + subj + '\n'
        
        # props = dict(boxstyle='round', facecolor='wheat')
        # ax.text(0.10, 0.25, textstr, transform=ax.transAxes, fontsize=14,
        #     bbox=props)
        
        #############################################################################
        # printing and saving of the figures
        #############################################################################



        # make sure to only do the double plot once - target handling
        if 'default' in target[0] and 'noassist' not in target[0]:
            pics = [ax_ind, ax_avg, ax_comb, both_ax_avg]
            figs = [fig_ind, fig_avg, fig_comb, both_fig_avg]
            targ = 0
            for pic in pics:
                pic.set_yticks(np.linspace(-25, 5, 7))
                pic.set_yticklabels(np.linspace(-25, 5, 7), fontsize=12)
                pic.set_xticks(x_pos)
                pic.set_xlim(x_pos[0]-1, x_pos[-1]+1)
                # pic.invert_yaxis()
                pic.set_xticklabels(['LOW','MED','HIGH','MAX'], fontsize=12)
                pic.set_title('Percent Change in Metabolic Rate', fontsize=16)
                pic.grid()
                pic.set_axisbelow(True)
                figs[pics.index(pic)].tight_layout()
                figs[pics.index(pic)].savefig(target[targ], dpi=600)
                figs[pics.index(pic)].savefig(target[targ + 1])
                pl.close(figs[pics.index(pic)])
                targ += 2

        elif 'raw' in target[0] and 'models' not in target[0]:
            print ('****need to check this****')
            pics = [rawboth_ax_avg, rawboth_ax_comb]
            figs = [rawboth_fig_avg, rawboth_fig_comb]
            targ = 0
            for pic in pics:
                pic.set_yticks(np.linspace(0, 6, 7))
                pic.set_yticklabels(np.linspace(0, 6, 7), fontsize=12)
                pic.set_ylabel('Metablic rate [W/kg]', fontsize=12)
                pic.set_xticks(raw_x_pos)
                pic.set_xlim(raw_x_pos[0]-1, raw_x_pos[-1]+1)
                # pic.invert_yaxis()
                pic.set_xticklabels(['SLACK','LOW','MED','HIGH','MAX'], fontsize=12)
                pic.set_xlabel('Assistance Condition', fontsize=12)
                # pic.set_title('Metabolic Rate', fontsize=16)
                pic.grid()
                pic.set_axisbelow(True)
                figs[pics.index(pic)].tight_layout()
                figs[pics.index(pic)].savefig(target[targ], dpi=600)
                figs[pics.index(pic)].savefig(target[targ + 1])
                pl.close(figs[pics.index(pic)])
                targ += 2
        
        elif 'models' in target[0]:
            
            pics = [models_ax_avg, rawmodels_ax_avg]
            figs = [models_fig_avg, rawmodels_fig_avg]
            targ = 0
            for pic in pics:
                if 'raw' in target[targ]:
                    print ('case 1')
                    pic.set_ylim(3, 6)
                    pic.set_yticks(np.linspace(3, 6, 4))
                    pic.set_yticklabels(np.linspace(3, 6, 4), fontsize=12)
                    pic.set_ylabel('Metablic rate [W/kg]', fontsize=12)
                    pic.set_xticks(raw_x_pos)
                    pic.set_xlim(raw_x_pos[0]-0.5, raw_x_pos[-1]+0.5)
                    # pic.invert_yaxis()
                    pic.set_xticklabels(['SLACK','LOW','MED','HIGH','MAX'], fontsize=12)
                    pic.set_xlabel('Assistance Condition', fontsize=12)
                    # pic.set_title('Metabolic Rate', fontsize=16)
                    pic.grid()
                    pic.set_axisbelow(True)
                    figs[pics.index(pic)].tight_layout()
                    figs[pics.index(pic)].savefig(target[targ], dpi=600)
                    figs[pics.index(pic)].savefig(target[targ + 1])
                    pl.close(figs[pics.index(pic)])
                    targ += 2
                else:
                    print ('case 2')
                    # pic.set_ylabel('Percent Metabolic Reduction [%]')
                    pic.set_yticks(np.linspace(-25, 5, 7))
                    pic.set_yticklabels(np.linspace(-25, 5, 7), fontsize=12)
                    pic.set_ylabel('Percent change in Metabolic rate [%]',fontsize=14)
                    pic.set_xticks(qualspos)
                    pic.set_xticklabels(pctplot, fontsize=12)
                    pic.set_xlabel('Metabolic Model',fontsize=14)
                    pic.grid(True)

                    
                    # pic.set_xticks(x_pos)
                    # pic.set_xlim(x_pos[0]-1, x_pos[-1]+1)
                    # pic.invert_yaxis()
                    # pic.set_xticklabels(['LOW','MED','HIGH','MAX'], fontsize=12)
                    # pic.set_title('Percent Change in Metabolic Rate', fontsize=16)
                    # pic.grid()
                    pic.set_axisbelow(True)
                    figs[pics.index(pic)].tight_layout()
                    figs[pics.index(pic)].savefig(target[targ], dpi=600)
                    figs[pics.index(pic)].savefig(target[targ + 1])
                    pl.close(figs[pics.index(pic)])
                    targ += 2

        else:
            pics = [ax_ind, ax_avg, ax_comb]
            figs = [fig_ind, fig_avg, fig_comb]
            targ = 0
            for pic in pics:
                pic.set_yticks(np.linspace(-25, 5, 7))
                pic.set_yticklabels(np.linspace(-25, 5, 7), fontsize=12)
                pic.set_ylabel('Percent change in Metabolic rate [%]',fontsize=12)
                pic.set_xticks(x_pos)
                pic.set_xlim(x_pos[0]-1, x_pos[-1]+1)
                pic.set_xlabel('Assistance Condition',fontsize=12)
                # pic.invert_yaxis()
                pic.set_xticklabels(['LOW','MED','HIGH','MAX'], fontsize=12)
                pic.set_title('Percent Change in Metabolic Rate', fontsize=16)
                pic.grid()
                pic.set_axisbelow(True)
                figs[pics.index(pic)].tight_layout()
                figs[pics.index(pic)].savefig(target[targ], dpi=600)
                figs[pics.index(pic)].savefig(target[targ + 1])
                pl.close(figs[pics.index(pic)])
                targ += 2


        # targ = 0
        # for pic in pics:
        # 	pic.set_yticks(np.linspace(-25, 5, 7))
        # 	pic.set_yticklabels(np.linspace(-25, 5, 7), fontsize=12)
        # 	pic.set_xticks(x_pos)
        # 	pic.set_xlim(x_pos[0]-1, x_pos[-1]+1)
        # 	# pic.invert_yaxis()
        # 	pic.set_xticklabels(['LOW','MED','HIGH','MAX'], fontsize=12)
        # 	pic.set_title('Percent Change in Metabolic Rate', fontsize=16)
        # 	pic.grid()
        # 	pic.set_axisbelow(True)
        # 	figs[pics.index(pic)].tight_layout()
        # 	figs[pics.index(pic)].savefig(target[targ], dpi=600)
        # 	figs[pics.index(pic)].savefig(target[targ + 1])
        # 	pl.close(figs[pics.index(pic)])
        # 	targ += 2

        # adding in the figure details
        # ax.set_yticks(np.linspace(-25,5,7))
        # ax.set_yticklabels(np.linspace(-25,5,7))
        # ax.set_xticks(x_pos)
        # ax.set_xlim(x_pos[0]-1, x_pos[-1]+1)
        # ax.invert_yaxis()
        # ax.set_xticklabels(['LOW','MED','HIGH','MAX'])
        # ax.set_title('Percent Change in Metabolic Rate')
        # ax.grid()
        # ax.set_axisbelow(True)
        # fig.tight_layout()
        # fig.savefig(target[0], dpi=600)
        # fig.savefig(target[1])
        # pl.close(fig)


def aggregate_moments(file_dep, target, cond_name, cycles):
    ## Percent gait cycle.
    num_time_points = 400
    pgc = np.linspace(0, 100, num_time_points)

    muscle_names = None

    subject_array = list()
    cycle_array = list()
    dof_array = list()
    muscle_array = list()
    all_data = list()
    for icycle, fpath in enumerate(file_dep):
        cycle = cycles[cond_name][icycle]
        #for cycle_info in cycles:
        cycle_df = pd.read_csv(fpath, index_col=0, header=[0, 1], skiprows=1)
        # Convert time to percent gait cycle.
        x = np.linspace(cycle.gl.cycle_start,
                cycle.gl.cycle_end, num_time_points)
        for column in cycle_df.columns:
            dof, actuator = column
            subject_array.append(cycle.subject.name)
            cycle_array.append(cycle.id)
            dof_array.append(dof)
            muscle_array.append(actuator)
            moment = np.interp(pgc, cycle_df.index, cycle_df[column])
            all_data.append(moment)
    # Convert from (n_cycles * n_dofs * n_muscles) x n_times
    #         to   n_times x (n_cycles * n_dofs * n_muscles)
    all_data_array = np.array(all_data).transpose()

    multiindex_arrays = [subject_array, cycle_array, dof_array, muscle_array]
    columns = pd.MultiIndex.from_arrays(multiindex_arrays,
            names=['subject', 'cycle', 'DOF', 'actuator'])

    all_data_df = pd.DataFrame(all_data_array, columns=columns, index=pgc)
    target_dir = os.path.dirname(target[0])
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    with file(target[0], 'w') as f:
        f.write('# all columns are moments normalized by subject '
                'mass (N-m/kg).\n')
        all_data_df.to_csv(f)
    # How to read this in: df.read_csv(..., index_col=0, header=[0, 1, 2, 3],
    #                                  skiprows=1)

class TaskAggregateMomentsExperiment(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, subject=None, conditions=['slack'],
                 cycle_nums=['cycle%02i' % i for i in range(1,11)]):
        super(TaskAggregateMomentsExperiment, self).__init__(study)
        self.name = 'aggregate_experiment_moments'
        self.doc = 'Aggregate no-mod actuator moments into a data file.'

        suffix = ''
        self.costdir = ''
        if not (study.costFunction == 'Default'):
            suffix += '_%s' % study.costFunction
            self.costdir = study.costFunction 

        if subject == None:
            subjects = [s.name for s in study.subjects]
        else:
            suffix += '_%s' % subject
            self.name += '_%s' % subject
            subjects = [subject]

        self.cycles = dict()
        for cond_name in conditions:
            self.cycles[cond_name] = list()
            deps = []
            for subject in study.subjects:
                if not subject.name in subjects: continue
                cond = subject.get_condition(cond_name)
                if not cond: continue
                # We know there is only one overground trial, but perhaps it
                # has not yet been added for this subject.
                assert len(cond.trials) <= 1
                if len(cond.trials) == 1:
                    trial = cond.trials[0]
                    for cycle in trial.cycles:
                        if cycle.name not in cycle_nums: continue
                        self.cycles[cond_name].append(cycle)

                        # Moment file paths
                        fpath = os.path.join(trial.results_exp_path, 'mrs',
                            cycle.name, self.costdir,
                             '%s_%s_mrs_moments.csv' % (study.name, cycle.id))
                        deps.append(fpath)

            self.add_action(deps,
                    [
                        os.path.join(study.config['results_path'], 
                            'experiments',
                            'experiment_%s_moments%s.csv' % (cond_name, suffix)),
                        ],
                    aggregate_moments, cond_name, self.cycles)


class TaskAggregateMomentsMod(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, mod_name, subject=None, 
                 cycle_nums=['cycle%02i' % i for i in range(1,11)]):
        super(TaskAggregateMomentsMod, self).__init__(study)
        self.name = 'aggregate_mod_moments_%s' % mod_name
        self.doc = 'Aggregate actuator moments into a data file.'
        self.mod_name = mod_name

        conditions = list()
        # conditions.append('slack') # uncomment for the slack kin with assistance
        conditions.append(mod_name)
        self.conditions = conditions

        suffix = ''
        self.costdir = ''
        if not (study.costFunction == 'Default'):
            suffix += '_%s' % study.costFunction
            self.costdir = study.costFunction 

        if subject == None:
            subjects = [s.name for s in study.subjects]
        else:
            suffix += '_%s' % subject
            self.name += '_%s' % subject
            subjects = [subject]

        self.cycles = dict()

        for cond_name in self.conditions:
            self.cycles[cond_name] = list()
            deps = []
            for subject in study.subjects:
                if not subject.name in subjects: continue
                cond = subject.get_condition(cond_name)
                if not cond: continue
                # We know there is only one overground trial, but perhaps it
                # has not yet been added for this subject.
                assert len(cond.trials) <= 1
                if len(cond.trials) == 1:
                    trial = cond.trials[0]
                    for cycle in trial.cycles:
                        if cycle.name not in cycle_nums: continue
                        self.cycles[cond_name].append(cycle)

                        deps.append(os.path.join(
                                self.study.config['results_path'],
                                'mrsmod_%s' % mod_name,
                                trial.rel_path, 'mrs', cycle.name, self.costdir,
                                '%s_%s_mrs_moments.csv' % (study.name,
                                    cycle.id))
                                )

            self.add_action(deps,
                    [os.path.join(study.config['results_path'], 
                        'mrsmod_%s' % mod_name,
                        'mod_%s_%s_moments%s.csv' % (mod_name,
                            cond_name, suffix)),
                        ],
                    aggregate_moments, cond_name, self.cycles)


class TaskPlotMoments(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, agg_task, mod_name=None, subject=None):
        super(TaskPlotMoments, self).__init__(study)
        task_name = 'experiment' if mod_name==None else mod_name
        self.name = 'plot_moment_breakdown_%s' % task_name
        self.doc = 'Plot joint moment breakdown by muscle and device moments'


        # print '\n\ngoing to get into things here to create a new task probably'


        if subject:
            self.name += '_%s' % subject

        for icond, agg_target in enumerate(agg_task.targets):
            # This assumes csv_task.targets and csv_task.cycles hold cycles in
            # the same order.

            # print 'icond'
            # print icond
            # print 'agg_target'
            # print agg_target

            # self.agg_target = agg_target
            self.add_action([agg_target], [],
                    self.plot_joint_moment_breakdown)

            # self.actions += [self.plot_joint_moment_breakdown]

    def plot_joint_moment_breakdown(self, file_dep, target):

        # print 'in the second function'
        # print 'file_dep'
        # print file_dep
        # print 'target'
        # print target

        # print '\nso it seems like we are looping through each subject to add all their data to a singular spreadsheet'
        # print '\nneed to figure out where the individual figures are generated and saved'

        # going to get indipendent data frames for exp and max
        slack_path = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/results/experiments/experiment_slack_moments_Met.csv'
        df_exp_all = pd.read_csv(slack_path, index_col=0,
                    header=[0,1,2,3], skiprows=1)
        df_exp_by_subj_dof_musc = df_exp_all.groupby(
                    level=['subject','DOF','actuator'], axis=1).mean()
        df_exp_mean = df_exp_by_subj_dof_musc.groupby(level=['DOF','actuator'],
                    axis=1).mean()
        df_exp_std = df_exp_by_subj_dof_musc.groupby(level=['DOF','actuator'],
                    axis=1).std()
        exp_pgc = df_exp_mean.index



        # getting the current dataframe
        df_all = pd.read_csv(file_dep[0], index_col=0,
                header=[0, 1, 2, 3], skiprows=1)
        # Average over cycles.
        # axis=1 for columns (not rows).
        
        print (df_all)

        df_by_subj_dof_musc = df_all.groupby(
                level=['subject', 'DOF', 'actuator'], axis=1).mean()
        
        print (df_by_subj_dof_musc)
        df_mean = df_by_subj_dof_musc.groupby(level=['DOF', 'actuator'],
                axis=1).mean()
        print (df_mean)
        df_std = df_by_subj_dof_musc.groupby(level=['DOF', 'actuator'],
                axis=1).std()
        pgc = df_mean.index


        import seaborn.apionly as sns
        palette2 = sns.color_palette('colorblind', 19) # 9
        palette = ['firebrick','mediumblue','forestgreen','gold','plum','orangered','c','red',
                    'dodgerblue','darkorchid','mediumseagreen','indianred','hotpink','orange','grey','sandybrown',
                    'slateblue','darkgreen','darkslategrey']
        ##TODO: the issue is that muscles doesn't have them all and then colors is too
        # short as a result

        # muscles = ['glut_max2_r', 'psoas_r', 'semimem_r', 'rect_fem_r',
        #         'bifemsh_r', 'vas_int_r', 'med_gas_r', 'soleus_r', 'tib_ant_r']
        
        muscles = ['add_brev_r',
                'add_long_r',
                'add_mag3_r',
                'add_mag4_r',
                'add_mag2_r',
                'add_mag1_r',
                'bifemlh_r',
                'bifemsh_r',
                'ext_dig_r',
                'ext_hal_r',
                'flex_dig_r',
                'flex_hal_r',
                'lat_gas_r',
                'med_gas_r',
                'glut_max1_r',
                'glut_max2_r',
                'glut_max3_r',
                'glut_med1_r',
                'glut_med2_r',
                'glut_med3_r',
                'glut_min1_r',
                'glut_min2_r',
                'glut_min3_r',
                'grac_r',
                'iliacus_r',
                'per_brev_r',
                'per_long_r',
                'peri_r',
                'psoas_r',
                'rect_fem_r',
                'sar_r',
                'semimem_r',
                'semiten_r',
                'soleus_r',
                'tfl_r',
                'tib_ant_r',
                'tib_post_r',
                'vas_int_r',
                'vas_lat_r',
                'vas_med_r']



        fig = pl.figure(figsize=(9, 3.75))
        dof_names = ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'hip_adduction_r']
        ylabels = ['hip extension', 'knee extension', 'ankle plantarflexion', 'hip abduction']


        # nice_act_names = {
        #         'glut_max2_r': 'glut. max.',
        #         'psoas_r': 'iliopsoas',
        #         'semimem_r': 'hamstrings',
        #         'rect_fem_r': 'rect. fem.',
        #         'bifemsh_r': 'bi. fem. s.h.',
        #         'vas_int_r': 'vasti',
        #         'med_gas_r': 'gastroc.',
        #         'soleus_r': 'soleus',
        #         'tib_ant_r': 'tib. ant.',
        #         'net': 'net',
        #         'active': 'active',
        #         'passive': 'passive',
        #         }

        ## TODO add in all the names of the muscles
        # nice_act_names = {
        #         'add_brev_r':'add. brev.',
        #         'add_long_r':'add. long.',
        #         'add_mag3_r':'add. mag. 3',
        #         'add_mag4_r':'add. mag. 4',
        #         'add_mag2_r':'add. mag. 2',
        #         'add_mag1_r':'add. mag. 1',
        #         'bifemlh_r':'bifem. lh.',
        #         'bifemsh_r':'bifem. sh.',
        #         'ext_dig_r':'ext. dig.',
        #         'ext_hal_r':'ext. hal.',
        #         'flex_dig_r':'flex. dig.',
        #         'flex_hal_r':'flex. hal.',
        #         'lat_gas_r':'lat. gastroc',
        #         'med_gas_r':'med. gastroc',
        #         'glut_max1_r':'glut. max. 1',
        #         'glut_max2_r':'glut. max 2',
        #         'glut_max3_r':'glut. max 3',
        #         'glut_med1_r':'glut. med 1',
        #         'glut_med2_r':'glut. med 2',
        #         'glut_med3_r':'glut. med 3',
        #         'glut_min1_r':'glut. min 1',
        #         'glut_min2_r':'glut. min 2',
        #         'glut_min3_r':'glut. min 3',
        #         'grac_r':'gracilis',
        #         'iliacus_r':'iliacus',
        #         'per_brev_r':'per. brev.',
        #         'per_long_r':'per. long.',
        #         'peri_r':'peri.',
        #         'psoas_r':'psoas',
        #         'rect_fem_r':'rec. fem.',
        #         'sar_r':'sartorius',
        #         'semimem_r':'semimem.',
        #         'semiten_r':'semiten',
        #         'soleus_r':'soleus',
        #         'tfl_r':'tfl.',
        #         'tib_ant_r':'tib. ant.',
        #         'tib_post_r':'tib. post.',
        #         'vas_int_r':'vas. int.',
        #         'vas_lat_r':'vas. lat.',
        #         'vas_med_r':'vas. med.',
        #         'net':'net',
        #         'active': 'active',
        #         'passive': 'passive',

        # }

        ## have to add this for the apoorva model
        nice_act_names = {
                'addbrev_r':'add. brev.',
                'addlong_r':'add. long.',
                'addmagDist_r':'add. mag. 3',
                'addmagIsch_r':'add. mag. 4',
                'addmagMid_r':'add. mag. 2',
                'addmagProx_r':'add. mag. 1',
                'bflh_r':'bifem. lh.',
                'bfsh_r':'bifem. sh.',
                'edl_r':'ext. dig.',
                'ehl_r':'ext. hal.',
                'fdl_r':'flex. dig.',
                'fhl_r':'flex. hal.',
                'gaslat_r':'lat. gastroc',
                'gasmed_r':'med. gastroc',
                'glmax1_r':'glut. max. 1',
                'glmax2_r':'glut. max 2',
                'glmax3_r':'glut. max 3',
                'glmed1_r':'glut. med 1',
                'glmed2_r':'glut. med 2',
                'glmed3_r':'glut. med 3',
                'glmin1_r':'glut. min 1',
                'glmin2_r':'glut. min 2',
                'glmin3_r':'glut. min 3',
                'grac_r':'gracilis',
                'iliacus_r':'iliacus',
                'perbrev_r':'per. brev.',
                'perlong_r':'per. long.',
                'piri_r':'peri.',
                'psoas_r':'psoas',
                'recfem_r':'rec. fem.',
                'sart_r':'sartorius',
                'semimem_r':'semimem.',
                'semiten_r':'semiten',
                'soleus_r':'soleus',
                'tfl_r':'tfl.',
                'tibant_r':'tib. ant.',
                'tibpost_r':'tib. post.',
                'vasint_r':'vas. int.',
                'vaslat_r':'vas. lat.',
                'vasmed_r':'vas. med.',
                'net':'net',
                'active': 'active',
                'passive': 'passive',
                'add_brev_r':'add. brev.',
                'add_long_r':'add. long.',
                'add_mag3_r':'add. mag. 3',
                'add_mag4_r':'add. mag. 4',
                'add_mag2_r':'add. mag. 2',
                'add_mag1_r':'add. mag. 1',
                'bifemlh_r':'bifem. lh.',
                'bifemsh_r':'bifem. sh.',
                'ext_dig_r':'ext. dig.',
                'ext_hal_r':'ext. hal.',
                'flex_dig_r':'flex. dig.',
                'flex_hal_r':'flex. hal.',
                'lat_gas_r':'lat. gastroc',
                'med_gas_r':'med. gastroc',
                'glut_max1_r':'glut. max. 1',
                'glut_max2_r':'glut. max 2',
                'glut_max3_r':'glut. max 3',
                'glut_med1_r':'glut. med 1',
                'glut_med2_r':'glut. med 2',
                'glut_med3_r':'glut. med 3',
                'glut_min1_r':'glut. min 1',
                'glut_min2_r':'glut. min 2',
                'glut_min3_r':'glut. min 3',
                'grac_r':'gracilis',
                'iliacus_r':'iliacus',
                'per_brev_r':'per. brev.',
                'per_long_r':'per. long.',
                'peri_r':'peri.',
                'psoas_r':'psoas',
                'rect_fem_r':'rec. fem.',
                'sar_r':'sartorius',
                'semimem_r':'semimem.',
                'semiten_r':'semiten',
                'soleus_r':'soleus',
                'tfl_r':'tfl.',
                'tib_ant_r':'tib. ant.',
                'tib_post_r':'tib. post.',
                'vas_int_r':'vas. int.',
                'vas_lat_r':'vas. lat.',
                'vas_med_r':'vas. med.',
                
        }


        act_names = df_mean.columns.levels[1]
        # print 'here are the act names'
        # print act_names

        act_exp_names = df_exp_mean.columns.levels[1]
        # print 'here are the exp act names'
        # print act_exp_names

        # print 'should be able to just create a list of muscles that I dont want to plot and check from it'
        plotlist = ['add_long_r','bifemsh_r','lat_gas_r','med_gas_r','glut_max2_r','glut_med1_r','glut_med2_r','glut_med3_r','iliacus_r',
                    'psoas_r','rect_fem_r','semimem_r','soleus_r','tfl_r','tib_ant_r','vas_lat_r','bifemlh_r','add_brev_r','add_mag1_r',
                    ]
        noplotlist = ['sar_r','peri_r','glut_max1_r','vas_int_r','vas_med_r','glut_min1_r','glut_max3_r']

        hipextnpl = ['tfl_r','glut_med2_r','glut_med3_r','add_long_r','add_brev_r']
        kneeextnpl = ['tfl_r','grac_r','med_gas_r']
        anklenpl = ['med_gas_r']
        hipaddnpl = ['psoas_r','add_long_r','rect_fem_r','glut_max2_r','semimem_r','add_brev_r','add_mag1_r']


        # colors = {plotlist[i]: palette[i] for i in range(19)}
        # colors['net'] = 'black' #(0.7,) * 3 # light gray
        # colors['active'] = 'orange'
        # colors['passive'] = 'blue'

        # print colors

        def plot(column_key, act_name):
            if column_key in df_mean.columns:
                y_mean = -df_mean[column_key]
                y_std = -df_std[column_key]
                
                if act_name == 'lat_gas_r':
                    # print('got it')
                    second = -df_mean[('ankle_angle_r','med_gas_r')]
                    y_mean = y_mean + second
                    # print second
                    # print y_mean

                if act_name == 'active' or act_name == 'passive':
                    ax.plot(pgc, y_mean, #color=colors[act_name],
                            label=nice_act_names[act_name],
                            linestyle='--')
                    # ax.fill_between(pgc, y_mean-y_std, y_mean+y_std,
                    #     alpha=0.3, color=colors[act_name])
                else:
                    if act_name == 'lat_gas_r':
                        ax.plot(pgc, y_mean, #color=colors[act_name],
                            label='gastroc_r')
                    else:
                        ax.plot(pgc, y_mean, #color=colors[act_name],
                            label=nice_act_names[act_name])


        def plotcompare(column_key, act_name, idof, starteridx):

            # print '\nprep yoself ho, its going down'
            # print 'column key'
            # print column_key
            # print column_key[0]
            # print 'act name'
            # print act_name
            # print 'idof'
            # print idof
            # print 'starter idx'
            # print starteridx

            if column_key in df_mean.columns:
                y_max_mean = -df_mean[column_key]
                y_max_std = -df_std[column_key]
                if column_key in df_exp_mean.columns:
                    y_exp_mean = -df_exp_mean[column_key]
                    y_exp_std = -df_exp_std[column_key]

                if act_name == 'lat_gas_r':
                    if column_key[0] == 'ankle_angle_r':
                        second_exp = -df_exp_mean[('ankle_angle_r','med_gas_r')]
                        y_exp_mean = y_exp_mean + second_exp
                        second_max = -df_mean[('ankle_angle_r','med_gas_r')]
                        y_max_mean = y_max_mean + second_max
                    elif column_key[0] == 'knee_angle_r':
                        second_exp = -df_exp_mean[('knee_angle_r','med_gas_r')]
                        y_exp_mean = y_exp_mean + second_exp
                        second_max = -df_mean[('knee_angle_r','med_gas_r')]
                        y_max_mean = y_max_mean + second_max

                if act_name == 'active' or act_name == 'passive':
                    axcomp[idof, starteridx].plot(exp_pgc, y_max_mean, #color=colors[act_name],
                            # label=nice_act_names[act_name],
                            linestyle='--', label='%s max' % act_name)
                    # ax.fill_between(pgc, y_mean-y_std, y_mean+y_std,
                    #     alpha=0.3, color=colors[act_name])
                    # axcomp[idof, starteridx].plot(exp_pgc, y_exp_mean, color=colors[act_name],
                    #         # label=nice_act_names[act_name],
                    #         linestyle='-', label='%s slack' % act_name)
                    axcomp[idof, starteridx].legend()
                else:
                    if act_name == 'lat_gas_r':
                        axcomp[idof, starteridx].plot(exp_pgc, y_max_mean, #color=colors[act_name],
                            # label='gastroc_r', 
                            linestyle='--', label='gastroc_r max')
                        axcomp[idof, starteridx].plot(exp_pgc, y_exp_mean, #color=colors[act_name],
                            # label='gastroc_r', 
                            linestyle='-', label='gastroc_r slack')
                        axcomp[idof, starteridx].legend()

                    else:
                        axcomp[idof, starteridx].plot(exp_pgc, y_max_mean, #color=colors[act_name],
                            # label=nice_act_names[act_name], 
                            linestyle='--', label='%s max' % act_name)
                        axcomp[idof, starteridx].plot(exp_pgc, y_exp_mean, #color=colors[act_name],
                            # label=nice_act_names[act_name], 
                            linestyle='-', label='%s slack' % act_name)
                        axcomp[idof, starteridx].legend()
            axcomp[idof, starteridx].legend(frameon=False, fontsize=10)
            axcomp[idof, starteridx].set_xlim(0,100)
            axcomp[idof, starteridx].set_ylim(-1.1,2.0)
            axcomp[idof, starteridx].set_ylabel('%s (N-m/kg)' % ylabels[idof])
            axcomp[idof, starteridx].set_xlabel('time (% gait cycle)')
            axcomp[idof, starteridx].spines['right'].set_visible(False)
            axcomp[idof, starteridx].spines['top'].set_visible(False)
            axcomp[idof, starteridx].xaxis.set_ticks_position('bottom')
            axcomp[idof, starteridx].yaxis.set_ticks_position('left')
        ### end 



        casecheck = file_dep[0][-19:-16]
        figcomp, axcomp  = pl.subplots(4,18,figsize=(40,15))


        for idof, dof_name in enumerate(dof_names):                 # going for each dof here
            ax = fig.add_subplot(1, len(dof_names), idof + 1)       # create a subplot for that dof
            ax.axhline(color='k', linewidth=0.5, zorder=0)          # creating the axis
            plot((dof_name, 'net'), 'net')
            
            # now to handle the individual plot for net
            starteridx = 0
            if casecheck == 'max':
                # axcompare = figcompare.add_subplot(4, 19, starteridx)
                # axcompare.axhline(color='k', linewidth=0.5, zorder=0)
                axcomp[idof, starteridx].axhline(color='k', linewidth=0.5, zorder=0)
                # starthope += 1
                plotcompare((dof_name, 'net'), 'net', idof, starteridx)
                starteridx += 1

            # looping through the muscles specified on no plot list
            for iact, act_name in enumerate(act_names):
                # if dof_name == 'hip_flexion_r':
                #     if act_name in hipextnpl: continue
                # if dof_name == 'knee_angle_r':
                #     if act_name in kneeextnpl: continue
                # if dof_name == 'ankle_angle_r':
                #     if act_name in anklenpl: continue
                # if dof_name == 'hip_adduction_r':
                    # if act_name in hipaddnpl: continue
                # if act_name in noplotlist: continue 
                if act_name == 'net': continue
                if act_name != 'active':
                    if act_name not in df_exp_mean[dof_name].columns or act_name not in df_mean[dof_name].columns:
                        continue

                # plot the combined figures for each actuator
                column_key = (dof_name, act_name)
                plot(column_key, act_name)

                # now need to create a new subplot for each actuator in the individual figure
                if casecheck == 'max':
                    # axcompare = figcompare.add_subplot(4,19, starteridx)
                    # axcompare.axhline(color='k', linewidth=0.5, zorder=0)
                    axcomp[idof, starteridx].axhline(color='k', linewidth=0.5, zorder=0)
                    plotcompare(column_key, act_name, idof, starteridx)
                    starteridx += 1


            # set up the general figure
            ax.legend(frameon=False, fontsize=6)
            ax.set_xlim(0, 100)
            ax.set_ylim(-1.1, 2.0)
            # if idof > 0:
            #     ax.set_yticklabels([])
            ax.set_ylabel('%s (N-m/kg)' % ylabels[idof])
            ax.set_xlabel('time (% gait cycle)')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            
        # finish and save the general figure
        fig.tight_layout()
        fig.savefig(file_dep[0].replace('.csv', '.pdf'))
        fig.savefig(file_dep[0].replace('.csv', '.png'), dpi=600)
        pl.close(fig)


        # finish and save the compare individual figure
        if casecheck == 'max':
            extrafig_target1 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/results/mrsmod_max/comparison_moments_plots.pdf'
            extrafig_target2 = 'C:/Users/JP/code/repos/Stanford/delplab/projects/AnkleHipExosuit/results/mrsmod_max/comparison_moments_plots.png'
            figcomp.tight_layout()
            # figcompare.savefig(extrafig_target1)
            # figcompare.savefig(extrafig_target2)
            figcomp.savefig(extrafig_target1)
            figcomp.savefig(extrafig_target2)
            pl.close(fig)


class TaskAggregateMuscleActivity(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, mod_name=None, subject=None,
                 cycle_nums=['cycle%02i' % i for i in range(1,11)]):
        super(TaskAggregateMuscleActivity, self).__init__(study)
        

        suffix = ''
        conditions = list()
        if mod_name == None:
            mod_fpath = 'experiments'
            mod_name = 'experiments'
            conditions.append('slack')
        else:
            mod_fpath = 'mrsmod_%s' % mod_name
            suffix += '_%s' % mod_name
            # conditions.append('slack')
            conditions.append(mod_name)


        self.costdir = ''   
        if not (study.costFunction == 'Default'):
            suffix += '_%s' % study.costFunction
            self.costdir = study.costFunction 

        if subject == None:
            subjects = [s.name for s in study.subjects]
        else:
            suffix += '_%s' % subject
            subjects = [subject]
        
        self.name = 'aggregate_muscle_activity%s' % suffix
        self.doc = 'Aggregate muscle activity into a data file.'

        self.cycles = dict()
        for cond_name in conditions:
            self.cycles[cond_name] = list()
            deps = []
            for subject in study.subjects:
                if not subject.name in subjects: continue
                cond = subject.get_condition(cond_name)
                if not cond: continue
                # We know there is only one overground trial, but perhaps it
                # has not yet been added for this subject.
                assert len(cond.trials) <= 1
                if len(cond.trials) == 1:
                    trial = cond.trials[0]
                    for cycle in trial.cycles:
                        if cycle.name not in cycle_nums: continue
                        self.cycles[cond_name].append(cycle)

                        # Results MAT file paths
                        fpath = os.path.join(study.config['results_path'], 
                            mod_fpath, subject.name, cond.name, 'mrs', 
                            cycle.name, self.costdir,
                            '%s_%s_mrs.mat' % (study.name, cycle.id))
                        deps.append(fpath)

            self.add_action(deps,
                    [
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_excitations%s.csv' % (
                                mod_name, cond_name, suffix)),
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_activations%s.csv' % (
                                mod_name, cond_name, suffix)),
                        ],
                    self.aggregate_muscle_activity, cond_name, self.cycles)

    def aggregate_muscle_activity(self, file_dep, target, cond_name, cycles):
        num_time_points = 400
        pgc = np.linspace(0, 100, num_time_points)

        muscle_names = None

        subject_array = list()
        cycle_array = list()
        muscle_array = list()
        all_exc = list()
        all_act = list()
        for icycle, fpath in enumerate(file_dep):
            cycle = cycles[cond_name][icycle]            
            muscle_names = util.hdf2list(fpath, 'MuscleNames', type=str)
            exc_df = util.hdf2pandas(fpath,
                'MExcitation', labels=muscle_names)
            act_df = util.hdf2pandas(fpath,
                'MActivation', labels=muscle_names)

            exc_index = np.linspace(0, 100, len(exc_df.index.values))
            act_index = np.linspace(0, 100, len(act_df.index.values))
            for muscle in exc_df.columns:
                subject_array.append(cycle.subject.name)
                cycle_array.append(cycle.id)
                muscle_array.append(muscle)
                exc = np.interp(pgc, exc_index, exc_df[muscle])
                act = np.interp(pgc, act_index, act_df[muscle])
                all_exc.append(exc)
                all_act.append(act)


        all_exc_array = np.array(all_exc).transpose()
        all_act_array = np.array(all_act).transpose()


        multiindex_arrays = [subject_array, cycle_array, muscle_array]
        columns = pd.MultiIndex.from_arrays(multiindex_arrays,
                names=['subject', 'cycle', 'muscle'])


        all_exc_df = pd.DataFrame(all_exc_array, columns=columns, index=pgc)
        # all_exc_df_transpose = all_exc_df.transpose()


        # print('all first')
        # print(all_exc_df)
        # print('transpose')
        # print(all_exc_df_transpose)

        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# all columns are muscle excitations.\n')
            all_exc_df.to_csv(f)


        all_act_df = pd.DataFrame(all_act_array, columns=columns, index=pgc)
        # all_act_df_transpose = all_act_df.transpose()
        target_dir = os.path.dirname(target[1])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[1], 'w') as f:
            f.write('# all columns are muscle activations.\n')
            all_act_df.to_csv(f)
            # all_act_df_transpose.to_csv(f)
        # How to read this in: df.read_csv(..., index_col=0, header=[0, 1, 2],
        #                                  skiprows=1)


class TaskAggregateNormalizedMuscleDynamics(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, mod_name=None, subject=None,
                 cycle_nums=['cycle%02i' % i for i in range(1,11)]):
        super(TaskAggregateNormalizedMuscleDynamics, self).__init__(study)
        
        suffix = ''
        conditions = list()
        if mod_name == None:
            mod_fpath = 'experiments'
            mod_name = 'experiments'
            conditions.append('slack')
        else:
            mod_fpath = 'mrsmod_%s' % mod_name
            suffix += '_%s' % mod_name
            # conditions.append('slack')
            conditions.append(mod_name)

        self.costdir = ''   
        if not (study.costFunction == 'Default'):
            suffix += '_%s' % study.costFunction
            self.costdir = study.costFunction 

        if subject == None:
            subjects = [s.name for s in study.subjects]
        else:
            suffix += '_%s' % subject
            subjects = [subject]
        
        self.name = 'aggregate_norm_muscle_dynamics%s' % suffix
        self.doc = 'Aggregate normalized muscle dynamics into a data file.'

        self.cycles = dict()
        for cond_name in conditions:
            self.cycles[cond_name] = list()
            deps = []
            for subject in study.subjects:
                if not subject.name in subjects: continue
                cond = subject.get_condition(cond_name)
                if not cond: continue
                # We know there is only one overground trial, but perhaps it
                # has not yet been added for this subject.
                assert len(cond.trials) <= 1
                if len(cond.trials) == 1:
                    trial = cond.trials[0]
                    for cycle in trial.cycles:
                        if cycle.name not in cycle_nums: continue
                        self.cycles[cond_name].append(cycle)

                        # Results MAT file paths
                        fpath = os.path.join(study.config['results_path'], 
                            mod_fpath, subject.name, cond.name, 'mrs', 
                            cycle.name, self.costdir,
                            '%s_%s_mrs.mat' % (study.name, cycle.id))
                        deps.append(fpath)

            self.add_action(deps,
                    [
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_norm_active_force%s.csv' % (
                                mod_name, cond_name, suffix)),
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_norm_passive_force%s.csv' % (
                                mod_name, cond_name, suffix)),
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_norm_tendon_force%s.csv' % (
                                mod_name, cond_name, suffix)),
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_norm_fiber_length%s.csv' % (
                                mod_name, cond_name, suffix)),
                        os.path.join(study.config['results_path'], 
                            mod_fpath,'%s_%s_norm_fiber_velocity%s.csv' % (
                                mod_name, cond_name, suffix)),
                        # TODO:
                        # os.path.join(study.config['results_path'], 
                        #     mod_fpath,'%s_%s_norm_tendon_length%s.csv' % (
                        #         mod_name, cond_name, suffix)),
                        ],
                    self.aggregate_muscle_activity, cond_name, self.cycles)

    def aggregate_muscle_activity(self, file_dep, target, cond_name, cycles):

        num_time_points = 400
        pgc = np.linspace(0, 100, num_time_points)

        muscle_names = None

        subject_array = list()
        cycle_array = list()
        muscle_array = list()
        all_fce = list()
        all_fpe = list()
        all_fse = list()
        all_lMtilde = list()
        all_vMtilde = list()
        # TODO all_lTtilde = list()
        for icycle, fpath in enumerate(file_dep):
            cycle = cycles[cond_name][icycle]

            muscle_names = util.hdf2list(fpath, 'MuscleNames', type=str)
            fce_df = util.hdf2pandas(fpath,
                'MuscleData/fce', labels=muscle_names)
            fpe_df = util.hdf2pandas(fpath,
                'MuscleData/fpe', labels=muscle_names)
            fse_df = util.hdf2pandas(fpath,
                'MuscleData/fse', labels=muscle_names)
            lMtilde_df = util.hdf2pandas(fpath,
                'MuscleData/lMtilde', labels=muscle_names)
            vMtilde_df = util.hdf2pandas(fpath,
                'MuscleData/vMtilde', labels=muscle_names)
            # lTtilde_df = util.hdf2pandas(fpath,
            #     'MuscleData/lTtilde', labels=muscle_names)

            fce_index = np.linspace(0, 100, len(fce_df.index.values))
            fpe_index = np.linspace(0, 100, len(fpe_df.index.values))
            fse_index = np.linspace(0, 100, len(fse_df.index.values))
            lMtilde_index = np.linspace(0, 100, len(lMtilde_df.index.values))
            vMtilde_index = np.linspace(0, 100, len(vMtilde_df.index.values))
            # lTtilde_index = np.linspace(0, 100, len(lTtilde_df.index.values))
            for muscle in fce_df.columns:
                subject_array.append(cycle.subject.name)
                cycle_array.append(cycle.id)
                muscle_array.append(muscle)
                fce = np.interp(pgc, fce_index, fce_df[muscle])
                fpe = np.interp(pgc, fpe_index, fpe_df[muscle])
                fse = np.interp(pgc, fse_index, fse_df[muscle])
                lMtilde = np.interp(pgc, lMtilde_index, lMtilde_df[muscle])
                vMtilde = np.interp(pgc, vMtilde_index, vMtilde_df[muscle])
                # lTtilde = np.interp(pgc, lTtilde_index, lTtilde_df[muscle])
                all_fce.append(fce)
                all_fpe.append(fpe)
                all_fse.append(fse)
                all_lMtilde.append(lMtilde)
                all_vMtilde.append(vMtilde)
                # all_lTtilde.append(lTtilde)

        all_fce_array = np.array(all_fce).transpose()
        all_fpe_array = np.array(all_fpe).transpose()
        all_fse_array = np.array(all_fse).transpose()
        all_lMtilde_array = np.array(all_lMtilde).transpose()
        all_vMtilde_array = np.array(all_vMtilde).transpose()
        # all_lTtilde_array = np.array(all_lTtilde).transpose()

        multiindex_arrays = [subject_array, cycle_array, muscle_array]
        columns = pd.MultiIndex.from_arrays(multiindex_arrays,
                names=['subject', 'cycle', 'muscle'])

        all_fce_df = pd.DataFrame(all_fce_array, columns=columns, index=pgc)
        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# all columns are normalized active force.\n')
            all_fce_df.to_csv(f)

        all_fpe_df = pd.DataFrame(all_fpe_array, columns=columns, index=pgc)
        target_dir = os.path.dirname(target[1])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[1], 'w') as f:
            f.write('# all columns are normalized passive force.\n')
            all_fpe_df.to_csv(f)

        all_fse_df = pd.DataFrame(all_fse_array, columns=columns, index=pgc)
        target_dir = os.path.dirname(target[2])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[2], 'w') as f:
            f.write('# all columns are normalized tendon force.\n')
            all_fse_df.to_csv(f)

        all_lMtilde_df = pd.DataFrame(all_lMtilde_array, columns=columns, 
            index=pgc)
        target_dir = os.path.dirname(target[3])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[3], 'w') as f:
            f.write('# all columns are normalized fiber length.\n')
            all_lMtilde_df.to_csv(f)

        all_vMtilde_df = pd.DataFrame(all_vMtilde_array, columns=columns, 
            index=pgc)
        target_dir = os.path.dirname(target[4])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[4], 'w') as f:
            f.write('# all columns are normalized fiber velocity.\n')
            all_vMtilde_df.to_csv(f)

        # all_lTtilde_df = pd.DataFrame(all_lTtilde_array, columns=columns, 
        #     index=pgc)
        # target_dir = os.path.dirname(target[4])
        # if not os.path.exists(target_dir):
        #     os.makedirs(target_dir)
        # with file(target[4], 'w') as f:
        #     f.write('# all columns are normalized tendon length.\n')
        #     all_lTtilde_df.to_csv(f) 
        # How to read this in: df.read_csv(..., index_col=0, header=[0, 1, 2],
        #                                  skiprows=1)


class TaskPlotMuscleData(osp.StudyTask):
    REGISTRY = []
    def __init__(self, study, agg_task, mod_name=None, 
                 agg_tasks_to_compare=None, mod_names_to_compare=None,  
                 subject=None):
        super(TaskPlotMuscleData, self).__init__(study)
        self.agg_task_targets = list(agg_task.targets)
        self.agg_tasks_to_compare_targets = None
        self.mod_names_to_compare = list()
        if agg_tasks_to_compare:
            self.agg_tasks_to_compare_targets = list()
            for agg_task_to_compare in agg_tasks_to_compare:
                self.agg_tasks_to_compare_targets.append(
                    list(agg_task_to_compare.targets))
        
            if mod_names_to_compare == None:
                for agg_task_to_compare in agg_tasks_to_compare:
                    self.mod_names_to_compare.append('other')
            else:
                self.mod_names_to_compare = mod_names_to_compare

        suffix = ''
        if mod_name == None:
            suffix += '_experiments'
            mod_name = 'experiment'
        else:
            suffix += '_%s' % mod_name
            for i in range(len(self.agg_tasks_to_compare_targets)):
                # Duplicate experiment tasks, since two mod tasks need to 
                # compare to one experiment task
                if 'experiment' in self.agg_tasks_to_compare_targets[i][0]:
                    self.agg_tasks_to_compare_targets[i] += \
                        self.agg_tasks_to_compare_targets[i]

            # print 'agg_task_targets:'
            # for target in self.agg_task_targets:
            #     print target
            # print ''
            # print 'agg_tasks_to_compare_targets:'
            # for compare_targets in self.agg_tasks_to_compare_targets:
            #     for target in compare_targets:
            #         print target
            # print ''
            
        self.costdir = ''
        if not (study.costFunction == 'Default'):
            suffix += '_%s' % study.costFunction
            self.costdir = study.costFunction 

        if subject == None:
            subjects = [s.name for s in study.subjects]
        else:
            suffix += '_%s' % subject
            subjects = [subject]

        def remove_prefix(text, prefix):
            if text.startswith(prefix):
                return text[len(prefix):]
            return text

        self.name = 'plot_%s' % remove_prefix(agg_task.name, 'aggregate_')
        self.doc = 'Plot muscle data for experiment and mod tasks'

        if agg_tasks_to_compare:
            for itable, table in enumerate(self.agg_task_targets):
                    
                # print 'table1: ', table
                # print 'table2: ', compare_tables[0]
                # print ''
                table_list = list()
                table_list.append(table)
                mod_name_list = list()
                mod_name_list.append(mod_name)
                for compare_targets, mod_name_to_compare in zip(
                    self.agg_tasks_to_compare_targets,
                    self.mod_names_to_compare):
                        table_list.append(compare_targets[itable])
                        mod_name_list.append(mod_name_to_compare)

                self.add_action(table_list,
                                [],
                                self.plot_muscle_data, 
                                mod_name_list)
        else: 
            for table in self.agg_task_targets:
                # print 'just table: ', table
                # print ''
                self.add_action([table],
                                [],
                                self.plot_muscle_data,
                                [mod_name])
 
    def plot_muscle_data(self, file_dep, target, mod_names):

        def plot_table(fig, table, label, color):

            df_all = pd.read_csv(table, index_col=0,
                    header=[0, 1, 2], skiprows=1)

            # Average over cycles.
            df_by_subj_musc = df_all.groupby(
                level=['subject', 'muscle'], axis=1).mean()
            df_mean = df_by_subj_musc.groupby(level=['muscle'], axis=1).mean()
            df_std = df_by_subj_musc.groupby(level=['muscle'], axis=1).std()

            pgc = df_mean.index
            muscles = self.study.muscle_names

            # print('\ndoing doing going doing this thing')
            
            # print(df_all.shape)
            # print(df_by_subj_musc.shape)
            # print(df_mean.shape)

            # print('TODO here: probably need to edit this to handle full model')

            # nice_act_names = {
            #         # 'glut_max2_r': 'glut. max.',
            #         # 'psoas_r': 'iliopsoas',
            #         # 'semimem_r': 'hamstrings',
            #         # 'rect_fem_r': 'rect. fem.',
            #        # 'bifemsh_r': 'bi. fem. s.h.',
            #         # 'vas_int_r': 'vasti',
            #         # 'med_gas_r': 'gastroc.',
            #         # 'soleus_r': 'soleus',
            #         # 'tib_ant_r': 'tib. ant.',

            #         'add_brev_r': 'add. long.',
            #         'add_long_r': 'add. long.',
            #         'add_mag3_r': 'add. mag. 3',
            #         'add_mag4_r': 'add. mag. 4',
            #         'add_mag2_r': 'add. mag. 2',
            #         'add_mag1_r': 'add. mag. 1',
            #         'bifemlh_r': 'bi. fem. l.h.',
            #         'bifemsh_r': 'bi. fem. s.h.',
            #         'ext_dig_r': 'ext. dig.',
            #         'ext_hal_r': 'ext. hal.',
            #         'flex_dig_r': 'flex. dig.',
            #         'flex_hal_r': 'flex. hal.',
            #         'lat_gas_r': 'lat. gastoc.',
            #         'med_gas_r': 'med. gastroc.',
            #         'glut_max1_r': 'glut. max. 1',
            #         'glut_max2_r': 'glut. max. 2',
            #         'glut_max3_r': 'glut. max. 3',
            #         'glut_med1_r': 'glut. med. 1',
            #         'glut_med2_r': 'glut. med. 2',
            #         'glut_med3_r': 'glut. med. 3',
            #         'glut_min1_r': 'glut. min. 1',
            #         'glut_min2_r': 'glut. min. 2',
            #         'glut_min3_r': 'glut. min. 3',
            #         'grac_r': 'gracilis',
            #         'iliacus_r': 'iliacus',
            #         'per_brev_r': 'per. brev.',
            #         'per_long_r': 'per. long.',
            #         'peri_r': 'peri',
            #         'psoas_r': 'iliopsoas',
            #         'rect_fem_r': 'rect. fem.',
            #         'sar_r': 'sartorius',
            #         'semimem_r': 'semimem.',
            #         'semiten_r': 'semiten.',
            #         'soleus_r': 'soleus',
            #         'tfl_r': 'tfl.',
            #         'tib_ant_r': 'tib. ant.',
            #         'tib_post_r': 'tib. post.',
            #         'vas_int_r': 'vas. int.',
            #         'vas_lat_r': 'vas. lat.',
            #         'vas_med_r': 'vas. med.'

            #         }
            ## have to adjust again for the new model 
            nice_act_names = {
                'addbrev_r':'add. brev.',
                'addlong_r':'add. long.',
                'addmagDist_r':'add. mag. 3',
                'addmagIsch_r':'add. mag. 4',
                'addmagMid_r':'add. mag. 2',
                'addmagProx_r':'add. mag. 1',
                'bflh_r':'bifem. lh.',
                'bfsh_r':'bifem. sh.',
                'edl_r':'ext. dig.',
                'ehl_r':'ext. hal.',
                'fdl_r':'flex. dig.',
                'fhl_r':'flex. hal.',
                'gaslat_r':'lat. gastroc',
                'gasmed_r':'med. gastroc',
                'glmax1_r':'glut. max. 1',
                'glmax2_r':'glut. max 2',
                'glmax3_r':'glut. max 3',
                'glmed1_r':'glut. med 1',
                'glmed2_r':'glut. med 2',
                'glmed3_r':'glut. med 3',
                'glmin1_r':'glut. min 1',
                'glmin2_r':'glut. min 2',
                'glmin3_r':'glut. min 3',
                'grac_r':'gracilis',
                'iliacus_r':'iliacus',
                'perbrev_r':'per. brev.',
                'perlong_r':'per. long.',
                'piri_r':'peri.',
                'psoas_r':'psoas',
                'recfem_r':'rec. fem.',
                'sart_r':'sartorius',
                'semimem_r':'semimem.',
                'semiten_r':'semiten',
                'soleus_r':'soleus',
                'tfl_r':'tfl.',
                'tibant_r':'tib. ant.',
                'tibpost_r':'tib. post.',
                'vasint_r':'vas. int.',
                'vaslat_r':'vas. lat.',
                'vasmed_r':'vas. med.',
                'net':'net',
                'active': 'active',
                'passive': 'passive',
                'add_brev_r':'add. brev.',
                'add_long_r':'add. long.',
                'add_mag3_r':'add. mag. 3',
                'add_mag4_r':'add. mag. 4',
                'add_mag2_r':'add. mag. 2',
                'add_mag1_r':'add. mag. 1',
                'bifemlh_r':'bifem. lh.',
                'bifemsh_r':'bifem. sh.',
                'ext_dig_r':'ext. dig.',
                'ext_hal_r':'ext. hal.',
                'flex_dig_r':'flex. dig.',
                'flex_hal_r':'flex. hal.',
                'lat_gas_r':'lat. gastroc',
                'med_gas_r':'med. gastroc',
                'glut_max1_r':'glut. max. 1',
                'glut_max2_r':'glut. max 2',
                'glut_max3_r':'glut. max 3',
                'glut_med1_r':'glut. med 1',
                'glut_med2_r':'glut. med 2',
                'glut_med3_r':'glut. med 3',
                'glut_min1_r':'glut. min 1',
                'glut_min2_r':'glut. min 2',
                'glut_min3_r':'glut. min 3',
                'grac_r':'gracilis',
                'iliacus_r':'iliacus',
                'per_brev_r':'per. brev.',
                'per_long_r':'per. long.',
                'peri_r':'peri.',
                'psoas_r':'psoas',
                'rect_fem_r':'rec. fem.',
                'sar_r':'sartorius',
                'semimem_r':'semimem.',
                'semiten_r':'semiten',
                'soleus_r':'soleus',
                'tfl_r':'tfl.',
                'tib_ant_r':'tib. ant.',
                'tib_post_r':'tib. post.',
                'vas_int_r':'vas. int.',
                'vas_lat_r':'vas. lat.',
                'vas_med_r':'vas. med.',
                
        }

            for imusc, musc_name in enumerate(muscles):
                side_len = np.ceil(np.sqrt(len(muscles)))
                ax = fig.add_subplot(side_len, side_len, imusc + 1)
                ax.axhline(color='k', linewidth=0.5, zorder=0)
                y_mean = df_mean[musc_name]
                y_std = df_std[musc_name]
                ax.plot(pgc, y_mean, color=color, linestyle='-', label=label)
                ax.fill_between(pgc, y_mean-y_std, y_mean+y_std,
                        color=color, alpha=0.3)
                if imusc == 0:
                    ax.legend(frameon=False, fontsize=6)
                ax.set_xlim(0, 100)
                if ('excitations' in table) or ('activations' in table):
                    ax.set_ylim(0, 1.0)
                elif ('norm_fiber_length' in table):
                    ax.set_ylim(0.4, 1.6)
                elif ('norm_fiber_velocity' in table):
                    ax.set_ylim(-0.7, 0.7)
                elif ('norm_tendon_force' in table):
                    ax.set_ylim(0, 1.5)
                elif ('norm_active_force' in table):
                    ax.set_ylim(0, 1.5)
                elif ('norm_passive_force' in table):
                    ax.set_ylim(0, 0.5)
                else:
                    ax.set_ylim(0, 1.0)
                ax.set_title(nice_act_names[musc_name])
                ax.set_xlabel('time (% gait cycle)')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')


        fig = pl.figure(figsize=(14, 14))
        import matplotlib.pyplot as plt 
        cmap = plt.get_cmap('bwr')
        dc = 2.0 / len(file_dep)

        for i, table in enumerate(file_dep):
            plot_table(fig, table, mod_names[i], cmap((i+1)*dc))
        fig.tight_layout()
        fig.savefig(file_dep[0].replace('.csv', '.pdf'))
        fig.savefig(file_dep[0].replace('.csv', '.png'), dpi=600)
        pl.close(fig)



## Adding in all the calibration tasks
class TaskScaleMuscleMaxIsometricForce(osp.SubjectTask):
    REGISTRY = []
    def __init__(self, subject):
        super(TaskScaleMuscleMaxIsometricForce, self).__init__(subject)
        self.subject = subject
        self.name = '%s_scale_max_force' % self.subject.name
        self.doc = 'Scale subject muscle Fmax parameters from Handsfield2014'
        self.generic_model_fpath = self.study.source_generic_model_fpath
        self.subject_model_fpath = os.path.join(self.subject.results_exp_path, 
            '%s.osim' % self.subject.name)
        self.scaled_param_model_fpath = os.path.join(
            self.subject.results_exp_path, 
            '%s_scaled_Fmax.osim' % self.subject.name)

        self.add_action([self.generic_model_fpath, self.subject_model_fpath],
                        [self.scaled_param_model_fpath],
                        self.scale_model_parameters)

    def scale_model_parameters(self, file_dep, target):
        """From Handsfields 2014 figure 5a and from Apoorva's muscle properties
       spreadsheet.
       
       v: volume fraction
       V: total volume
       F: max isometric force
       l: optimal fiber length

       F = v * sigma * V / l

       *_g: generic model.
       *_s: subject-specific model.

       F_g = v * sigma * V_g / l_g
       F_s = v * sigma * V_s / l_s

       F_s = (F_g * l_g / V_g) * V_s / l_s
           = F_g * (V_s / V_g) * (l_g / l_s)

        Author: Chris Dembia 
        Borrowed from mrsdeviceopt GitHub repo:
        https://github.com/chrisdembia/mrsdeviceopt          
       """

        print("\nMuscle force scaling: "
              "total muscle volume and optimal fiber length.")
        # print('\n' + str(self.subject.mass))
        def total_muscle_volume_regression(mass):
            return 91.0*mass + 588.0

        generic_TMV = total_muscle_volume_regression(75.337)
        subj_TMV = total_muscle_volume_regression(self.subject.mass)

        import opensim as osm
        generic_model = osm.Model(file_dep[0])
        # print file_dep[0]
        # print file_dep[1]
        subj_model = osm.Model(file_dep[1])


        generic_mset = generic_model.getMuscles()
        subj_mset = subj_model.getMuscles()


        for im in range(subj_mset.getSize()):
            muscle_name = subj_mset.get(im).getName()
            generic_muscle = generic_mset.get(muscle_name)
            subj_muscle = subj_mset.get(muscle_name)
            generic_OFL = generic_muscle.get_optimal_fiber_length()
            subj_OFL = subj_muscle.get_optimal_fiber_length()

            scale_factor = (subj_TMV / generic_TMV) * (generic_OFL / subj_OFL)
            print("Scaling '%s' muscle force by %f." % (muscle_name,
                scale_factor))

            generic_force = generic_muscle.getMaxIsometricForce()
            scaled_force = generic_force * scale_factor
            subj_muscle.setMaxIsometricForce(scaled_force)

        subj_model.printToXML(target[0])



class TaskCalibrateParametersSetup(osp.SetupTask):
    REGISTRY = []
    def __init__(self, trial, param_dict, cost_dict, passive_precalibrate=False,
            **kwargs):
        super(TaskCalibrateParametersSetup, self).__init__('calibrate', trial,
            **kwargs)
        self.doc = "Create a setup file for a parameter calibration tool."
        self.kinematics_file = os.path.join(self.trial.results_exp_path, 'ik',
                '%s_%s_ik_solution.mot' % (self.study.name, self.trial.id))
        self.rel_kinematics_file = os.path.relpath(self.kinematics_file,
                self.path)
        self.kinetics_file = os.path.join(self.trial.results_exp_path,
                'id', 'results', '%s_%s_id_solution.sto' % (self.study.name,
                    self.trial.id))
        self.rel_kinetics_file = os.path.relpath(self.kinetics_file,
                self.path)
        self.emg_file = os.path.join(self.trial.results_exp_path, 
                'expdata', 'emg.sto')
        self.rel_emg_file = os.path.relpath(self.emg_file, self.path)
        self.results_setup_fpath = os.path.join(self.path, 'setup.m')
        self.results_output_fpath = os.path.join(self.path, 
            '%s_%s_calibrate.mat' % (self.study.name, self.tricycle.id))

        self.param_dict = param_dict
        self.cost_dict = cost_dict
        self.passive_precalibrate = passive_precalibrate

        # Fill out setup.m template and write to results directory
        self.create_setup_action()

    def create_setup_action(self): 
        self.add_action(
                    ['templates/%s/setup.m' % self.tool],
                    [self.results_setup_fpath],
                    self.fill_setup_template,  
                    init_time=self.init_time,
                    final_time=self.final_time,      
                    )


    def fill_setup_template(self, file_dep, target,
                            init_time=None, final_time=None):

        print('\nup here as well\n')
        print(file_dep)

        with open(file_dep[0]) as ft:
            content = ft.read()

            possible_params = ['optimal_fiber_length', 'tendon_slack_length',
                               'pennation_angle', 'muscle_strain']
            paramstr = ''
            for param in possible_params:
                if param in self.param_dict:
                    paramstr += param + ' = true;\n'
                else:
                    paramstr += param + ' = false;\n'

            possible_costs = ['emg']
            coststr = ''
            for cost in possible_costs:
                if cost in self.cost_dict:
                    coststr += cost + ' = true;\n'
                else:
                    coststr += cost + ' = false;\n'


            pass_cal = ''
            if self.passive_precalibrate:
                pass_cal = 'Misc.passive_precalibrate = true;\n'

            content = content.replace('Misc = struct();',
                'Misc = struct();\n' + paramstr + coststr + pass_cal + '\n')

            content = content.replace('@STUDYNAME@', self.study.name)
            content = content.replace('@NAME@', self.tricycle.id)
            # TODO should this be an RRA-adjusted model?
            content = content.replace('@MODEL@', os.path.relpath(
                self.subject.scaled_model_fpath, self.path))
            content = content.replace('@REL_PATH_TO_TOOL@', os.path.relpath(
                self.study.config['optctrlmuscle_path'], self.path))
            # TODO provide slop on either side? start before the cycle_start?
            # end after the cycle_end?
            content = content.replace('@INIT_TIME@',
                    '%.5f' % init_time)
            content = content.replace('@FINAL_TIME@', 
                    '%.5f' % final_time)
            content = content.replace('@IK_SOLUTION@',
                    self.rel_kinematics_file)
            content = content.replace('@ID_SOLUTION@',
                    self.rel_kinetics_file)
            content = content.replace('@SIDE@',
                    self.trial.primary_leg[0])
            content = content.replace('@EMG_PATH@', self.rel_emg_file)
            if 'optimal_fiber_length' in self.param_dict:
                content = content.replace('@lMo_MUSCLES@',
                        ','.join(self.param_dict['optimal_fiber_length']))
            if 'tendon_slack_length' in self.param_dict:
                content = content.replace('@lTs_MUSCLES@',
                        ','.join(self.param_dict['tendon_slack_length']))
            if 'pennation_angle' in self.param_dict:
                content = content.replace('@alf_MUSCLES@',
                        ','.join(self.param_dict['pennation_angle']))
            if 'muscle_strain' in self.param_dict:
                content = content.replace('@e0_MUSCLES@',
                        ','.join(self.param_dict['muscle_strain']))
            if 'emg' in self.cost_dict:
                content = content.replace('@emg_MUSCLES@',
                        ','.join(self.cost_dict['emg']))

        print('\nhere we go again\n')
        print(target[0])

        with open(target[0], 'w') as f:
            f.write(content)

class TaskCalibrateParameters(osp.ToolTask):
    REGISTRY = []
    def __init__(self, trial, calibrate_setup_task, **kwargs):
        super(TaskCalibrateParameters, self).__init__(calibrate_setup_task, 
            trial, opensim=False, **kwargs)
        self.doc = "Run parameter calibration tool via DeGroote MRS solver."
        self.results_setup_fpath = calibrate_setup_task.results_setup_fpath
        self.results_output_fpath = calibrate_setup_task.results_output_fpath

        self.file_dep += [
                self.results_setup_fpath,
                self.subject.scaled_model_fpath,
                calibrate_setup_task.kinematics_file,
                calibrate_setup_task.kinetics_file,
                calibrate_setup_task.emg_file,
                ]

        self.actions += [
                self.run_parameter_calibration,
                self.delete_muscle_analysis_results,
                ]

        self.targets += [
                self.results_output_fpath
                ]

    def run_parameter_calibration(self):
        with util.working_directory(self.path):
            # On Mac, CmdAction was causing MATLAB ipopt with GPOPS output to
            # not display properly.

            status = os.system('matlab %s -logfile matlab_log.txt -wait -r "try, '
                    "run('%s'); disp('SUCCESS'); "
                    'catch ME; disp(getReport(ME)); exit(2), end, exit(0);"\n'
                    % ('-automation' if os.name == 'nt' else '',
                        self.results_setup_fpath)
                    )
            if status != 0:
                # print 'Non-zero exist status. Continuing....'
                raise Exception('Non-zero exit status.')

            # Wait until output mat file exists to finish the action
            import time
            while True:
                time.sleep(3.0)

                mat_exists = os.path.isfile(self.results_output_fpath)
                if mat_exists:
                    break

    def delete_muscle_analysis_results(self):
        if os.path.exists(os.path.join(self.path, 'results')):
            import shutil
            shutil.rmtree(os.path.join(self.path, 'results'))


def get_muscle_parameters_as_list(fpath):

    import h5py
    hdf_output = h5py.File(fpath, 'r')
    params = hdf_output['OptInfo']['paramCal']

    lMo = list()
    lTs = list()
    alf = list()
    e0 = list()
    musc_names = list()
    for musc in params:
        musc_names.append(musc)
        for param in params[musc]:
            if param == 'optimal_fiber_length':
                lMo.append(params[musc][param][0][0])
            elif param == 'tendon_slack_length':
                lTs.append(params[musc][param][0][0])
            elif param == 'pennation_angle':
                alf.append(params[musc][param][0][0])
            elif param == 'muscle_strain':
                e0.append(params[musc][param][0][0])

    return lMo, lTs, alf, e0, musc_names

def get_muscle_parameters_as_dict(fpath):

    import h5py
    hdf_output = h5py.File(fpath, 'r')
    params = hdf_output['OptInfo']['paramCal']

    lMo = dict()
    lTs = dict()
    alf = dict()
    e0 = dict()
    musc_names = list()
    for musc in params:
        musc_names.append(musc)
        for param in params[musc]:
            if param == 'optimal_fiber_length':
                lMo[musc] = params[musc][param][0][0]
            elif param == 'tendon_slack_length':
                lTs[musc] = params[musc][param][0][0]
            elif param == 'pennation_angle':
                alf[musc] = params[musc][param][0][0]
            elif param == 'muscle_strain':
                e0[musc] = params[musc][param][0][0]

    return lMo, lTs, alf, e0, musc_names



class TaskCalibrateParametersPost(osp.PostTask):
    REGISTRY = []
    def __init__(self, trial, calibrate_setup_task, **kwargs):
        super(TaskCalibrateParametersPost, self).__init__(calibrate_setup_task,
            trial, **kwargs)
        self.doc = 'Postprocessing of parameter calibration results.'
        self.setup_task = calibrate_setup_task
        self.results_output_fpath = self.setup_task.results_output_fpath

        self.emg_fpath = os.path.join(trial.results_exp_path, 'expdata', 
            'emg_with_headers.sto')

        self.add_action([self.emg_fpath,
                         self.results_output_fpath],
                        [os.path.join(self.path, 'muscle_activity'),
                         os.path.join(self.path, 'reserve_activity.pdf')],
                        self.plot_muscle_and_reserve_activity)

        self.add_action([self.results_output_fpath],
                        [os.path.join(self.path, 'optimal_fiber_length.pdf'),
                         os.path.join(self.path, 'tendon_slack_length.pdf'),
                         os.path.join(self.path, 'pennation_angle.pdf'),
                         os.path.join(self.path, 'muscle_strain.pdf')],
                        self.plot_muscle_parameters)

    def plot_muscle_and_reserve_activity(self, file_dep, target):

        emg = util.storage2numpy(file_dep[0])
        time = emg['time']

        def min_index(vals):
            idx, val = min(enumerate(vals), key=lambda p: p[1])
            return idx

        start_idx = min_index(abs(time-self.setup_task.cycle.start))
        end_idx = min_index(abs(time-self.setup_task.cycle.end))

        # Load mat file fields
        muscle_names = util.hdf2list(file_dep[1], 'MuscleNames', type=str)
        df_exc = util.hdf2pandas(file_dep[1], 'MExcitation', labels=muscle_names)
        df_act = util.hdf2pandas(file_dep[1], 'MActivation', labels=muscle_names)
        dof_names = util.hdf2list(file_dep[1], 'DatStore/DOFNames', 
            type=str)

        pgc_emg = np.linspace(0, 100, len(time[start_idx:end_idx]))
        pgc_exc = np.linspace(0, 100, len(df_exc.index))
        pgc_act = np.linspace(0, 100, len(df_act.index))

        muscles = self.study.muscle_names
        fig = pl.figure(figsize=(12, 12))
        nice_act_names = {
                'glut_max2_r': 'glut. max.',
                'psoas_r': 'iliopsoas',
                'semimem_r': 'hamstrings',
                'rect_fem_r': 'rect. fem.',
                'bifemsh_r': 'bi. fem. s.h.',
                'vas_int_r': 'vasti',
                'med_gas_r': 'gastroc.',
                'soleus_r': 'soleus',
                'tib_ant_r': 'tib. ant.',
                }

        emg_map = {
                'med_gas_r': 'gasmed_r',
                'glut_max2_r': 'glmax2_r',
                'rect_fem_r': 'recfem_r',
                'semimem_r': 'semimem_r',
                'soleus_r': 'soleus_r',
                'tib_ant_r': 'tibant_r',
                'vas_int_r': 'vasmed_r', 
        }

        emg_muscles = ['bflh_r', 'gaslat_r', 'gasmed_r', 'glmax1_r', 'glmax2_r',
                       'glmax3_r', 'glmed1_r', 'glmed2_r', 'glmed3_r', 
                       'recfem_r', 'semimem_r', 'semiten_r', 'soleus_r',
                       'tibant_r', 'vaslat_r', 'vasmed_r']

        for imusc, musc_name in enumerate(muscles):
            side_len = np.ceil(np.sqrt(len(muscles)))
            ax = fig.add_subplot(side_len, side_len, imusc + 1)
            ax.axhline(color='k', linewidth=0.5, zorder=0)
            y_exc = df_exc[musc_name]
            y_act = df_act[musc_name]
            exc_plot, = ax.plot(pgc_exc, y_exc, color='blue', 
                linestyle='--')
            act_plot, = ax.plot(pgc_act, y_act, color='red', 
                linestyle='--')
            handles = [exc_plot, act_plot   ]
            labels = ['%s exc.' % nice_act_names[musc_name],
                      '%s act.' % nice_act_names[musc_name]]
            ax.legend(handles, labels)
            
            if emg_map.get(musc_name):
                y_emg = emg[emg_map[musc_name]]
                ax.plot(pgc_emg, y_emg[start_idx:end_idx], color='black', 
                    linestyle='-')

            # ax.legend(frameon=False, fontsize=6)
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 1.0)
            ax.set_title(nice_act_names[musc_name])
            ax.set_xlabel('time (% gait cycle)')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        fig.tight_layout()
        fig.savefig(target[0]+'.pdf')
        fig.savefig(target[0]+'.png', dpi=600)
        pl.close(fig)

        # Plot reserve activity
        df_res = util.hdf2pandas(file_dep[1],'RActivation', labels=dof_names)
        pp.plot_reserve_activity(target[1], df_res)

    def plot_muscle_parameters(self, file_dep, target):

        lMo, lTs, alf, e0, musc_names = get_muscle_parameters_as_dict(
            file_dep[0])

        # Plot muscle parameters
        def param_barplot(param_name, param, output_fpath):

            fig = pl.figure(figsize=(11,6))
            ax = fig.add_subplot(1,1,1)
            musc_names = list(param.keys())
            pos = np.arange(len(musc_names))

            param_to_plot = [(val - 1.0)*100 for val in param.values()]
            bar = ax.bar(pos, param_to_plot, color='green')
            ax.set_xticks(pos)
            ax.set_xticklabels(musc_names, fontsize=10)
            ax.set_ylabel('Percent Change in ' + param_name, fontsize=12)
            ax.set_ylim([-30, 30])
            ax.set_yticks(np.linspace(-30, 30, 13))
            ax.grid(which='both', axis='both', linestyle='--')
            ax.set_axisbelow(True)

            fig.tight_layout()
            fig.savefig(output_fpath)
            pl.close(fig)

        if lMo:
            param_barplot('Optimal Fiber Length', lMo, target[0])

        if lTs:
            param_barplot('Tendon Slack Length', lTs, target[1])

        if alf:
            param_barplot('Pennation Angle', alf, target[2])

        if e0:
            param_barplot('Muscle Strain', e0, target[3])



def construct_multiindex_tuples_for_subject(subject, conditions, 
    muscles=None, cycles_to_exclude=None):
    ''' Construct multiindex tuples and list of cycles for DataFrame indexing.
    '''
    
    multiindex_tuples = list()
    cycles = list()

    for cond_name in conditions:
        cond = subject.get_condition(cond_name)
        if not cond: continue
        # We know there is only one overground trial, but perhaps it
        # has not yet been added for this subject.
        assert len(cond.trials) <= 1
        if len(cond.trials) == 1:
            trial = cond.trials[0]
            for cycle in trial.cycles:
                if cycle.name in cycles_to_exclude: continue
                cycles.append(cycle)
                if not muscles:
                    multiindex_tuples.append((
                        cycle.condition.name,
                        # This must be the full ID, not just the cycle
                        # name, because 'cycle01' from subject 1 has
                        # nothing to do with 'cycle01' from subject 2
                        # (whereas the 'walk2' condition for subject 1 is
                        # related to 'walk2' for subject 2).
                        cycle.id))
                if muscles:
                    for mname in muscles:
                        multiindex_tuples.append((
                            cycle.condition.name,
                            cycle.id,
                            mname))

    return multiindex_tuples, cycles




class TaskAggregateMuscleParameters(osp.SubjectTask):
    """Aggregate calibrated muscle parameters for a given subject."""
    REGISTRY = []
    def __init__(self, subject, param_dict, 
            conditions=['walk1','walk2','walk3'], 
            cycles_to_exclude=['cycle03']):
        super(TaskAggregateMuscleParameters, self).__init__(subject)
        self.name = '%s_aggregate_muscle_parameters' % self.subject.name
        self.csv_fpaths = list()
        self.csv_params = list()
        self.param_dict = param_dict
        self.conditions = conditions
        self.cycles_to_exclude = cycles_to_exclude
        
        for param in param_dict:
            muscle_names = param_dict[param]
            multiindex_tuples, cycles = construct_multiindex_tuples_for_subject( 
                self.subject, conditions, muscle_names, self.cycles_to_exclude)

            # Prepare for processing simulations of experiments.
            deps = list()
            for cycle in cycles:
                if cycle.name in cycles_to_exclude: continue
                deps.append(os.path.join(
                        cycle.trial.results_exp_path, 'calibrate', cycle.name,
                        '%s_%s_calibrate.mat' % (self.study.name, cycle.id))
                        )

            csv_path = os.path.join(self.subject.results_exp_path, 
                '%s_agg.csv' % param)
            self.csv_params.append(param)
            self.csv_fpaths.append(csv_path)
            self.add_action(deps,
                    [csv_path],
                    self.aggregate_muscle_parameters, param, multiindex_tuples)

    def aggregate_muscle_parameters(self, file_dep, target, param,
            multiindex_tuples):
        from collections import OrderedDict
        muscle_params = OrderedDict()
        for ifile, fpath in enumerate(file_dep):
            lMo, lTs, alf, e0, musc_names = get_muscle_parameters_as_dict(fpath)
            if not param in muscle_params:
                muscle_params[param] = list()
            for musc in self.param_dict[param]:
                if param == 'optimal_fiber_length':
                    muscle_params[param].append(lMo[musc])
                elif param == 'tendon_slack_length':
                    muscle_params[param].append(lTs[musc])
                elif param == 'pennation_angle':
                    muscle_params[param].append(alf[musc])
                elif param == 'muscle_strain':
                    muscle_params[param].append(e0[musc])
       
        # http://pandas.pydata.org/pandas-docs/stable/advanced.html#advanced-hierarchical
        index = pd.MultiIndex.from_tuples(multiindex_tuples,
                names=['condition', 'cycle', 'muscle'])

        df = pd.DataFrame(muscle_params, index=index)

        target_dir = os.path.dirname(target[0])
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        with file(target[0], 'w') as f:
            f.write('# columns contain calibrated muscle parameters for ' 
                    '%s \n' % self.subject.name)
            df.to_csv(f)



class TaskPlotMuscleParameters(osp.SubjectTask):
    REGISTRY = []
    def __init__(self, subject, agg_task, cycles_to_exclude=None, **kwargs):
        super(TaskPlotMuscleParameters, self).__init__(subject)
        self.name = '%s_plot_muscle_parameters' % subject.name
        self.doc = 'Plot aggregated muscle parameters from calibration.'
        self.agg_task = agg_task
        self.cycles_to_exclude = cycles_to_exclude

        for csv_fpath, param in zip(agg_task.csv_fpaths, agg_task.csv_params):
            self.add_action([csv_fpath],
                            [os.path.join(self.subject.results_exp_path, param)],
                            self.plot_muscle_parameters, param)

    def plot_muscle_parameters(self, file_dep, target, param):

        # Process muscle parameters
        df = pd.read_csv(file_dep[0], index_col=[0, 1, 2], skiprows=1)
        if self.cycles_to_exclude:
            for cycle in self.cycles_to_exclude:
                for cond in self.agg_task.conditions:
                    cycle_id = '%s_%s_%s' % (self.subject.name, cond, cycle)
                    df.drop(cycle_id, level='cycle', inplace=True)

        df_mean = df.mean(level='muscle')
        df_std = df.std(level='muscle')

        def param_barplot(param_name, param_mean, param_std, musc_names, 
                output_fpath):

            fig = pl.figure(figsize=(11,6))
            ax = fig.add_subplot(1,1,1)
            pos = np.arange(len(musc_names))

            param_mean_to_plot = [(val - 1.0)*100 for val in param_mean]
            param_std_to_plot = [val*100 for val in param_std]
            bar = ax.bar(pos, param_mean_to_plot, yerr=param_std_to_plot,
                color='green')
            ax.set_xticks(pos)
            ax.set_xticklabels(musc_names, fontsize=10)
            ax.set_ylabel('Percent Change in ' + param_name, fontsize=12)
            ax.set_ylim([-30, 30])
            ax.set_yticks(np.linspace(-30, 30, 13))
            ax.grid(which='both', axis='both', linestyle='--')
            ax.set_axisbelow(True)

            fig.tight_layout()
            fig.savefig(output_fpath)
            pl.close(fig)

        if param == 'optimal_fiber_length':
            lMo_mean = df_mean['optimal_fiber_length']
            lMo_std = df_std['optimal_fiber_length']
            param_barplot('Optimal Fiber Length', lMo_mean, lMo_std, 
                lMo_mean.index, target[0]+'.pdf')
            lMo_mean.to_csv(target[0]+'.csv')

        if param == 'tendon_slack_length':
            lTs_mean = df_mean['tendon_slack_length']
            lTs_std = df_std['tendon_slack_length']
            param_barplot('Tendon Slack Length', lTs_mean, lTs_std, 
                lTs_mean.index, target[0]+'.pdf')
            lTs_mean.to_csv(target[0]+'.csv')

        if param == 'pennation_angle':
            alf_mean = df_mean['pennation_angle']
            alf_std = df_std['pennation_angle']
            param_barplot('Pennation Angle', alf_mean, alf_std, 
                alf_mean.index, target[0]+'.pdf')
            alf_mean.to_csv(target[0]+'.csv')

        if param == 'muscle_strain':
            e0_mean = df_mean['muscle_strain']
            e0_std = df_std['muscle_strain']
            param_barplot('Muscle Strain', e0_mean, e0_std, 
                e0_mean.index, target[0]+'.pdf')
            e0_mean.to_csv(target[0]+'.csv')




class TaskCalibrateParametersMultiPhaseSetup(osp.SetupTask):
    REGISTRY = []
    # TODO: the specific trial is ignored in this task, since this is 
    # really a subject task. 
    def __init__(self, trial, conditions, cycles, param_dict, cost_dict, 
            passive_precalibrate=False, **kwargs):
        super(TaskCalibrateParametersMultiPhaseSetup, self).__init__(
            'calibrate_multiphase', trial, **kwargs)

        self.doc = "Create a setup file for a mulit-phase parameter " \
                   "calibration tool."
        self.path = os.path.join(self.subject.results_exp_path, 'calibrate')
        self.name = '%s_calibrate_multiphase_setup' % self.subject.name

        self.kinematics_files = list()
        self.rel_kinematics_files = list()
        self.kinetics_files = list()
        self.rel_kinetics_files = list()
        self.emg_files = list()
        self.rel_emg_files = list()
        self.init_times = list()
        self.final_times = list()
        self.cycle_ids = list()
        self.trial_results_exp_paths = list()

        for cond in self.subject.conditions:
            if not (cond.name in conditions): continue
            # Assume only trial per condition
            trial = cond.trials[0]
            for cycle in trial.cycles:
                if not(cycle.name in cycles): continue

                kinematics_file = os.path.join(trial.results_exp_path, 'ik', 
                    '%s_%s_ik_solution.mot' % (self.study.name, trial.id))
                self.kinematics_files.append(kinematics_file)
                self.rel_kinematics_files.append(
                        os.path.relpath(kinematics_file, self.path)
                    )

                kinetics_file = os.path.join(trial.results_exp_path, 'id', 
                    'results', '%s_%s_id_solution.sto' % (self.study.name, 
                        trial.id))
                self.kinetics_files.append(kinetics_file)
                self.rel_kinetics_files.append(
                        os.path.relpath(kinetics_file, self.path)
                    )

                emg_file = os.path.join(trial.results_exp_path, 'expdata', 
                    'emg.sto')
                self.emg_files.append(emg_file)
                self.rel_emg_files.append(
                        os.path.relpath(emg_file, self.path)
                    )
                self.init_times.append(cycle.start)
                self.final_times.append(cycle.end)

                self.cycle_ids.append(cycle.id)
                self.trial_results_exp_paths.append(trial.results_exp_path)

        self.results_setup_fpath = os.path.join(self.path, 'setup.m')
        self.results_output_fpath = os.path.join(self.path, 
            '%s_%s_calibrate.mat' % (self.study.name, self.subject.name))
        self.param_dict = param_dict
        self.cost_dict = cost_dict
        self.passive_precalibrate = passive_precalibrate

        # Fill out setup.m template and write to results directory
        self.create_setup_action()

    def create_setup_action(self): 
        if not os.path.exists(self.path): os.makedirs(self.path)
        self.add_action(
                    ['templates/%s/setup.m' % self.tool],
                    [self.results_setup_fpath],
                    self.fill_setup_template,       
                    )

    def fill_setup_template(self, file_dep, target):
        with open(file_dep[0]) as ft:
            content = ft.read()

            possible_params = ['optimal_fiber_length', 'tendon_slack_length',
                               'pennation_angle', 'muscle_strain']
            paramstr = ''
            for param in possible_params:
                if param in self.param_dict:
                    paramstr += param + ' = true;\n'
                else:
                    paramstr += param + ' = false;\n'

            possible_costs = ['emg']
            coststr = ''
            for cost in possible_costs:
                if cost in self.cost_dict:
                    coststr += cost + ' = true;\n'
                else:
                    coststr += cost + ' = false;\n'


            pass_cal = ''
            if self.passive_precalibrate:
                pass_cal = 'Misc.passive_precalibrate = true;\n'

            content = content.replace('Misc = struct();',
                'Misc = struct();\n' + paramstr + coststr + pass_cal + '\n')

            content = content.replace('@STUDYNAME@', self.study.name)
            content = content.replace('@NAME@', self.subject.name)
            # TODO should this be an RRA-adjusted model?
            content = content.replace('@MODEL@', os.path.relpath(
                self.subject.scaled_model_fpath, self.path))
            content = content.replace('@REL_PATH_TO_TOOL@', os.path.relpath(
                self.study.config['optctrlmuscle_path'], self.path))
            # TODO provide slop on either side? start before the cycle_start?
            # end after the cycle_end?
            content = content.replace('@INIT_TIMES@',
                    ','.join('%0.5f' % x for x in self.init_times))
            content = content.replace('@FINAL_TIMES@', 
                    ','.join('%0.5f' % x for x in self.final_times))
            content = content.replace('@IK_SOLUTIONS@',
                    ','.join(self.rel_kinematics_files))
            content = content.replace('@ID_SOLUTIONS@',
                    ','.join(self.rel_kinetics_files))
            content = content.replace('@SIDE@',
                    self.trial.primary_leg[0])
            content = content.replace('@EMG_PATHS@', 
                    ','.join(self.rel_emg_files))
            content = content.replace('@CYCLE_IDS@',
                    ','.join(self.cycle_ids))
            if 'optimal_fiber_length' in self.param_dict:
                content = content.replace('@lMo_MUSCLES@',
                        ','.join(self.param_dict['optimal_fiber_length']))
            if 'tendon_slack_length' in self.param_dict:
                content = content.replace('@lTs_MUSCLES@',
                        ','.join(self.param_dict['tendon_slack_length']))
            if 'pennation_angle' in self.param_dict:
                content = content.replace('@alf_MUSCLES@',
                        ','.join(self.param_dict['pennation_angle']))
            if 'muscle_strain' in self.param_dict:
                content = content.replace('@e0_MUSCLES@',
                        ','.join(self.param_dict['muscle_strain']))
            if 'emg' in self.cost_dict:
                content = content.replace('@emg_MUSCLES@',
                        ','.join(self.cost_dict['emg']))

        with open(target[0], 'w') as f:
            f.write(content)





class TaskCalibrateParametersMultiPhase(osp.ToolTask):
    REGISTRY = []
    def __init__(self, trial, calibrate_setup_task, **kwargs):
        super(TaskCalibrateParametersMultiPhase, self).__init__(
            calibrate_setup_task, trial, opensim=False, **kwargs)
        self.doc = "Run multi-phase parameter calibration tool via DeGroote " \
                   "MRS solver."
        self.name = '%s_calibrate_multiphase' % self.subject.name
        self.results_setup_fpath = calibrate_setup_task.results_setup_fpath
        self.results_output_fpath = calibrate_setup_task.results_output_fpath

        self.file_dep += [
                self.results_setup_fpath,
                self.subject.scaled_model_fpath,
                ] + calibrate_setup_task.kinematics_files \
                + calibrate_setup_task.kinetics_files \
                + calibrate_setup_task.emg_files

        self.actions += [
                self.run_parameter_calibration,
                # self.delete_muscle_analysis_results,
                ]

        self.targets += [
                self.results_output_fpath
                ]

    def run_parameter_calibration(self):
        with util.working_directory(self.path):
            # On Mac, CmdAction was causing MATLAB ipopt with GPOPS output to
            # not display properly.

            status = os.system('matlab %s -logfile matlab_log.txt -wait -r "try, '
                    "run('%s'); disp('SUCCESS'); "
                    'catch ME; disp(getReport(ME)); exit(2), end, exit(0);"\n'
                    % ('-automation' if os.name == 'nt' else '',
                        self.results_setup_fpath)
                    )
            if status != 0:
                # print 'Non-zero exist status. Continuing....'
                raise Exception('Non-zero exit status.')

            # Wait until output mat file exists to finish the action
            import time
            while True:
                time.sleep(3.0)

                mat_exists = os.path.isfile(self.results_output_fpath)
                if mat_exists:
                    break

    def delete_muscle_analysis_results(self):
        if os.path.exists(os.path.join(self.path, 'results')):
            import shutil
            shutil.rmtree(os.path.join(self.path, 'results'))




class TaskCalibrateParametersMultiPhasePost(osp.PostTask):
    REGISTRY = []
    def __init__(self, trial, calibrate_setup_task, **kwargs):
        super(TaskCalibrateParametersMultiPhasePost, self).__init__(
            calibrate_setup_task, trial, **kwargs)
        self.doc = 'Postprocessing of multi-phase parameter calibration results.'
        self.name = '%s_calibrate_multiphase_post' % self.subject.name
        self.setup_task = calibrate_setup_task
        self.results_output_fpath = self.setup_task.results_output_fpath

        numPhases = len(calibrate_setup_task.kinematics_files)
        for i in range(numPhases):
            trial_results_exp_path = \
                calibrate_setup_task.trial_results_exp_paths[i]
            cycle_id = calibrate_setup_task.cycle_ids[i]
            init_time = calibrate_setup_task.init_times[i]
            final_time = calibrate_setup_task.final_times[i]

            emg_fpath = os.path.join(trial_results_exp_path, 'expdata', 
                'emg_with_headers.sto')

            self.add_action([emg_fpath, self.results_output_fpath],
                [os.path.join(self.path, '%s_muscle_activity' % cycle_id),
                 os.path.join(self.path, '%s_reserve_activity.pdf' % cycle_id)],
                self.plot_muscle_and_reserve_activity,
                cycle_id, init_time, final_time)

        self.add_action([self.results_output_fpath],
            [os.path.join(self.path, 'optimal_fiber_length.pdf'),
             os.path.join(self.path, 'tendon_slack_length.pdf'),
             os.path.join(self.path, 'pennation_angle.pdf'),
             os.path.join(self.path, 'muscle_strain.pdf')],
            self.plot_muscle_parameters)

    def plot_muscle_and_reserve_activity(self, file_dep, target, cycle_id, 
        init_time, final_time):

        emg = util.storage2numpy(file_dep[0])
        time = emg['time']

        def min_index(vals):
            idx, val = min(enumerate(vals), key=lambda p: p[1])
            return idx

        start_idx = min_index(abs(time-init_time))
        end_idx = min_index(abs(time-final_time))

        # Load mat file fields
        muscle_names = util.hdf2list(file_dep[1], 'MuscleNames', type=str)
        df_exc = util.hdf2pandas(file_dep[1], 'MExcitation/%s' % cycle_id, 
            labels=muscle_names)
        df_act = util.hdf2pandas(file_dep[1], 'MActivation/%s' % cycle_id, 
            labels=muscle_names)
        dof_names = util.hdf2list(file_dep[1], 'DOFNames', 
            type=str)

        pgc_emg = np.linspace(0, 100, len(time[start_idx:end_idx]))
        pgc_exc = np.linspace(0, 100, len(df_exc.index))
        pgc_act = np.linspace(0, 100, len(df_act.index))

        muscles = self.study.muscle_names
        fig = pl.figure(figsize=(12, 12))
        nice_act_names = {
                'glut_max2_r': 'glut. max.',
                'psoas_r': 'iliopsoas',
                'semimem_r': 'hamstrings',
                'rect_fem_r': 'rect. fem.',
                'bifemsh_r': 'bi. fem. s.h.',
                'vas_int_r': 'vasti',
                'med_gas_r': 'gastroc.',
                'soleus_r': 'soleus',
                'tib_ant_r': 'tib. ant.',
                }

        emg_map = {
                'med_gas_r': 'gasmed_r',
                'glut_max2_r': 'glmax2_r',
                'rect_fem_r': 'recfem_r',
                'semimem_r': 'semimem_r',
                'soleus_r': 'soleus_r',
                'tib_ant_r': 'tibant_r',
                'vas_int_r': 'vasmed_r', 
        }

        emg_muscles = ['bflh_r', 'gaslat_r', 'gasmed_r', 'glmax1_r', 'glmax2_r',
                       'glmax3_r', 'glmed1_r', 'glmed2_r', 'glmed3_r', 
                       'recfem_r', 'semimem_r', 'semiten_r', 'soleus_r',
                       'tibant_r', 'vaslat_r', 'vasmed_r']

        for imusc, musc_name in enumerate(muscles):
            side_len = np.ceil(np.sqrt(len(muscles)))
            ax = fig.add_subplot(side_len, side_len, imusc + 1)
            ax.axhline(color='k', linewidth=0.5, zorder=0)
            y_exc = df_exc[musc_name]
            y_act = df_act[musc_name]
            exc_plot, = ax.plot(pgc_exc, y_exc, color='blue', 
                linestyle='--')
            act_plot, = ax.plot(pgc_act, y_act, color='red', 
                linestyle='--')
            handles = [exc_plot, act_plot   ]
            labels = ['%s exc.' % nice_act_names[musc_name],
                      '%s act.' % nice_act_names[musc_name]]
            ax.legend(handles, labels)
            
            if emg_map.get(musc_name):
                y_emg = emg[emg_map[musc_name]]
                ax.plot(pgc_emg, y_emg[start_idx:end_idx], color='black', 
                    linestyle='-')

            # ax.legend(frameon=False, fontsize=6)
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 1.0)
            ax.set_title(nice_act_names[musc_name])
            ax.set_xlabel('time (% gait cycle)')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        fig.tight_layout()
        fig.savefig(target[0]+'.pdf')
        fig.savefig(target[0]+'.png', dpi=600)
        pl.close(fig)

        # Plot reserve activity
        df_res = util.hdf2pandas(file_dep[1],'RActivation/%s' % cycle_id, 
            labels=dof_names)
        pp.plot_reserve_activity(target[1], df_res)

    def plot_muscle_parameters(self, file_dep, target):

        lMo, lTs, alf, e0, musc_names = get_muscle_parameters_as_dict(
            file_dep[0])

        # Plot muscle parameters
        def param_barplot(param_name, param, output_fpath):

            fig = pl.figure(figsize=(11,6))
            ax = fig.add_subplot(1,1,1)
            musc_names = list(param.keys())
            pos = np.arange(len(musc_names))

            param_to_plot = [(val - 1.0)*100 for val in param.values()]
            bar = ax.bar(pos, param_to_plot, color='green')
            ax.set_xticks(pos)
            ax.set_xticklabels(musc_names, fontsize=10)
            ax.set_ylabel('Percent Change in ' + param_name, fontsize=12)
            ax.set_ylim([-30, 30])
            ax.set_yticks(np.linspace(-30, 30, 13))
            ax.grid(which='both', axis='both', linestyle='--')
            ax.set_axisbelow(True)

            fig.tight_layout()
            fig.savefig(output_fpath)
            pl.close(fig)

        if lMo:
            param_barplot('Optimal Fiber Length', lMo, target[0])
            fpath = os.path.join(self.subject.results_exp_path, 
                    os.path.basename(target[0]))
            param_barplot('Optimal Fiber Length', lMo, fpath)
            df_lMo = pd.DataFrame.from_dict(lMo, orient='index')
            df_lMo.to_csv(fpath.replace('.pdf', '.csv'), header=False)

        if lTs:
            param_barplot('Tendon Slack Length', lTs, target[1])
            fpath = os.path.join(self.subject.results_exp_path, 
                    os.path.basename(target[1]))
            param_barplot('Tendon Slack Length', lTs, fpath)
            df_lTs = pd.DataFrame.from_dict(lTs, orient='index')
            df_lTs.to_csv(fpath.replace('.pdf', '.csv'), header=False)

        if alf:
            param_barplot('Pennation Angle', alf, target[2])
            fpath = os.path.join(self.subject.results_exp_path, 
                    os.path.basename(target[2]))
            param_barplot('Pennation Angle', alf, fpath)
            df_alf = pd.DataFrame.from_dict(alf, orient='index')
            df_alf.to_csv(fpath.replace('.pdf', '.csv'), header=False)

        if e0:
            param_barplot('Muscle Strain', e0, target[3])
            fpath = os.path.join(self.subject.results_exp_path, 
                    os.path.basename(target[3]))
            param_barplot('Muscle Strain', e0, fpath)
            df_e0 = pd.DataFrame.from_dict(e0, orient='index')
            df_e0.to_csv(fpath.replace('.pdf', '.csv'), header=False)


