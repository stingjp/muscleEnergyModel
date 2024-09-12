function [subj_model] = scaleModelMaxIsometricForces(subjectmass, subjectheight)
    
    %{
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
    %}

    import org.opensim.modeling.*

    % create struct for if we want height or not based on subject


    workdir = pwd;
    repodir = 'G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\';
%     cd(repodir)
%     % get the subject name and mass
%     load subjectmass.mat;
    %     # load height
    
    
%     cd(workdir);
    [~,trialname,~] = fileparts(pwd);
    cd ../
    [~,condname,~] = fileparts(pwd);
    cd ../
    [~,subjectname,~] = fileparts(pwd);
    experimentname = subjectname(1:4);
    cd(workdir);
    
    % # set the specific tension N/cm^2
    % specific_tension = 60; % # kPa?
    subj_mass = subjectmass.(genvarname(subjectname)); % kg
    
    % dont have height for all
    try
        subj_height = subjectheight.(genvarname(subjectname)); % m?
        haveheight = true;
    catch
        haveheight = false;
    end
    
    % get the generic model file.
    generic_model = Model('G:\Shared drives\Exotendon\muscleModel\muscleEnergyModel\Rajagopal_40_abdChg_passiveCalib_hippathsadjust_welkrecollect.osim');
    % subj_model = Model('./subject_updated.osim');
    subj_model = Model('./subject_redoarms.osim');


    function [output] = total_muscle_volume_regression_mass(mass)
        % """cm^3"""
        output = 91.0*mass + 588.0;
    end

    function [output] = total_muscle_volume_regression_massandheight(mass, height)
        % """cm^3"""
         output = 47.05*mass*height + 1289.6;
     end

    % # get the set of muscles
    generic_mset = generic_model.getMuscles();
    subj_mset = subj_model.getMuscles();
    
    % # scale using height and mass or just mass
    if haveheight
        % # TMV: total muscle volume (cm^3).
        generic_TMV_massandheight = total_muscle_volume_regression_massandheight(75.337, 1.7);
        subj_TMV_massandheight = total_muscle_volume_regression_massandheight(subj_mass, subj_height);

        % # go through and scale each muscle
        for im = 0:subj_mset.getSize()-1
            subj_muscle_name = subj_mset.get(im).getName();
            subj_muscle = subj_mset.get(subj_muscle_name);
            generic_muscle = generic_mset.get(subj_muscle_name);

            % # OFL: optimal fiber length (cm).
            generic_OFL = generic_muscle.get_optimal_fiber_length() * 100;
            subj_OFL = subj_muscle.get_optimal_fiber_length() * 100;

            scale_factor = (subj_TMV_massandheight/generic_TMV_massandheight) * (generic_OFL/subj_OFL);
%             print("Scaling '%s' muscle force by %f." % (subj_muscle_name, scale_factor))

            generic_force = generic_muscle.get_max_isometric_force();
            scaled_force = generic_force*scale_factor;
            subj_muscle.set_max_isometric_force(scaled_force);
        end
    else
        % # TMV: total muscle volume (cm^3).
        generic_TMV_mass = total_muscle_volume_regression_mass(75.337);
        subj_TMV_mass = total_muscle_volume_regression_mass(subj_mass);

        % # go through and scale each muscle
        for im = 0:subj_mset.getSize()-1
            subj_muscle_name = subj_mset.get(im).getName();
%             subj_muscle = subj_model.get(subj_muscle_name);
%             generic_muscle = generic_mset.get(subj_muscle_name);
            subj_muscle = subj_mset.get(subj_muscle_name);
            generic_muscle = generic_mset.get(subj_muscle_name);

            % # OFL: optimal fiber length (cm).
            generic_OFL = generic_muscle.get_optimal_fiber_length() * 100;
            subj_OFL = subj_muscle.get_optimal_fiber_length() * 100;

            scale_factor = (subj_TMV_mass/generic_TMV_mass) * (generic_OFL/subj_OFL);
            % print("Scaling '%s' muscle force by %f." % (subj_muscle_name, scale_factor))

            generic_force = generic_muscle.get_max_isometric_force();
            scaled_force = generic_force*scale_factor;
            subj_muscle.set_max_isometric_force(scaled_force);
        end
    end
    
    subj_model.print('subject_updated_Fscaled_redoarms.osim');
end


