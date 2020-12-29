def scaleModelMaxIsometricForces(model_file, haveheight):

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

    import opensim

    generic_model = opensim.Model('./Rajagopal_40_abdChg_passiveCalib_hippathsadjust_dembiaoldmarkerset.osim')
    subj_model = opensim.Model(model_file)


    def total_muscle_volume_regression_mass(mass):
        """cm^3"""
        return 91.0*mass + 588.0

    def total_muscle_volume_regression_massandheight(mass, height):
        """cm^3"""
        return 47.05*mass*height + 1289.6

    # get the set of muscles
    generic_mset = generic_model.getMuscles()
    subj_mset = subj_model.getMuscles()

    # set the specific tension N/cm^2
    specific_tension = 60 # kPa?

    # scale using height and mass or just mass
    if haveheight:
        # TMV: total muscle volume (cm^3).
        generic_TMV_massandheight = total_muscle_volume_regression_massandheight(75.337, 1.7)
        subj_TMV_massandheight = total_muscle_volume_regression_massandheight(...,...)

        # go through and scale each muscle
        for im in range(subj_model.getMuscles().getSize()):
            subj_muscle_name = subj_mset.get(im).getName()
            subj_muscle = subj_model.get(subj_muscle_name)
            generic_muscle = generic_mset.get(subj_muscle_name)

            # OFL: optimal fiber length (cm).
            generic_OFL = generic_muscle.get_optimal_fiber_length() * 100
            subj_OFL = subj_muscle.get_optimal_fiber_length() * 100

            scale_factor = (subj_TMV_massandheight/generic_TMV_massandheight) * (generic_OFL/subj_OFL)
            print("Scaling '%s' muscle force by %f." % (subj_muscle_name, scale_factor))

            generic_force = generic_muscle.get_max_isometric_force()
            scaled_force = generic_force*scale_factor
            subj_muscle.set_max_isometric_force(scaled_force)

    else:
        # TMV: total muscle volume (cm^3).
        generic_TMV_mass = total_muscle_volume_regression_mass(75.337)
        subj_TMV_mass = total_muscle_volume_regression_mass(...)

        # go through and scale each muscle
        for im in range(subj_model.getMuscles().getSize()):
            subj_muscle_name = subj_mset.get(im).getName()
            subj_muscle = subj_model.get(subj_muscle_name)
            generic_muscle = generic_mset.get(subj_muscle_name)

            # OFL: optimal fiber length (cm).
            generic_OFL = generic_muscle.get_optimal_fiber_length() * 100
            subj_OFL = subj_muscle.get_optimal_fiber_length() * 100

            scale_factor = (subj_TMV_massandheight/generic_TMV_massandheight) * (generic_OFL/subj_OFL)
            print("Scaling '%s' muscle force by %f." % (subj_muscle_name, scale_factor))

            generic_force = generic_muscle.get_max_isometric_force()
            scaled_force = generic_force*scale_factor
            subj_muscle.set_max_isometric_force(scaled_force)


    # need to save the model then. 
    # subj_model.printToXML(target[0])




    # this remains from chris' original loaded walking version of the script?
        # # MIF: max isometric force (N).
        # MIF = muscle.get_max_isometric_force()
        # vfst = MIF * OFL / generic_TMV
        # # unitless.
        # vf = vfst / specific_tension
        # print('%s: %s, %s' % (muscle.getName(), vfst, vf))

