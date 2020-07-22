# post matlab solve
import opensim as osim

analyze = osim.AnalyzeTool()
analyze.setName("analyze")
analyze.setModelFilename("hanging_muscle.osim")
analyze.setStatesFileName("exampleHangingMuscle_states.sto")
analyze.updAnalysisSet().cloneAndAppend(osim.MuscleAnalysis())
analyze.updAnalysisSet().cloneAndAppend(osim.ProbeReporter())
analyze.updControllerSet().cloneAndAppend(
    osim.PrescribedController("exampleHangingMuscle_controls.sto"))
analyze.printToXML("exampleHangingMuscle_AnalyzeTool_setup.xml")