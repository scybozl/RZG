# Set up for pp -> ttbar  process
cmds.append("cd /Herwig/MatrixElements")
cmds.append("insert SimpleQCD:MatrixElements[0] MEHeavyQuark")

cmds.append("set /Herwig/EventHandlers/LHCHandler:Weighted On")

cmds.append("insert /Herwig/Generators/LHCGenerator:AnalysisHandlers 0 /Herwig/Analysis/HepMCFile")
cmds.append("set /Herwig/Analysis/HepMCFile:PrintEvent 1000000")
cmds.append("set /Herwig/Analysis/HepMCFile:Format GenEvent")
cmds.append("set /Herwig/Analysis/HepMCFile:Units GeV_mm")

cmds.append("set /Herwig/EventHandlers/LHCHandler:CascadeHandler NULL")
cmds.append("set /Herwig/EventHandlers/LHCHandler:HadronizationHandler NULL")

