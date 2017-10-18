import os
import subprocess

process = []
process.append(subprocess.Popen('sframe_main LQToTopMuAnalysis_JEC_down.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuAnalysis_JEC_up.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuAnalysis_JER_down.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuAnalysis_JER_up.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuSidebandAnalysis_EE_JEC_down.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuSidebandAnalysis_EE_JEC_up.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuSidebandAnalysis_EE_JER_down.xml', shell=True))
process.append(subprocess.Popen('sframe_main LQToTopMuSidebandAnalysis_EE_JER_up.xml', shell=True))

for thisproc in process:
    thisproc.wait()

srproc = subprocess.Popen('python sframe_syst.py LQToTopMuAnalysis.xml', shell=True)
srproc.wait()
crproc = subprocess.Popen('python sframe_syst_Sideband.py LQToTopMuSidebandAnalysis.xml', shell=True)
crproc.wait()
