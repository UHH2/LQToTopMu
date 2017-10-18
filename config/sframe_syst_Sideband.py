#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     'NOMINAL'               : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'SCALE_upup'            : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_upup/NewExtrapolationFunction_TTbarDY.root'},
     'SCALE_upnone'          : {'ScaleVariationMuR':'up','ScaleVariationMuF':'none','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_upnone/NewExtrapolationFunction_TTbarDY.root'},
     'SCALE_noneup'          : {'ScaleVariationMuR':'none','ScaleVariationMuF':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_noneup/NewExtrapolationFunction_TTbarDY.root'},
     'SCALE_nonedown'        : {'ScaleVariationMuR':'none','ScaleVariationMuF':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_nonedown/NewExtrapolationFunction_TTbarDY.root'},
     'SCALE_downnone'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'none','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_downnone/NewExtrapolationFunction_TTbarDY.root'},
     'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/SCALE_downdown/NewExtrapolationFunction_TTbarDY.root'},
     'PU_up'                 : {'Systematic_PU':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/PU_up/NewExtrapolationFunction_TTbarDY.root'},
     'PU_down'               : {'Systematic_PU':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/PU_down/NewExtrapolationFunction_TTbarDY.root'},
     'MUID_up'               : {'Systematic_MuonID':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUID_up/NewExtrapolationFunction_TTbarDY.root'},
     'MUID_down'             : {'Systematic_MuonID':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUID_down/NewExtrapolationFunction_TTbarDY.root'},
     'MUTR_up'               : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUTR_up/NewExtrapolationFunction_TTbarDY.root'}, 
     'MUTR_down'             : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUTR_down/NewExtrapolationFunction_TTbarDY.root'}, 
     'BTAG_bc_up'            : {'Systematic_BTag':'up_bc','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/BTAG_bc_up/NewExtrapolationFunction_TTbarDY.root'},
     'BTAG_bc_down'          : {'Systematic_BTag':'down_bc','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/BTAG_bc_down/NewExtrapolationFunction_TTbarDY.root'},
     'BTAG_udsg_up'          : {'Systematic_BTag':'up_udsg','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/BTAG_udsg_up/NewExtrapolationFunction_TTbarDY.root'},
     'BTAG_udsg_down'        : {'Systematic_BTag':'down_udsg','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/BTAG_udsg_down/NewExtrapolationFunction_TTbarDY.root'},
     'MUISO_up'              : {'Systematic_MuonIso':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUISO_up/NewExtrapolationFunction_TTbarDY.root'},
     'MUISO_down'            : {'Systematic_MuonIso':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUISO_down/NewExtrapolationFunction_TTbarDY.root'},
     'MUTRK_up'              : {'Systematic_MuonTrk':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUTRK_up/NewExtrapolationFunction_TTbarDY.root'},
     'MUTRK_down'            : {'Systematic_MuonTrk':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUTRK_down/NewExtrapolationFunction_TTbarDY.root'},
     'ELEID_up'              : {'Systematic_EleID':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEID_up/NewExtrapolationFunction_TTbarDY.root'},
     'ELEID_down'            : {'Systematic_EleID':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEID_down/NewExtrapolationFunction_TTbarDY.root'},
     'ELEREC_up'             : {'Systematic_EleReco':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEREC_up/NewExtrapolationFunction_TTbarDY.root'},
     'ELEREC_down'           : {'Systematic_EleReco':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEREC_down/NewExtrapolationFunction_TTbarDY.root'},
     'ELETR_up'              : {'Systematic_EleTrigger':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELETR_up/NewExtrapolationFunction_TTbarDY.root'},
     'ELETR_down'            : {'Systematic_EleTrigger':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELETR_down/NewExtrapolationFunction_TTbarDY.root'},     
     'ELEFAKE_up'            : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEFAKE_up/NewExtrapolationFunction_TTbarDY.root'}, 
     'ELEFAKE_down'          : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ELEFAKE_down/NewExtrapolationFunction_TTbarDY.root'},      
     'MUFAKE_up'             : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUFAKE_up/NewExtrapolationFunction_TTbarDY.root'}, 
     'MUFAKE_down'           : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/MUFAKE_down/NewExtrapolationFunction_TTbarDY.root'}, 
     'PDF_up'                : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/PDF_up/NewExtrapolationFunction_TTbarDY.root'},
     'PDF_down'              : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/PDF_down/NewExtrapolationFunction_TTbarDY.root'},
     'ALPHA_up'              : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ALPHA_up/NewExtrapolationFunction_TTbarDY.root'},
     'ALPHA_down'            : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/ALPHA_down/NewExtrapolationFunction_TTbarDY.root'},
     #'NORMALIZATION_up'      : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/NORMALIZATION_up/NewExtrapolationFunction_TTbarDY.root'},
     #'NORMALIZATION_down'    : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/NORMALIZATION_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_ttbar_up'         : {'Systematic_TTbar':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_ttbar_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_ttbar_down'       : {'Systematic_TTbar':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_ttbar_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_dy_up'            : {'Systematic_DY':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_dy_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_dy_down'          : {'Systematic_DY':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_dy_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_db_up'            : {'Systematic_DB':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_db_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_db_down'          : {'Systematic_DB':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_db_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_st_up'            : {'Systematic_ST':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_st_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_st_down'          : {'Systematic_ST':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_st_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_wj_up'            : {'Systematic_WJ':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_wj_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_wj_down'          : {'Systematic_WJ':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_wj_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_qcd_up'           : {'Systematic_QCD':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_qcd_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_qcd_down'         : {'Systematic_QCD':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_qcd_down/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_ttv_up'           : {'Systematic_TTV':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_ttv_up/NewExtrapolationFunction_TTbarDY.root'},
     'RATE_ttv_down'         : {'Systematic_TTV':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/RATE_ttv_down/NewExtrapolationFunction_TTbarDY.root'},


     #'JEC_up'                     : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/JEC_up/NewExtrapolationFunction_TTbarDY.root'},
     #'JEC_down'                   : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/JEC_down/NewExtrapolationFunction_TTbarDY.root'},
     #'JER_up'                     : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/JER_up/NewExtrapolationFunction_TTbarDY.root'},
     #'JER_down'                   : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/35867fb_JetEta24/JER_down/NewExtrapolationFunction_TTbarDY.root'},
}
start_all_parallel = True

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude!'
    exit(-1)


def set_uncert_func(uncert_name):
    uncert = sys_uncerts[uncert_name]
    def do_set_uncert(element_tree):
        cycle = element_tree.getroot().find('Cycle')
        user_config = cycle.find('UserConfig')
        output_dir = cycle.get('OutputDirectory')

        cycle.set('OutputDirectory', os.path.join(output_dir, uncert_name+'/'))

        for name, value in uncert.iteritems():
            uc_item = list(i for i in user_config if i.get('Name') == name)
            assert uc_item, 'could not find item with name: %s' % name
            uc_item[0].set('Value', value)

    return do_set_uncert


from varial.extensions.sframe import SFrame
from varial import tools
if start_all_parallel:
    ToolChain = tools.ToolChainParallel
else:
    ToolChain = tools.ToolChain


class MySFrameBatch(SFrame):

    def configure(self):
        self.xml_doctype = self.xml_doctype + """
<!--
   <ConfigParse NEventsBreak="100000" FileSplit="0" AutoResubmit="0" />
   <ConfigSGE RAM ="2" DISK ="2" Mail="heiner@cern.de" Notification="as" Workdir="workdir"/>
-->
"""
        if os.path.exists(self.cwd + 'workdir'):
            opt = ' -rl --exitOnQuestion'
        else:
            opt = ' -sl --exitOnQuestion'

        self.exe = 'sframe_batch.py' + opt


sframe_tools = ToolChain(
    'SFrameUncertsCR',
    list(
        SFrame(
            cfg_filename=sys.argv[1],
            xml_tree_callback=set_uncert_func(uncert),
            name='SFrame_' + uncert,
            halt_on_exception=False,
        ) 
        for uncert in sys_uncerts
    )
)


if __name__ == '__main__':
    varial.tools.Runner(sframe_tools)
