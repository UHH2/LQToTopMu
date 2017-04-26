#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     'SCALE_upup'            : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_upup/ExtrapolationFunction.root'},
     'SCALE_upnone'          : {'ScaleVariationMuR':'up','ScaleVariationMuF':'none','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_upnone/ExtrapolationFunction.root'},
     'SCALE_noneup'          : {'ScaleVariationMuR':'none','ScaleVariationMuF':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_noneup/ExtrapolationFunction.root'},
     'SCALE_nonedown'        : {'ScaleVariationMuR':'none','ScaleVariationMuF':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_nonedown/ExtrapolationFunction.root'},
     'SCALE_downnone'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'none','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_downnone/ExtrapolationFunction.root'},
     'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/SCALE_downdown/ExtrapolationFunction.root'},
     'PU_up'                 : {'Systematic_PU':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/PU_up/ExtrapolationFunction.root'},
     'PU_down'               : {'Systematic_PU':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/PU_down/ExtrapolationFunction.root'},
     'MUID_up'               : {'Systematic_MuonID':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUID_up/ExtrapolationFunction.root'},
     'MUID_down'             : {'Systematic_MuonID':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUID_down/ExtrapolationFunction.root'},
     'MUTR_up'               : {'Systematic_MuonTrigger':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUTR_up/ExtrapolationFunction.root'},
     'MUTR_down'             : {'Systematic_MuonTrigger':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUTR_down/ExtrapolationFunction.root'},
     'BTAG_bc_up'            : {'Systematic_BTag':'up_bc','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/BTAG_bc_up/ExtrapolationFunction.root'},
     'BTAG_bc_down'          : {'Systematic_BTag':'down_bc','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/BTAG_bc_down/ExtrapolationFunction.root'},
     'BTAG_udsg_up'          : {'Systematic_BTag':'up_udsg','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/BTAG_udsg_up/ExtrapolationFunction.root'},
     'BTAG_udsg_down'        : {'Systematic_BTag':'down_udsg','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/BTAG_udsg_down/ExtrapolationFunction.root'},
     'MUISO_up'              : {'Systematic_MuonIso':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUISO_up/ExtrapolationFunction.root'},
     'MUISO_down'            : {'Systematic_MuonIso':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/MUISO_down/ExtrapolationFunction.root'},
     'ELEID_up'              : {'Systematic_EleID':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ELEID_up/ExtrapolationFunction.root'},
     'ELEID_down'            : {'Systematic_EleID':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ELEID_down/ExtrapolationFunction.root'},
     'ELEREC_up'             : {'Systematic_EleReco':'up','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ELEREC_up/ExtrapolationFunction.root'},
     'ELEREC_down'           : {'Systematic_EleReco':'down','filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ELEREC_down/ExtrapolationFunction.root'},
     'PDF_up'                : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/PDF_up/ExtrapolationFunction.root'},
     'PDF_down'              : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/PDF_down/ExtrapolationFunction.root'},
     'ALPHA_up'              : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ALPHA_up/ExtrapolationFunction.root'},
     'ALPHA_down'            : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/ALPHA_down/ExtrapolationFunction.root'},
     'NORMALIZATION_up'      : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/NORMALIZATION_up/ExtrapolationFunction.root'},
     'NORMALIZATION_down'    : {'filepath_alpha':'/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TTbarSideband/EE_Channel/Optimization/36773fb/NORMALIZATION_down/ExtrapolationFunction.root'},
     #'JEC_up'               : {'jecsmear_direction':'up'},
     #'JEC_down'             : {'jecsmear_direction':'down'},
     #'JER_up'               : {'jersmear_direction':'up'},
     #'JER_down'             : {'jersmear_direction':'down'},
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
    'SFrameUncerts',
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
