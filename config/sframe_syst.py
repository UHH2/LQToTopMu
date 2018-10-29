#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     'NOMINAL'               : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'ALPHA_up'              : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'ALPHA_down'            : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'PDF_up'                : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'PDF_down'              : {'Systematic_TTbar':'nominal'}, #{...} is a dummy
     'SCALE_upup'            : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up'},
     'SCALE_upnone'          : {'ScaleVariationMuR':'up','ScaleVariationMuF':'none'},
     'SCALE_noneup'          : {'ScaleVariationMuR':'none','ScaleVariationMuF':'up'},
     'SCALE_nonedown'        : {'ScaleVariationMuR':'none','ScaleVariationMuF':'down'},
     'SCALE_downnone'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'none'},
     'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down'},
     'PU_up'                 : {'Systematic_PU':'up'},
     'PU_down'               : {'Systematic_PU':'down'},
     'MUID_up'               : {'Systematic_MuonID':'up'},
     'MUID_down'             : {'Systematic_MuonID':'down'},
     'MUTR_up'               : {'Systematic_MuonTrigger':'up'},
     'MUTR_down'             : {'Systematic_MuonTrigger':'down'},
     'BTAG_bc_up'            : {'Systematic_BTag':'up_bc'},
     'BTAG_bc_down'          : {'Systematic_BTag':'down_bc'},
     'BTAG_udsg_up'          : {'Systematic_BTag':'up_udsg'},
     'BTAG_udsg_down'        : {'Systematic_BTag':'down_udsg'},
     'MUISO_up'              : {'Systematic_MuonIso':'up'},
     'MUISO_down'            : {'Systematic_MuonIso':'down'},
     'MUTRK_up'              : {'Systematic_MuonTrk':'up'},
     'MUTRK_down'            : {'Systematic_MuonTrk':'down'},
     'ELEID_up'              : {'Systematic_EleID':'up'},
     'ELEID_down'            : {'Systematic_EleID':'down'},
     'ELEREC_up'             : {'Systematic_EleReco':'up'},
     'ELEREC_down'           : {'Systematic_EleReco':'down'},     
     'ELETR_up'              : {'Systematic_TTbar':'nominal'}, #dummy
     'ELETR_down'            : {'Systematic_TTbar':'nominal'}, #dummy
     'ELEFAKE_up'            : {'Systematic_EleFake':'up'},
     'ELEFAKE_down'          : {'Systematic_EleFake':'down'},
     'MUFAKE_up'             : {'Systematic_MuFake':'up'},
     'MUFAKE_down'           : {'Systematic_MuFake':'down'},
     'RATE_ttbar_up'         : {'Systematic_TTbar':'up'},
     'RATE_ttbar_down'       : {'Systematic_TTbar':'down'},
     'RATE_dy_up'            : {'Systematic_DY':'up'},
     'RATE_dy_down'          : {'Systematic_DY':'down'},
     'RATE_st_up'            : {'Systematic_ST':'up'},
     'RATE_st_down'          : {'Systematic_ST':'down'},
     'RATE_wj_up'            : {'Systematic_WJ':'up'},
     'RATE_wj_down'          : {'Systematic_WJ':'down'},
     'RATE_db_up'            : {'Systematic_DB':'up'},
     'RATE_db_down'          : {'Systematic_DB':'down'},
     'RATE_qcd_up'           : {'Systematic_QCD':'up'},
     'RATE_qcd_down'         : {'Systematic_QCD':'down'},
     'RATE_ttv_up'           : {'Systematic_TTV':'up'},
     'RATE_ttv_down'         : {'Systematic_TTV':'down'},
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
    'SFrameUncertsSR',
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
