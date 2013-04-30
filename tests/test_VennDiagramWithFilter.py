from unittest2 import TestCase, skip
from bsPlugins.VennDiagramWithFilter import VennDiagramWithFilterPlugin
import os

path = 'testing_files/'


class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramWithFilterPlugin()

    def test_VennDiagramPlugin(self):
        infile=path+"Dup_vs_Ctrl_resSelectedFrags_fromSmoothed_KCTD13_part.txt"
        s_cols="5,7,8"
        s_filters=">=3,<7,>15"
        self.plugin(**{'table':infile,'id_columns':s_cols,'filters':s_filters,'format':'png'})


    def tearDown(self):
        if 1:
            for f in os.listdir('.'):
                if f.startswith('tmp'):
                    os.system("rm -rf %s" % f)


"""
# nosetests --logging-filter=-tw2 test_VennDiagramWithFilter.py

"""

