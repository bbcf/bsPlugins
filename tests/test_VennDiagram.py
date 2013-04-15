from unittest2 import TestCase, skip
from bsPlugins.VennDiagram import VennDiagramPlugin
import os

path = 'testing_files/'

class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramPlugin()

    def test_venn_diagram(self):
        self.plugin(**{'files':[path+'test1.bedGraph',path+'test2.bedGraph'],
                       'assembly':'mm9', 'format':'png'})
        #raise

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

