from unittest2 import TestCase, skip
from bsPlugins.Pca import PcaPlugin
import os

path = 'testing_files/'

class Test_PcaPlugin(TestCase):
    def setUp(self):
        self.plugin = PcaPlugin()

    #def test_pca_select(self):
    #    self.plugin(**{'table': path+'genes_expression_100_Leanne.txt','columns': "1,2,3"})

    #def test_pca_all(self):
    #    self.plugin(**{'table': path+'genes_expression_100_Leanne.txt','columns': "all"})

    def test_pca_rpkm(self):
        self.plugin(**{'table': path+'genes_expression_100_Leanne.txt','columns': "rpkm"})

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --nocapture --logging-filter=-tw2 test_Pca.py

