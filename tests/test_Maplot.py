from unittest2 import TestCase, skip
from bsPlugins.Maplot import MaplotPlugin
import os

path = 'testing_files/'

class Test_MaplotPlugin(TestCase):
    def setUp(self):
        self.plugin = MaplotPlugin()

    def test_with_signals(self):
        self.plugin(**{'input_type':'Signal','signals1':[path+'KO50.bedGraph'], 'signals2':[path+'WT50.bedGraph'],
                               'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})

    def test_with_table(self):
        self.plugin(**{'input_type':'Table','table':path+'genes_table.tab', 'assembly':'mm9'})

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

