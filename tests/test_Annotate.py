from unittest2 import TestCase, skip
from bsPlugins.Annotate import AnnotatePlugin
import os

path = 'testing_files/'

class Test_AnnotatePlugin(TestCase):
    def setUp(self):
        self.plugin = AnnotatePlugin()

    def test_with_signals(self):
        self.plugin(**{'track':path+'KO50.bedGraph', 'assembly':'mm9',
                       'promoter':0, 'intergenic':0, 'UTR':0})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),50)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_Annotate.py
