from unittest2 import TestCase, skip
from bsPlugins.DESeq import DESeqPlugin
import os

path = 'testing_files/'

class Test_DESeqPlugin(TestCase):
    def setUp(self):
        self.plugin = DESeqPlugin()

    def test_with_signals(self):
        self.plugin(**{'input_type':'Signal',
                       'Group1':{'signals1':[path+'KO50.bedGraph']},
                       'Group2':{'signals2':[path+'WT50.bedGraph']},
                       'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),9)

    def test_with_table(self):
        self.plugin(**{'input_type':'Table','table':path+'genes_table.tab', 'assembly':'mm9'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),20)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_DESeq.py
