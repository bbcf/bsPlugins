from bsPlugins.DESeq import DESeqPlugin
import os
from unittest2 import TestCase, skip

path = 'testing_files/DESeq/'

class Test_DESeqPlugin(TestCase):
    def setUp(self):
        self.out_signals = ''
        self.out_table = ''

    def test_with_signals(self):
        self.out_signals = DESeqPlugin()(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                               'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
        with open(self.out_signals,'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),9)

    def test_with_table(self):
        self.out_table = DESeqPlugin()(**{'input_type':'Table','table':path+'genes_table.tab', 'assembly':'mm9'})
        with open(self.out_table,'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),20)

    def tearDown(self):
        if os.path.exists(self.out_signals): os.remove(self.out_signals)
        if os.path.exists(self.out_table): os.remove(self.out_table)


