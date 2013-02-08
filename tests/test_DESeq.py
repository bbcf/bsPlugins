from bsPlugins.DESeq import DESeqPlugin
import os
import unittest2

path = 'testing_files/DESeq/'

class Test_DESeqPlugin(unittest2.TestCase):
    def setUp(self):
        self.out_signals = None
        self.out_table = None

    def test_with_signals(self):
        self.out_signals = DESeqPlugin()(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                               'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
        with open(self.out_signals,'rb') as f:
            print f.read()

    def test_with_table(self):
        self.out_table = DESeqPlugin()(**{'input_type':'Table','table':path+'genes_table.tab', 'assembly':'mm9'})
        with open(self.out_table,'rb') as f:
            print f.read()

    def tearDown(self):
        os.remove(self.out_signals)
        os.remove(self.out_table)


