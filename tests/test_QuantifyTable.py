from unittest2 import TestCase, skip
from bsPlugins.QuantifyTable import QuantifyTablePlugin
import os, shutil

path = 'testing_files/DESeq/'

class Test_DESeqPlugin(TestCase):
    def setUp(self):
        self.out_signals = ''

    def test_quantify_table(self):
        self.out_signals = QuantifyTablePlugin()(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                               'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
        from bbcflib.btrack import track
        with track(self.out_signals) as t:
            for line in t.read():
                print line
            raise
            #self.assertEqual(len(content),9)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)
                #shutil.rmtree(f)

