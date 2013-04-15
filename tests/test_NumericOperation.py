from unittest2 import TestCase, skip
from bsPlugins.NumericOperation import NumericOperationPlugin
from bbcflib.btrack import track
import os
path = './testing_files/'
class Test_NumericOperationPlugin(TestCase):
    def setUp(self):
        self.plugin = NumericOperationPlugin()
    def test_with_signals(self):
        self.plugin(**{'track':[path+'test3.bedGraph'],
        'assembly':'mm9', "format":"sql",'function':'sqrt'})
        self.plugin(**{'track':[path+'test1.bedGraph',path+'CTCF_deconv.sql', path+'test2.bedGraph", path+'WT50.bedGraph', path+'KO50.bedGraph'],
        'assembly':'mm9', "format":"",'function':'log2'})
        self.plugin(**{'track':[path+'bigWigExamplehg19.bw'],'assembly':'hg19', "format":"bedGraph",'function':'log2'})
    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)
