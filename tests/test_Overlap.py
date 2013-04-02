from unittest2 import TestCase, skip
from bbcflib.btrack import track
from bsPlugins.Overlap import OverlapPlugin
import os

path = 'testing_files/'

class Test_OverlapPlugin(TestCase):
    def setUp(self):
        self.plugin = OverlapPlugin()

    def test_overlap(self):
        self.plugin(**{'input_type':'Signal','filter':path+'peaks.bedGraph',
                       'features':path+'features.bed','feature_type':3,'assembly':'mm9','format':'sql'})
        with track(self.plugin.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            print content
            self.assertEqual(len(content),3)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

