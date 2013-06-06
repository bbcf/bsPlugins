from unittest2 import TestCase, skip
from bbcflib.btrack import track
from bsPlugins.Ratios import RatiosPlugin
import os

path = 'testing_files/'

class Test_RatiosPlugin(TestCase):
    def setUp(self):
        self.plugin = RatiosPlugin()

    def test_ratios(self):
        self.plugin(**{'numerator':path+'KO50.bedGraph', 'denominator':path+'WT50.bedGraph', 'format':'bedGraph'})
        with track(self.plugin.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
        #print content
        #raise

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_Ratios.py
