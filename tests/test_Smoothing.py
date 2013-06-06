from unittest2 import TestCase, skip
from bsPlugins.Smoothing import SmoothingPlugin
from bbcflib.btrack import track
import os

path = 'testing_files/'

class Test_SmoothingPlugin(TestCase):
    def setUp(self):
        self.plugin = SmoothingPlugin()

    def test_smoothing(self):
        self.plugin(**{'track':path+'KO50.bedGraph', 'assembly':'mm9', 'format':'bedGraph'})
        with track(self.plugin.output_files[0][0]) as t:
            content = list(t.read())
            self.assertEqual(len(content),501)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_Smoothing.py
