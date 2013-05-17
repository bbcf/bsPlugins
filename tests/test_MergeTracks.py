from unittest2 import TestCase, skip
from bsPlugins.MergeTracks import MergeTracksPlugin
import os

path = 'testing_files/'

class Test_MergeTracksPlugin(TestCase):
    def setUp(self):
        self.plugin = MergeTracksPlugin()

    def test_with_signals(self):
        self.plugin(**{'forward':path+'KO50.bedGraph', 'reverse':path+'WT50.bedGraph',
                         'assembly':'mm9', 'shift':0, 'format':'bedGraph'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),95)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_MergeTracks.py
