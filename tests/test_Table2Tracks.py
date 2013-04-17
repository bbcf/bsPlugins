from unittest2 import TestCase, skip
from bsPlugins.Table2Tracks import Table2TracksPlugin
import os

path = 'testing_files/'

class Test_Table2TracksPlugin(TestCase):
    def setUp(self):
        self.plugin = Table2TracksPlugin()

    def test_A(self):
        self.plugin(**{'tableFile':path+'Dup_vs_Ctrl_resSelectedFrags_fromSmoothed_KCTD13_part.txt', 'id_columns':'6,7,9,10', 'assembly':'hg19', 'format':'bedGraph'})
        self.assertEqual(len(self.plugin.output_files[0]),4)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

