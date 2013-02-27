from unittest2 import TestCase, skip
from bsPlugins.Normalize import NormalizePlugin
import os

path = 'testing_files/'

class Test_NormalizePlugin(TestCase):

    def setUp(self):
        self.plugin = NormalizePlugin()
    def test_with_table(self):
        self.plugin(**{'table':path+'genes_table.tab', 'method':'total'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),20)
    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

