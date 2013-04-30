from unittest2 import TestCase, skip
from bsPlugins.List2Track import List2TrackPlugin
import os

path = 'testing_files/'

class Test_List2TrackPlugin(TestCase):
    def setUp(self):
        self.plugin = List2TrackPlugin()

    #@skip('')
    def test_with_genes_list(self):
        self.plugin(**{'ids_list':path+'list_yeast_genes.txt', 'feature_type':'genes',
                       'assembly':'sacCer2', 'format':'bed'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            first = f.readline().split()
            self.assertEqual(first,['chrII','45974','46367','YBL092W|RPL32','0','+'])

    #@skip('')
    def test_with_type_only(self):
        self.plugin(**{'feature_type':'genes','assembly':'sacCer2', 'format':'bed'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            first = f.readline().split()
            self.assertEqual(first,['chrI','136915','137512','YAL008W|FUN14','0','+'])

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

