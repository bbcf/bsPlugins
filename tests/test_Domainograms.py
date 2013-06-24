from unittest2 import TestCase, skip
from bsPlugins.Domainograms import DomainogramsPlugin
import os, tarfile

path = 'testing_files/'

class Test_DomainogramsPlugin(TestCase):
    def setUp(self):
        self.plugin = DomainogramsPlugin()

    def test_domainograms(self):
        self.plugin(**{'sample':path+'segToFrag_part.bedGraph.gz', 'name':'test_domainogram', 'region':'chr16', 'wmax_domainograms':'50', 'wmax_BRICKS':'500'})
        res_tar = tarfile.open(self.plugin.output_files[0][0], "r:gz")
        self.assertEqual(len(res_tar.getmembers()),6, msg="at leat one expected result file is missing") # the archive contains 6 files
        self.assertTrue([x.size>0 for x in res_tar], msg="at least one of the output file is empty") # no file is empty

    def test_domainograms_upDown(self):
        self.plugin(**{'sample':path+'segToFrag_part.bedGraph.gz', 'name':'test_domainogram', 'region':'chr16:29927651-29928325', 'wmax_domainograms':'50', 'wmax_BRICKS':'500'})
        res_tar = tarfile.open(self.plugin.output_files[0][0], "r:gz")
        self.assertEqual(len(res_tar.getmembers()),12, msg="at leat one expected result file is missing") # the archive contains 12 files
        self.assertTrue([x.size>0 for x in res_tar], msg="at least one of the output file is empty")

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_Smoothing.py
