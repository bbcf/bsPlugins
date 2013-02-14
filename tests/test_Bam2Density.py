from unittest2 import TestCase, skip
from bsPlugins.Bam2Density import Bam2DensityPlugin
import os, shutil

path = 'testing_files/'

class Test_DESeqPlugin(TestCase):
    def setUp(self):
        self.out = ''

    def test_bam2density_bedGraph_nomerge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'bedGraph', 'merge_strands':-1})
        with open(self.out[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),447)

    def test_bam2density_bedGraph_merge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'bedGraph', 'merge_strands':1})
        with open(self.out[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),322)

    def test_bam2density_sql_nomerge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'sql', 'merge_strands':-1})
        self.assertEqual(len(self.out),2)

    def test_bam2density_sql_merge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'sql', 'merge_strands':1})
        self.assertEqual(len(self.out),1)

    def test_bam2density_other_nomerge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'wig', 'merge_strands':-1})
        with open(self.out[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),613)

    def test_bam2density_other_merge(self):
        self.out = Bam2DensityPlugin()(**{'sample':path+'Gapdh.bam', 'format':'wig', 'merge_strands':1})
        with open(self.out[0][0],'rb') as f:
            content = f.readlines()
            self.assertEqual(len(content),437)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)
                #shutil.rmtree(f)

