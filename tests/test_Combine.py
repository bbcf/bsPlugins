from unittest2 import TestCase, skip
from bbcflib.btrack import track
from bsPlugins.Combine import IntersectPlugin,UnionPlugin,ComplementPlugin,SubtractPlugin
import os

path = 'testing_files/'

class Test_CombinePlugin(TestCase):
    def setUp(self):
        self.intersect = IntersectPlugin()
        self.union = UnionPlugin()
        self.complement = ComplementPlugin()
        self.subtract = SubtractPlugin()
        self.kw = {'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                   'assembly':'mm9', 'format':'sql'}

    def test_intersect(self):
        self.intersect(**self.kw)
        with track(self.intersect.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            self.assertEqual(len(content),24)

    def test_union(self):
        self.union(**self.kw)
        with track(self.union.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            self.assertEqual(len(content),95)

    def test_complement(self):
        self.complement(**self.kw)
        with track(self.complement.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            self.assertEqual(len(content),28)

    def test_subtract(self):
        self.subtract(**self.kw)
        with track(self.subtract.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            self.assertEqual(len(content),39)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

