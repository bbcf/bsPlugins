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
        self.fields = ['chr','start','end','score']
        self.kw = {'signals':[path+'test1.bedgraph', path+'test2.bedGraph'],
                   'assembly':'mm9', 'format':'bedGraph'}

    def test_intersect(self):
        self.intersect(**self.kw)
        with track(self.intersect.output_files[0][0]) as t:
            s = t.read(fields=self.fields)
            content = list(s)
            expected = [('chr1',10,15,17.0),('chr1',24,35,107.0)]
            self.assertListEqual(content,expected)

    def test_union(self):
        self.union(**self.kw)
        with track(self.union.output_files[0][0]) as t:
            s = t.read(fields=self.fields)
            content = list(s)
            expected = [('chr1',8,10,12.0),('chr1',10,15,17.0),('chr1',15,19,12.0),
                        ('chr1',21,24,17.0),('chr1',24,35,107.0),('chr1',35,39,90.0)]
            self.assertListEqual(content,expected)

    def test_complement(self):
        self.complement(**self.kw)
        with track(self.complement.output_files[0][0]) as t:
            s = t.read(fields=self.fields)
            content = list(s)
            expected = [('chr1',19,21,0.0)]
            self.assertListEqual(content,expected)

    def test_subtract(self):
        self.subtract(**self.kw)
        with track(self.subtract.output_files[0][0]) as t:
            s = t.read(fields=self.fields)
            content = list(s)
            expected = [('chr1',21,24,17.0)]
            self.assertListEqual(content,expected)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

