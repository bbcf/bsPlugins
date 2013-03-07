from unittest2 import TestCase, skip
from bbcflib.btrack import track
from bsPlugins.TopGo import TopGoPlugin 
import os

path = 'testing_files/'

class Test_TopGoPlugin(TestCase):
    def setUp(self):
        self.plugin = TopGoPlugin()

    def test_TopGo_Single(self):
        self.plugin(gene_list=path+'ensemblids',assembly="NCBIM37")
	self.assertEqual([x[1] for x in self.plugin.output_files], ['pdf','txt'])

    def test_TopGo_Multiple(self):
        self.plugin(gene_list=path+'list_3clusters.txt', assembly="GRCh37.p5")
	self.assertEqual([x[1] for x in self.plugin.output_files], ['file','file'])

