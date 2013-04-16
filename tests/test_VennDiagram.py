from unittest2 import TestCase, skip
from bsPlugins.VennDiagram import VennDiagramPlugin
import os

path = 'testing_files/'

class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramPlugin()

    def test_venn_diagram(self):
        self.plugin(**{'files':[path+'test1.bedGraph',path+'test2.bedGraph'],
                       'assembly':'mm9', 'format':'png'})
        #raise

    def tearDown(self):
        if 1:
            for f in os.listdir('.'):
                if f.startswith('tmp'):
                    os.system("rm -rf %s" % f)


"""
# nosetests --logging-filter=-tw2 test_VennDiagram.py

A
('chr1', 10, 15, 5.0)
('chr1', 21, 35, 17.0)
B
('chr1', 8, 19, 12.0)
('chr1', 24, 39, 90.0)

   10-----15   21---------------35      -> 29 tot coverage, 3/29 A, 10/29 B, 16/29 A|B
  8-----------19  24----------------39                      0.103   0.345    0.552
                                                      total 65% A,  90% B,   55% A|B   of tot coverage
"""

