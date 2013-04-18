from unittest2 import TestCase, skip
from bsPlugins.VennDiagram import VennDiagramPlugin
import os

path = 'testing_files/'

class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramPlugin()

    #@skip('')
    def test_venn_diagram(self):
        self.plugin(**{'files':[path+'test1.bedGraph',path+'test2.bedGraph'],
                       'assembly':'mm9', 'format':'png'})
        #raise

    def test_allkinds(self):
        from bbcflib.bFlatMajor.figure import venn
        D1 = {'Group1':126}
        D2 = {'Group1':126, 'Group2':247, 'Group1|Group2':50}
        D31 = {'A':521, 'B':14, 'C':290, 'A|B':11, 'A|C':100, 'B|C':4, 'A|B|C':1}
        D32 = {'A':521, 'B':300, 'C':290, 'A|B':11, 'A|C':100, 'B|C':44, 'A|B|C':5}
        D4 = {'A':230, 'B':230, 'C':230, 'D':230,
              'A|B':80, 'A|C':80, 'A|D':80, 'B|C':80, 'B|D':80, 'C|D':80,
              'A|B|C':30, 'A|B|D':30, 'A|C|D':30, 'B|C|D':30, 'A|B|C|D':10}
        D5 = {'A':321, 'B':540, 'C':490, 'D':1000, 'E':500,
              'A|B':200, 'A|C':200, 'A|D':200, 'A|E':200,
              'B|C':200, 'B|D':200, 'B|E':200, 'C|D':200, 'C|E':200, 'D|E':200,
              'A|B|C':130, 'A|B|D':130, 'A|B|E':130, 'A|C|D':130, 'A|C|E':130, 'A|D|E':130,
              'B|C|D':130, 'B|C|E':130, 'B|D|E':130, 'C|D|E':130,
              'A|B|C|D':60, 'A|B|C|E':60, 'A|B|D|E':60, 'A|C|D|E':60, 'B|C|D|E':60, 'A|B|C|D|E':20}
        venn(D1,output='tmp/d1.pdf',legend=['/adsf/qwer/zxcv/file1.bed'])
        venn(D2,output='tmp/d2.pdf',legend=['file1,file2'])
        venn(D31,output='tmp/d3.1.pdf',legend=['file1','file2','file3'])
        venn(D32,output='tmp/d3.2.pdf',legend=['file1','file2','file3'])
        venn(D4,output='tmp/d4.pdf',legend=['file1','file2','file3','file4'])
        #venn(D5,output='tmp/d5.pdf')

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

