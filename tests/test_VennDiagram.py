from unittest2 import TestCase, skip
from bsPlugins.VennDiagram import VennDiagramPlugin
import os

path = 'testing_files/'


class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramPlugin()

    #@skip('')
    def test_venn_diagram_2way(self):
        self.plugin(**{'SigMulti':{'files':[path+'test1.bedGraph',path+'test2.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'coverage %'})
    @skip('')
    def test_venn_diagram_3way(self):
        self.plugin(**{'SigMulti':{'files':[path+'venn/test1.bedGraph',
                                            path+'venn/test2.bedGraph',path+'venn/test3.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'coverage %'})
    @skip('')
    def test_venn_diagram_4way(self):
        self.plugin(**{'SigMulti':{'files':[path+'venn/test1.bedGraph',path+'venn/test2.bedGraph',
                                            path+'venn/test3.bedGraph',path+'venn/test4.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'tag count'})

    @skip('Set par() options so that it looks like it should in pdf format')
    def test_allkinds(self):
        format = 'png'
        from bbcflib.bFlatMajor.figure import venn
        D1 = {'A':126}
        D2 = {'A':126, 'B':247, 'A|B':50}
        D31 = {'A':521, 'B':14, 'C':290, 'A|B':11, 'A|C':100, 'B|C':4, 'A|B|C':1}
        D32 = {'A':521, 'B':300, 'C':290, 'A|B':11, 'A|C':100, 'B|C':44, 'A|B|C':5}
        D4 = {'A':230, 'B':230, 'C':230, 'D':230,
              'A|B':80, 'A|C':80, 'A|D':80, 'B|C':80, 'B|D':80, 'C|D':80,
              'A|B|C':30, 'A|B|D':30, 'A|C|D':30, 'B|C|D':30, 'A|B|C|D':10}
        venn(D1,output='temp/d1.'+format,legend=['file1.bed'],format=format)
        venn(D2,output='temp/d2.'+format,legend=['file1','file2'],format=format)
        venn(D31,output='temp/d3.1.'+format,legend=['file1','file2','file3'],format=format)
        venn(D32,output='temp/d3.2.'+format,legend=['file1','file2','file3'],format=format)
        venn(D4,output='temp/d4.'+format,legend=['file1','file2','file3','file4'],format=format)
        raise

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

