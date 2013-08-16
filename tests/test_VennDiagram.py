from unittest2 import TestCase, skip
from bsPlugins.VennDiagram import VennDiagramPlugin
import os

path = 'testing_files/venn/'


class Test_VennDiagramPlugin(TestCase):
    def setUp(self):
        self.plugin = VennDiagramPlugin()

    #@skip('')
    def test_venn_diagram_2way(self):
        self.plugin(**{'SigMulti':{'files':[path+'../test1.bedGraph',path+'../test2.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'coverage %', 'names':''})
    @skip('')
    def test_venn_diagram_3way(self):
        self.plugin(**{'SigMulti':{'files':[path+'test1.bedGraph',
                                            path+'test2.bedGraph',path+'test3.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'coverage %',
                       'names':"Group1 Group2 Group3"})
    @skip('')
    def test_venn_diagram_4way(self):
        self.plugin(**{'SigMulti':{'files':[path+'test1.bedGraph',path+'test2.bedGraph',
                                            path+'test3.bedGraph',path+'test4.bedGraph']},
                       'assembly':'mm9', 'format':'png', 'type':'tag count', 'names':''})

    #@skip('Set par() options so that it looks like it should in pdf format')
    def test_allkinds(self):
        format = 'pdf'
        from bbcflib.gfminer.figure import venn
        D1 = {'A':126}
        D2 = {'Albert':126, 'Barthur':247, 'Albert|Barthur':50}
        D31 = {'Ar':521, 'Bi':14, 'Co':290, 'Ar|Bi':11, 'Ar|Co':100, 'Bi|Co':4, 'Ar|Bi|Co':1}
        D32 = {'A':521, 'B':300, 'C':290, 'A|B':11, 'A|C':100, 'B|C':44, 'A|B|C':5}
        D4 = {'A':210, 'B':220, 'C':230, 'D':240,
              'A|B':80, 'A|C':80, 'A|D':80, 'B|C':80, 'B|D':80, 'C|D':80,
              'A|B|C':30, 'A|B|D':30, 'A|C|D':30, 'B|C|D':30, 'A|B|C|D':10}
        venn(D1,output=path+'d1.'+format,legend=['file1.bed'],format=format)
        venn(D2,output=path+'d2.'+format,format=format)
        venn(D31,output=path+'d3.1.'+format,format=format)
        venn(D32,output=path+'d3.2.'+format,legend=['file1','file2','file3'],format=format)
        venn(D4,output=path+'d4.'+format,legend=['file1','file2','file3','file4'],format=format)

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

