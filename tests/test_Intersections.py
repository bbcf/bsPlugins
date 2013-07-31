from unittest2 import TestCase, skip
from bsPlugins.Intersections import IntersectionsPlugin
import os

path = 'testing_files/intersections/'

class Test_IntersectionsPlugin(TestCase):
    def setUp(self):
        self.plugin = IntersectionsPlugin()

    def test_intersections(self):
        self.plugin(**{'SigMulti':{'files':[path+'genes_differential_KE1-KE2.txt_UP',
                                            path+'genes_differential_KH1-KH2.txt_UP',
                                            path+'genes_differential_R1-R3.txt_UP']
                                  },
                       'column':1})

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

# nosetests --logging-filter=-tw2 test_Intersections.py
