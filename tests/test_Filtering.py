from unittest2 import TestCase, skip
from bsPlugins.Filtering import FilteringPlugin
from bbcflib.track import track
import os
path = './testing_files/'
class Test_FilteringPlugin(TestCase):
    def setUp(self):
        self.plugin = FilteringPlugin()
    def test_with_signals(self):
#        self.plugin(**{'track':[path+'test3.bedGraph'],
#        'assembly':'mm10',
#        'minimum_score': 0 ,'maximum_score': 3 ,
#        'minimum_length':20,'maximum_length':30,
#        'chromosome_name':["chr1","chr2", "chrX"]})
#        self.plugin(**{'track':[path+'test3.bedGraph'],
#        'assembly':'mm10',
#        'minimum_score': 0 ,'maximum_score':3,
#        'minimum_length':20,'maximum_length':'',
#        'chromosome_name':["chr1","chr2", "chr3"]})
        self.plugin(**{'track':[path+'test3.bedGraph'],
        'assembly':'mm10',
        'minimum_score':"0",'maximum_score':"10",
        'minimum_length':'9','maximum_length':'500000',
        'chromosome_name':["chr1","chrX"]})
#    def tearDown(self):
#        for f in os.listdir('.'):
#            if f.startswith('tmp'):
#                os.system("rm -rf %s" % f)
