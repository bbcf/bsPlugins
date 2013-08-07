from unittest2 import TestCase, skip
from bsPlugins.Statistics import FileStatisticsPlugin
import os

path = 'testing_files/'


class Test_StatisticsPlugin(TestCase):
    def setUp(self):
        self.plugin = FileStatisticsPlugin()

    #@skip('')
    def test_stats(self):
        self.plugin(**{'sample':path+'test1.bedGraph', 'assembly':'mm9'})
        with open(self.plugin.output_files[0][0],'rb') as f:
            content = f.readlines()
            for x in content: print x.strip()
            raise

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)


# nosetests --logging-filter=-tw2 test_Statistics.py

