from unittest2 import TestCase, skip
from bsPlugins.QuantifyTable import QuantifyTablePlugin
import os

path = 'testing_files/'

class Test_QuantifyTablePlugin(TestCase):
    def setUp(self):
        self.plugin = QuantifyTablePlugin()

    def test_quantify_table(self):
        self.plugin(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                               'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
        from bbcflib.btrack import track
        with track(self.plugin.output_files[0][0]) as t:
            content = list(t.read())
            for x in content: print x
            print "fields", t.fields
            self.assertListEqual(t.fields, ["chr","start","end","name","score0","score1"])
            self.assertEqual(len(content),9)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

