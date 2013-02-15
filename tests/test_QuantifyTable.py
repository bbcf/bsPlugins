from unittest2 import TestCase, skip
from bbcflib.btrack import track
from bsPlugins.QuantifyTable import QuantifyTablePlugin
import os

path = 'testing_files/'

class Test_QuantifyTablePlugin(TestCase):
    def setUp(self):
        self.plugin = QuantifyTablePlugin()

    def test_quantify_table_sql(self):
        self.plugin(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                       'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9', 'format':'sql'})
        with track(self.plugin.output_files[0][0]) as t:
            s = t.read()
            content = list(s)
            self.assertListEqual(s.fields, ["chr","start","end","name","score0","score1"])
            self.assertEqual(len(content),9)

    def test_quantify_table_text(self):
        self.plugin(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'],
                       'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9', 'format':'txt'})
        with track(self.plugin.output_files[0][0], fields=["chr","start","end","name","score0","score1"]) as t:
            s = t.read()
            content = list(s)
            self.assertEqual(len(content),9)

    def tearDown(self):
        for f in os.listdir('.'):
            if f.startswith('tmp'):
                os.system("rm -rf %s" % f)

