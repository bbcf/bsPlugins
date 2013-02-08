from bsPlugins.DESeq import DESeqPlugin

path = 'testing_files/DESeq/'
DESeqPlugin()(**{'signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'], 'features':path+'features.bed', 'feature_type':3})

#DESeqPlugin()(**{'table':'tests/DESeq/table.tab'})#...
