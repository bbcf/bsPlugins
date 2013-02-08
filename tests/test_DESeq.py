from bsPlugins.DESeq import DESeqPlugin

path = 'testing_files/DESeq/'

DESeqPlugin()(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'], 'features':path+'features.bed', 'feature_type':3})

#DESeqPlugin()(**{'input_type':'Table','table':path+'genes_table.tab'})

