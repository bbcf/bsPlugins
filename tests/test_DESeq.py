from bsPlugins.DESeq import DESeqPlugin
import os

path = 'testing_files/DESeq/'

out = DESeqPlugin()(**{'input_type':'Signal','signals':[path+'KO50.bedGraph', path+'WT50.bedGraph'], 'features':path+'features.bed', 'feature_type':3, 'assembly':'mm9'})
with open(out,'rb') as f:
    print f.read()
os.remove(out)


#DESeqPlugin()(**{'input_type':'Table','table':path+'genes_table.tab', 'assembly':'mm9'})

