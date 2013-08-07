from bsPlugins import *
from bbcflib.track import track,stats
from bbcflib import genrep
import os

output_list = ['txt','pdf']

class FileStatisticsForm(BaseForm):
    child = twd.HidingTableLayout()
    sample = twb.BsFileField(label='Input file: ',
        help_text='Select the file to examine',
        validator=twb.BsFileFieldValidator(required=True))
    output = twf.SingleSelectField(label='Output: ',
                                   options=output_list,
                                   prompt_text=None,
                                   help_text='Type of report')
    by_chrom = twf.CheckBox(label='By chromosome: ',
                            value=False,
                            help_text='Split statistics by chromosome (default: whole genome)')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'sample', 'type':'track', 'required':True},
        {'id':'output', 'type':'list', 'required':True},
        {'id':'by_chrom', 'type':'boolean', 'required':True}]
out_parameters = [{'id':'stats', 'type':'file'},
                  {'id':'pdf', 'type':'file'}]


class FileStatisticsPlugin(BasePlugin):
    """Calculates diverse statistics from a track file,
    such as a distribution of scores and feature lengths, and prints them to the output file."""
    info = {
        'title': 'Basic track statistics',
        'description': __doc__,
        'path': ['Analysis', 'File Statistics'],
        'output': FileStatisticsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def _plot_pdf(self,filename,stats,title=""):
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri as numpy2ri
        robjects.r('pdf("%s",paper="a4",height=8,width=8)' %filename)
        for chrom,st in stats.iteritems():
            if chrom:
                _title = title+":"+chrom
            else:
                _title = title
            if 'feat_stats' in st:
                fst = st['feat_stats']
                robjects.r.assign('len',numpy2ri.numpy2ri(fst[1].keys()))
                robjects.r.assign('num',numpy2ri.numpy2ri(fst[1].values()))
                robjects.r.assign('ylim',max(10,fst[0]))
                robjects.r.assign('med',fst[2][5])
                robjects.r.assign('men',fst[2][3])
                robjects.r.assign('sdv',fst[2][4])
                robjects.r("""
ypos=1
len=as.numeric(len)
num=as.numeric(num)
par(lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=.8,mar=c(5,5,1,1),las=1,pch=20)
plot(len,num,type='h',main='%s',xlab='Feature Length',ylab='Frequency',ylim=c(1,ylim),log='y')
abline(v=med,col='red')
text(med,ylim,paste("median",med,sep="="),col='red',pos=4)
abline(h=ylim[1],col='green')
mtext(paste(ylim[1],"features"),side=2,at=10,col='green',las=1)
arrows(men-sdv,ypos,men+sdv,ypos,angle=90,code=3,length=.15,col='blue')
points(men,ypos,pch=19,col='blue')
"""%_title)
            if 'score_stats' in st:
                sst = st['score_stats']
                robjects.r.assign('score',numpy2ri.numpy2ri(sst[0].keys()))
                robjects.r.assign('num',numpy2ri.numpy2ri(sst[0].values()))
                robjects.r.assign('med',sst[1][5])
                robjects.r.assign('men',sst[1][3])
                robjects.r.assign('sdv',sst[1][4])
                robjects.r("""
ypos=1
score=as.numeric(score)
num=as.numeric(num)
par(lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=0.8,mar=c(5,5,1,1),las=1,pch=20)
plot(score,num,type='h',main='%s',xlab='Score',ylab='Frequency',log='y')
abline(v=med,col='red')
text(med,ylim[1],paste("median",med,sep="="),col='red',pos=4)
arrows(men-sdv,ypos,men+sdv,ypos,angle=90,code=3,length=.15,col='blue')
points(men,ypos,pch=19,col='blue')
"""%_title)
        robjects.r("dev.off()")
        return None

    def __call__(self, **kw):
        self.debug(**kw)
        sample = track(kw['sample'],chrmeta="guess")
        by_chrom = kw.get('by_chrom',False)
        if isinstance(by_chrom, basestring):
            by_chrom = (by_chrom.lower() in ['1', 'true', 't'])
        outf = kw.get('output')
        if outf not in output_list:
            outf = output_list[0]
        output = self.temporary_path(fname=sample.name+'_stats.'+outf)
        if outf == 'txt':
            out = open(output,"w")
        else:
            out = {}
        if by_chrom:
            chromlist = sample.chrmeta.keys()
        else:
            chromlist = [None]
        for chrom in chromlist:
            if outf == 'txt':
                if chrom:
                    out.write("Chromosome %s\n--------------------\n"%chrom)
                stats(sample,out=out,selection=chrom)
            else:
                out[chrom] = {}
                stats(sample,out=out[chrom],selection=chrom)
            if outf == 'txt' and chrom:
                out.write("\n--------------------\n")
        if outf == 'txt':
            out.close()
            self.new_file(output, 'stats')
        else:
            self._plot_pdf(output,out,sample.name)
            self.new_file(output, 'pdf')
        return self.display_time()

