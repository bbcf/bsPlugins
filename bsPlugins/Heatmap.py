from bsPlugins import *
from bbcflib.gfminer.figure import heatmap
from bbcflib.track import track
from numpy import array, log2, median, abs
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri

nb_colors_def = 10

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'track'},
                 {'id': 'log', 'type': 'boolean'},
                 {'id': 'list', 'type': 'boolean'},
                 {'id': 'nb_colors', 'type': 'int'}]
out_parameters = [{'id': 'Heatmap', 'type': 'pdf'},
                  {'id': 'List', 'type': 'file'}]


class HeatmapForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
                                help_text='Enter a table in which the first column contains the IDs and the others the values',
                                validator=twb.BsFileFieldValidator(required=True))
    log = twf.CheckBox(label='Log: ',
                                value=False,
                                help_text='Take the log2(1+x) of each value x')
    list = twf.CheckBox(label='List: ',
                                value=False,
                                help_text='Number the rows in the heatmap and make a file with the corresponding IDs')
    nb_colors = twf.TextField(label='Number of colors: ',
                                validator=twc.IntValidator(required=False),
                                value=nb_colors_def,
                                help_text='Number of colors between blue and red (default: 10)')
    submit = twf.SubmitButton(id="submit", value="Plot")


class HeatmapPlugin(BasePlugin):
    """Creates a heatmap of the table using *rows* as row labels and *columns* as column labels. The values are assumed to be equal to log2(raw data). If not, you have to select the option *Log*. Selecting the option *List* will print numbers beside the rows in the heatmap and make a file with the corresponding IDs."""
    info = {
        'title': 'Make heatmap of the table',
        'description': __doc__,
        'path': ['Graphics', 'Heatmap'],
        'output': HeatmapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }


    def __call__(self, **kw):

        def heatmapID(M,names,output1,output2,**kwargs):
            robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
            robjects.r.assign('names',numpy2ri.numpy2ri(names))
            if kwargs.get('nb_colors') is not None:
                robjects.r("ncol=%i" %kwargs['nb_colors'])
            else:
                robjects.r("ncol=10")

            robjects.r("""
              library(gplots)
              library(RColorBrewer)
              ymin=floor(min(Mdata,na.rm=T))
              ymax=ceiling(max(Mdata,na.rm=T))
              myBreaks = seq(ymin,ymax,length.out=ncol+1)
              myColors = rev(colorRampPalette(brewer.pal(ncol,"RdYlBu"))(ncol))
        #      myCor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
              data = as.matrix(Mdata)
              h=heatmap.2(data,
                        col=myColors, trace="none", breaks=myBreaks, #distfun=myCor,
                        na.rm=TRUE, dendrogram='both', density.info='none')
              nb_names=ceiling(nrow(data)/30)
              rownames=rep('',nrow(data))
              rownames[h$rowInd[seq(1,length(h$rowInd),nb_names)]]=match(seq(1,nrow(data)),h$rowInd)[h$rowInd[seq(1,length(h$rowInd),nb_names)]]
              n=c(); for (k in 1:length(names)){n[k]=names[[k]]}
              listnames=n[h$rowInd]""")
            output1.append(array(robjects.r("rownames")))
            output2.append(array(robjects.r("listnames")))

        table = open(kw.get('table'))
        col_labels = table.readline().strip("\n").strip("\r").split("\t")[1:]
        nb_colors = int(kw.get('nb_colors') or nb_colors_def)
        names = []; values = []
        for line in table:
            newline = line.strip("\n").strip("\r").split("\t")
            names.append(newline[0])
            values.append([float(newline[i]) for i in range(1,len(newline))])
        logscale = kw.get('log',False)
        if isinstance(logscale, basestring):
            logscale = (logscale.lower() in ['1', 'true', 't','on'])
        if logscale:
            for i in range(0,len(values)):
                values[i] = log2([x+1 for x in values[i]])
        values = array(values)
        med = median(values,axis=0)
        mad = median(abs(values-med),axis=0)
        values = ((values-med)/mad)
        make_list = kw.get('list',False)
        if isinstance(make_list, basestring):
            make_list = (make_list.lower() in ['1', 'true', 't','on'])
        if make_list:
            v1 = []; v2 = []
            heatmapID(array(values),names,v1,v2,nb_colors=nb_colors)
            names = v1[0].tolist()
            listnames = v2[0].tolist()
            List_ID = self.temporary_path(fname='List_ID.txt')
            with open(List_ID,"w") as L:
                for position, name in enumerate(listnames):
                    L.write(str(position+1)+"\t"+name+"\n")
            self.new_file(List_ID, 'List')
        title = table.name
        pdf = self.temporary_path(fname='Heatmap.pdf')
        heatmap(array(values), output=pdf, rows=names, columns=col_labels, main=title, nb_colors=nb_colors)
        self.new_file(pdf, 'Heatmap')
        return self.display_time()
