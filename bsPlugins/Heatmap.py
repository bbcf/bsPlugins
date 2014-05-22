from bsPlugins import *
from bbcflib.gfminer.figure import heatmap
from bbcflib.track import track
from numpy import array, log2, median, abs

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

    def heatmapID(M,names,output=None,**kwargs):
        """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
        If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""
        robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
        robjects.r.assign('names',numpy2ri.numpy2ri(names))
        if orderCols and orderRows:
            plotopt += ",dendrogram='both',lhei=c(2,10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(0,2,0,0,3,1,0,4),ncol=2)"
        elif orderCols:
            plotopt += ",Rowv=F,dendrogram='column',lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)"
        elif orderRows:
            plotopt += ",Colv=F,dendrogram='row',lhei=c(10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(2,0,3,1,0,4),ncol=2)"
        else:
            plotopt += ",Colv=F,Rowv=F,dendrogram='none',lhei=c(10,1,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(1,2,3,4),ncol=1)"
        if kwargs.get('ymin') is not None:
            robjects.r("ymin=%f" %float(kwargs['ymin']))
        else:
            robjects.r("ymin=floor(min(Mdata,na.rm=T))")
        if kwargs.get('ymax') is not None:
            robjects.r("ymax=%f" %float(kwargs['ymax']))
        else:
            robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
        if kwargs.get('nb_colors') is not None:
            robjects.r("ncol=%i" %kwargs['nb_colors'])
        else:
            robjects.r("ncol=10")

        robjects.r("""
          library(gplots)
          library(RColorBrewer)
          myBreaks = seq(ymin,ymax,length.out=ncol+1)
          myColors = rev(colorRampPalette(brewer.pal(ncol,"RdYlBu"))(ncol))
    #      myCor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
          data = as.matrix(Mdata)
          h=heatmap.2(data,
                    col=myColors, trace="none", breaks=myBreaks, #distfun=myCor,
                    na.rm=TRUE, density.info='none'%s)
          names = rep('',nrow(data))
          names[h$rowInd[seq(1,length(h$rowInd),200)]]= match(seq(1,nrow(dmat)),h$rowInd)[h$rowInd[seq(1,length(h$rowInd),200)]]
          genes = gene_names[h$rowInd]
          output = h$rowInd""" %plotopt)
        return output

    def __call__(self, **kw):
        table = open(kw.get('table'))
        names = []; values = []
        for line in table:
            newline = line.strip("\n").strip("\r").split("\t")
            name = newline[0]
            names.append(name)
            value = [float(newline[i]) for i in range(1,len(newline))]
            values.append(value)
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
            h=heatmapID(array(values),names)
            List_ID = self.temporary_path(fname='List_ID.txt')
            with open(List_ID,"w") as L:
                L.write("\t".join(names)+"\n")
            self.new_file(List_ID, 'List')
        nb_colors = int(kw.get('nb_colors') or nb_colors_def)
        title = "List of genes"
        pdf = self.temporary_path(fname='Heatmap.pdf')
        heatmap(array(values), output=pdf, rows=names, main=title, nb_colors=nb_colors)
        self.new_file(pdf, 'Heatmap')
        return self.display_time()
