"""Still need to (g?)zip the output folder and  use track to read input files"""

from bsPlugins import *
from bbcflib.track import track
from itertools import combinations
import os
import zipfile


class IntersectionsForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signals: '
        signals = twb.BsFileField(label=' ',
                                help_text='Select signal files (e.g. bedgraph)',
                                validator=twb.BsFileFieldValidator(required=True))
    column = twf.TextField(label='Column: ',
        prompt_text='0',
        value = 0,
        help_text='Column number (1-based) or column name')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': 'SigMulti', 'required': True},
                 {'id': 'column', 'type': 'text'}]
out_parameters = [{'id': 'intersections', 'type': 'track'}]


class IntersectionsPlugin(BasePlugin):
    """Returns the elements that are common to a set of files,
for instance the list of genes common to several lists of genes or annotation files.<br /><br />

In the case when more that two files are given, all the possible combinations of intersections
are performed (2-by-2, 3-by-3, etc.), in the manner of a Venn diagram.
If the elements to intersect are not in the first column, one can specify the column to consider,
either by its column index (first column is 1) or by the field name.<br />
The output is a zipped folder containing a summary file and a sub-folder with all the possible
intersections, i.e. for each intersection one text file with the list of common elements.
    """
    info = {
        'title': 'Intersections',
        'description': __doc__,
        'path': ['Signal', 'Intersections'],
        'output': IntersectionsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def zip_file(self,target):
        target_dir = os.path.abspath(target)
        target_name = os.path.splitext(os.path.basename(target))[0]
        parent_dir = os.path.dirname(target)
        zipped = os.path.join(parent_dir,target_name+'.zip')
        zip = zipfile.ZipFile(zipped, 'w', zipfile.ZIP_DEFLATED)
        rootlen = len(target_dir) + 1
        for base, dirs, files in os.walk(parent_dir):
           for file in files:
              fn = os.path.join(base, file)
              zip.write(fn, fn[rootlen:])
        return zipped

    def intersect(self, files_list, idx=0):
        common = set()
        for n,f in enumerate(files_list):
            col = set(line.strip().split()[idx] for line in open(f))
            if n==0:
                common = col
            else:
                common &= col
        return common

    def compare(self, files_list, output, idx=0):
        if not os.path.exists(output):
            os.mkdir(output)
        counts = {}
        counts = open(os.path.join(output,"summary.txt"), 'wb')
        counts.write("# Legend:\n")
        for i,f in enumerate(files_list):
            counts.write("%d\t%s\n" % (i,f))
        counts.write("\n### Files\tnb_elements\n")
        counts.write("\n# Self\n\n")
        for i,f in enumerate(files_list):
            counts.write("%d\t%d\n" % (i,len(open(f).readlines())))
        for k in range(2,len(files_list)+1):
            counts.write("\n# %d-by-%d\n\n" % (k,k))
            path = os.path.join(output,"%s-by-%s/" % (k,k))
            if not os.path.exists(path):
                os.mkdir(path)
            combs = combinations(range(len(files_list)), k)
            for cb in combs:
                names = sorted([str(x) for x in cb])
                name = "-".join(names)
                out = open(os.path.join(path,"%s.txt" % name), 'wb')
                common = self.intersect([files_list[i] for i in cb], idx)
                counts.write("%s\t%s\n" % (name,len(common)))
                for x in common:
                    out.write(x+'\n')
                out.close()
        counts.close()

    def __call__(self,**kw):
        files_list = kw['signals']
        column = kw['column']
        t = None # later
        try: column = int(column)
        except ValueError: column = t.fields.index(column)
        output = 'intersections'
        self.compare(files_list, output, column)
        #output = self.zip_file(output)
        self.new_file(output, 'intersections')
        return self.display_time()

