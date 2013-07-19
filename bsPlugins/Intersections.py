from bsPlugins import *
from itertools import combinations
import os
import tarfile


class IntersectionsForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signals: '
        signals = twb.BsFileField(label=' ',
            help_text='Select signal files (e.g. bedgraph)',
            validator=twb.BsFileFieldValidator(required=True))
    columns = twf.TextField(label='Column(s): ',
        prompt_text='0',
        value = 0,
        help_text='Column(s) number (1-based).')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': 'SigMulti', 'required': True},
                 {'id': 'column', 'type': 'text'}]
out_parameters = [{'id': 'intersections', 'type': 'track'}]


class IntersectionsPlugin(BasePlugin):
    """Returns the elements that are common to a set of text files,
for instance the list of genes common to several lists of genes or annotation files.<br /><br />

In the case when more that two files are given, all possible combinations of intersections
are performed (2-by-2, 3-by-3, etc.), in the manner of a Venn diagram.
If the elements to intersect are not in the first column, one can specify the column to consider
by its index (first column is 1), or a list of column indices separated by commas if the columns
are different between files (e.g. '5,1,2' for column 5 of file 1, column 1 of file 2, etc.).<br /><br />

Since the number of comparisons is approximately 2^(number of files), it is unadvised to compare more
that a dozen of files (10 input files -> 2^10-11=1013 comparisons).

The output is a compressed folder containing a summary file and a sub-folder with all the possible
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
        return output

    def __call__(self,**kw):
        files_list = kw['signals']
        column = kw['column']
        output = self.temporary_path(fname='intersections.')
        self.compare(files_list, output, column)
        # compress
        output_targz = self.temporary_path(fname=output+'tar.gz')
        tar = tarfile.open(output_targz, 'w:gz')
        tar.add(output)
        tar.close()
        # send
        self.new_file(output+'.tar.gz', 'intersections')
        return self.display_time()

