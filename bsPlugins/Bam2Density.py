from bsPlugins import *
from bein import execution
from bbcflib.btrack import track, convert
from bbcflib.mapseq import bam_to_density


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'sample', 'type': 'bam', 'required': True},
                 {'id': 'control', 'type': 'bam'},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'normalization', 'type': 'int'},
                 {'id': 'merge_strands', 'type': 'int'},
                 {'id': 'read_extension', 'type': 'int'}]
out_parameters = [{'id': 'density_merged', 'type': 'track'},
                  {'id': 'density_fwd', 'type': 'track'},
                  {'id': 'density_rev', 'type': 'track'}]

__requires__ = ["pysam"]


class Bam2DensityForm(BaseForm):
    sample = twf.FileField(label='Test BAM: ',
                           help_text='Select main bam file',
                           validator=twf.FileValidator(required=True))
    control = twf.FileField(label='Control BAM: ',
                            help_text='Select control bam file to compute enrichment')
    format = twf.SingleSelectField(label='Output format: ',
                                   options=["sql", "bedGraph", "bigWig"],
                                   prompt_text=None,
                                   help_text='Format of the output file')
    normalization = twf.TextField(label='Normalization: ',
                                  validator=twc.IntValidator(required=False),
                                  help_text='Normalization factor, default is total number of reads')
    merge_strands = twf.TextField(label='Shift and merge strands: ',
                                  validator=twc.IntValidator(required=False),
                                  help_text='Enter shift value (in bp) if you want to merge strand-specific densities')
    read_extension = twf.TextField(label='Read extension: ',
                                   validator=twc.IntValidator(required=False),
                                   help_text='Enter read extension (in bp) to be applied when constructing densities')
    submit = twf.SubmitButton(id="submit", value='bam2density')


class Bam2DensityPlugin(BasePlugin):

    info = {
        'title': 'Bam2density',
        'description': 'Bam2density generates genome-wide densities from bam files',
        'path': ['Files', 'Bam2density'],
        'output': Bam2DensityForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        b2wargs = []
        control = None
        sample = kw.get("sample")
        assert os.path.exists(str(sample)), "Bam file not found: '%s'." % sample
        if kw.get('control'):
            control = kw['control']
            b2wargs = ["-c", str(control)]
            assert os.path.exists(str(control)), "Control file not found: '%s'." % control
            control = os.path.abspath(control)
        sample = os.path.abspath(sample)
        nreads = int(kw.get('normalization') or -1)
        bamfile = track(sample, format='bam')
        if nreads < 0:
            if control is None:
                nreads = len(set((t[4] for t in bamfile.read())))
            else:
                b2wargs += ["-r"]
        merge_strands = int(kw.get('merge_strands') or -1)
        read_extension = int(kw.get('read_extension') or -1)
        output = self.temporary_path(fname='density_')
        format = kw.get("format", "sql")
        with execution(None) as ex:
            files = bam_to_density(ex, sample, output,
                                   nreads=nreads, merge=merge_strands,
                                   read_extension=read_extension,
                                   sql=True, args=b2wargs)
        if merge_strands >= 0:
            suffixes = ["merged"]
        else:
            suffixes = ["fwd", "rev"]
        for n, x in enumerate(files):
            tsql = track(x, format='sql', fields=['start', 'end', 'score'],
                         chrmeta=bamfile.chrmeta,
                         info={'datatype': 'quantitative'})
            tsql.save()
            if format == "sql":
                outname = x
            else:
                outname = os.path.splitext(x)[0]+"."+format
                convert(x, outname, mode="overwrite")
            self.new_file(outname, 'density_'+suffixes[n])
        return self.display_time()
