from bsPlugins import *
from bein import execution
from bbcflib.btrack import track, convert
from bbcflib.mapseq import bam_to_density


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'sample', 'type': 'bam', 'required': True},
                 {'id': 'control', 'type': 'bam'},
                 {'id': 'format', 'type': 'txt'},
                 {'id': 'normalization', 'type': 'int'},
                 {'id': 'merge_strands', 'type': 'int'},
                 {'id': 'read_extension', 'type': 'int'}]
out_parameters = [{'id': 'density_merged', 'type': 'track'},
                  {'id': 'density_fwd', 'type': 'track'},
                  {'id': 'density_rev', 'type': 'track'}]

__requires__ = ["pysam"]


class Bam2DensityForm(BaseForm):
    sample = twf.FileField(label_text='Test BAM: ',
        help_text='Select main bam file',
        validator=twf.FileValidator(required=True))
    control = twf.FileField(label_text='Control BAM: ',
        help_text='Select control bam file to compute enrichment')
    format = twf.SingleSelectField(label='Output format: ',
        options=["wig","sql","bedGraph","bigWig"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    normalization = twf.TextField(label_text='Normalization: ',
        validator=twc.IntValidator(required=False),
        help_text='Normalization factor, default is total number of reads')
    merge_strands = twf.TextField(label_text='Shift and merge strands',
        validator=twc.IntValidator(required=False),
        help_text='Enter shift value (in bp) if you want to merge strand-specific densities')
    read_extension = twf.TextField(label_text='Read extension: ',
        validator=twc.IntValidator(required=False),
        help_text='Enter read extension (in bp) to be applied when constructing densities')
    submit = twf.SubmitButton(id="submit", value='bam2density')


class Bam2DensityPlugin(OperationPlugin):

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
        if kw.get('control'):
            control = kw['control']
            b2wargs = ["-c", str(control)]
        bamfile = track(kw['sample'], format='bam')
        nreads = int(kw.get('normalization') or -1)
        if nreads < 0:
            if control is None:
                nreads = len(set((t[4] for t in bamfile.read())))
            else:
                b2wargs += ["-r"]
        merge_strands = int(kw.get('merge_strands') or -1)
        read_extension = int(kw.get('read_extension') or -1)
        output = self.temporary_path(fname='density_')
        format = kw.get("format","sql")
        with execution(None) as ex:
            files = bam_to_density(ex, kw['sample'], output,
                                    nreads=nreads, merge=merge_strands,
                                    read_extension=read_extension,
                                    sql=(format=="sql"), args=b2wargs)
        if isinstance(files,list): # sql format, cf the mess with bam_to_density
            if merge_strands >= 0: suffixes = ["_merged"]
            else: suffixes = ["_fwd", "_rev"]
            for n, x in enumerate(files):
                if format == "sql":
                    tsql = track(x, format='sql', fields=['start', 'end', 'score'],
                                  chrmeta=bamfile.chrmeta, info={'datatype': 'quantitative'})
                    tsql.save()
                    self.new_file(x, 'density' + suffixes[n])
        elif format == "bedGraph":
            suffix = "_merged" if merge_strands >= 0 else ""
            self.new_file(x, 'density'+suffix+".bedGraph")
        else:
            suffix = "_merged" if merge_strands >= 0 else ""
            convert(x, (y,format), chrmeta=bamfile.chrmeta, info={'datatype': 'quantitative'}, mode="overwrite")
            self.new_file(y, 'density'+suffix+'.'+format)
        return 1
