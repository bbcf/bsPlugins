from bsPlugins import *
from bein import execution
from bbcflib.track import track, convert
from bbcflib.mapseq import bam_to_density
from bbcflib.gfminer.stream import merge_scores
import os

__requires__ = ["pysam"]


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'sample', 'type': 'bam', 'required': True, 'multiple': 'BamMulti'},
                 {'id': 'control', 'type': 'bam'},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'normalization', 'type': 'int'},
                 {'id': 'merge_strands', 'type': 'int'},
                 {'id': 'read_extension', 'type': 'int'}]
out_parameters = [{'id': 'density_merged', 'type': 'track'},
                  {'id': 'density_fwd', 'type': 'track'},
                  {'id': 'density_rev', 'type': 'track'}]


class Bam2DensityForm(BaseForm):
    class BamMulti(twb.BsMultiple):
        label='Test BAMs: '
        sample = twb.BsFileField(label=' ',
                                 help_text='Select main bam file(s)',
                                 validator=twb.BsFileFieldValidator(required=True))
    control = twb.BsFileField(label='Control BAM: ',
                              help_text='Select control bam file to compute enrichment',
                              validator=twb.BsFileFieldValidator(required=False))
    format = twf.SingleSelectField(label='Output format: ',
                                   options=["sql", "bedGraph", "bigWig"],
                                   prompt_text=None,
                                   help_text='Format of the output file')
    normalization = twf.TextField(label='Normalization: ',
                                  validator=twc.IntValidator(),
                                  help_text='Normalization factor, default is total number of reads')
    merge_strands = twf.TextField(label='Shift and merge strands: ',
                                  validator=twc.IntValidator(),
                                  value=-1,
                                  help_text='Shift value (in bp) if you want to merge strand-specific densities (will not merge if negative)')
    read_extension = twf.TextField(label='Read extension: ',
                                   validator=twc.IntValidator(),
                                   value=-1,
                                   help_text='Read extension (in bp) to be applied when constructing densities (will use read length if negative)')
    submit = twf.SubmitButton(id="submit", value='bam2density')


class Bam2DensityPlugin(BasePlugin):
    """From a BAM file, creates a track file of the read count/density along the whole genome,
in the chosen format.

Read counts are divided by 10^-7 times the normalization factor (which is total number of reads by default).
Positive and negative strand densities are generated and optionally merged (averaged) if a
shift value >=0 is given. The read extension is the number of basepairs a read will cover,
starting from its most 5' position (e.g. with a read extension of 1, only the starting position of
each alignment will be considered, default is read length).
"""
    info = {
        'title': 'Genome-wide reads density from BAM',
        'description': __doc__,
        'path': ['Files', 'Bam2density'],
        'output': Bam2DensityForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        b2wargs = []
        control = None
        samples = kw.get('BamMulti',{}).get('sample', [])
        if not isinstance(samples, list): samples = [samples]
        samples = [os.path.abspath(s) for s in samples if os.path.exists(s)]
        if kw.get('control'):
            control = kw['control']
            b2wargs = ["-c", str(control)]
            assert os.path.exists(str(control)), "Control file not found: '%s'." % control
            control = os.path.abspath(control)
        try: 
            nreads = int(kw.get('normalization'))
        except ValueError:
            nreads = -1
        bamfiles = [track(s, format='bam') for s in samples]
        if nreads < 0:
            _nreads = [0]*len(samples)
            if control is not None:
                b2wargs += ["-r"]
        else:
            _nreads = [nreads for s in samples]
        try: 
            merge_strands = int(kw.get('merge_strands'))
        except ValueError:
            merge_strands = -1
        try: 
            read_extension = int(kw.get('read_extension'))
        except ValueError:
            read_extension = -1
        output = [self.temporary_path(fname=b.name+'_density_') for b in bamfiles]
        format = kw.get("format", "sql")
        with execution(None) as ex:
            files = [bam_to_density( ex, s, output[n], nreads=_nreads[n], 
                                     merge=merge_strands,
                                     read_extension=read_extension,
                                     sql=True, args=b2wargs ) 
                     for n,s in enumerate(samples)]
        info = {'datatype': 'quantitative',
                'read_extension': read_extension}
        if merge_strands >= 0:
            suffixes = ["merged"]
            info['shift'] = merge_strands
        else:
            suffixes = ["fwd", "rev"]
        chrmeta = bamfiles[0].chrmeta
        for suf in suffixes:
            all_s_files = [x for y in files for x in y if x.endswith(suf+".sql")]
            if len(all_s_files) > 1:
                x = self.temporary_path(fname="Density_average_"+suf+".sql")
                tsql = track( x, fields=['start', 'end', 'score'],
                              chrmeta=chrmeta, info={'datatype': 'quantitative'} )
                insql = []
                for f in all_s_files:
                    t = track(f, format='sql', chrmeta=chrmeta)
                    t.save()
                    insql.append(t)
                for c in tsql.chrmeta:
                    tsql.write(merge_scores([t.read(c) for t in insql]),chrom=c)
            else:
                x = all_s_files[0]
                tsql = track( x, format='sql', fields=['start', 'end', 'score'],
                              chrmeta=chrmeta, info=info )
                tsql.save()
            if format in [None,"sql"]:
                outname = x
            else:
                outname = os.path.splitext(x)[0]+"."+format
                convert(x, outname, mode="overwrite")
            self.new_file(outname, 'density_'+suf)
        return self.display_time()
