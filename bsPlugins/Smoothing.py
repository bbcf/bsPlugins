from bsPlugins import *
from bbcflib.gfminer import stream as gm_stream
from bbcflib.track import track
from bbcflib import genrep

size_def = 11
step_def = 1

meta = {'version': "1.0.0",
            'author': "BBCF",
            'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'track', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'window_size', 'type': 'int', 'required': True},
                 {'id': 'window_step', 'type': 'int', 'required': True},
                 {'id': 'by_feature', 'type': 'boolean'},
                 {'id': 'format', 'type': 'list'}]
out_parameters = [{'id': 'smoothed_track', 'type': 'track'}]


class SmoothingForm(BaseForm):
    track = twb.BsFileField(
        label='Signal: ',
        help_text='Select signal file (e.g. bedgraph)',
        validator=twb.BsFileFieldValidator(required=True))
    assembly = twf.SingleSelectField(
        label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    window_size = twf.TextField(
        label='Window size: ',
        validator=twc.IntValidator(required=True),
        value=size_def,
        help_text='Size of the sliding window')
    window_step = twf.TextField(
        label='Window step: ',
        validator=twc.IntValidator(required=True),
        value=step_def,
        help_text='Size of steps between windows')
    by_feature = twf.CheckBox(
        label='Window size in features (not basepairs): ',
        value=False,
        help_text='Will count size and step parameters in number of features, not in basepairs')
    format = twf.SingleSelectField(
        label='Output format: ',
        options=['sql','bedGraph','wig','bigWig','sga'],
        prompt_text=None,
        help_text='Format of the output file', )
    submit = twf.SubmitButton(id="submit", value="Submit")

class SmoothingPlugin(BasePlugin):
    """Applies a moving average transformation to smooth the signal of a quantitative track. """

    info = {
        'title': 'Window smoothing',
        'description': __doc__,
        'path': ['Signal', 'Smoothing'],
        'output': SmoothingForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        tinput = track(kw.get('track'), chrmeta=kw.get('assembly') or None)
        outformat = kw.get('format',tinput.format)
        wsize = int(kw.get('window_size', size_def))
        wstep = int(kw.get('window_step', step_def))
        featurewise = kw.get('by_feature', False)
        if isinstance(featurewise, basestring):
            featurewise = (featurewise.lower() in ['1', 'true', 't'])
        output = self.temporary_path(fname=tinput.name+'_smoothed', ext=outformat)
        if featurewise:
            outfields = tinput.fields
            datatype = "qualitative"
        else:
            outfields = ["chr","start", "end", "score"]
            datatype = "quantitative"
        tout = track(output, format=outformat, fields=outfields, chrmeta=tinput.chrmeta, info={'datatype': datatype})
        for chrom in tout.chrmeta.keys():
            s = gm_stream.window_smoothing(
                    tinput.read(selection=chrom, fields=outfields),
                    window_size=wsize, step_size=wstep,
                    featurewise=featurewise)
            tout.write(s, chrom=chrom)
        tout.close()
        self.new_file(output, 'smoothed_track')
        return self.display_time()
