from bsPlugins import *
from bbcflib.bFlatMajor import stream as gm_stream
from bbcflib import btrack as track
from bbcflib import genrep

size_def = 11
step_def = 1

class SmoothingForm(BaseForm):
    track = twf.FileField(label_text='Signal: ',
        help_text='Select signal file (e.g. bedgraph)',
        validator=twf.FileValidator(required=True))
    assembly = twf.SingleSelectField(label_text='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    window_size = twf.TextField(label_text='Window size: ',
        validator=twc.IntValidator(required=True),
        value=size_def,
        help_text='Size of window')
    window_step = twf.TextField(label_text='Window step: ',
        validator=twc.IntValidator(required=True),
        value=step_def,
        help_text='Size of steps between windows')
    by_feature = twf.CheckBox(label_text='Window size in features (not basepairs): ',
        value=False,
        help_text='Will count size and step parameters in number of features, not in basepairs')
    submit = twf.SubmitButton(id="submit", value="Smooth")

meta = {'version': "1.0.0",
            'author': "BBCF",
            'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'track', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'window_size', 'type': 'int', 'required': True},
                 {'id': 'window_step', 'type': 'int', 'required': True},
                 {'id': 'by_feature', 'type': 'boolean'}]
out_parameters = [{'id': 'smoothed_track', 'type': 'track'}]


class SmoothingPlugin(OperationPlugin):

    info = {
        'title': 'Window smoothing',
        'description': 'Window smoothing',
        'path': ['Signal', 'Smoothing'],
        'output': SmoothingForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        tinput = track.track(kw.get('track'), chrmeta=kw.get('assembly') or None)
        wsize = int(kw.get('window_size', size_def))
        wstep = int(kw.get('window_step', step_def))
        featurewise = kw.get('by_feature', False)
        if isinstance(featurewise, basestring):
            featurewise = (featurewise.lower() in ['1', 'true', 't'])
        output = self.temporary_path(fname='smoothed_track', ext='sql')
        if featurewise:
            outfields = tinput.fields
            datatype = "qualitative"
        else:
            outfields = ["start", "end", "score"]
            datatype = "quantitative"
        tout = track.track(output, fields=outfields,
                           chrmeta=tinput.chrmeta,
                           info={'datatype': datatype})
        for chrom in tout.chrmeta.keys():
            tout.write(gm_stream.window_smoothing(
                    tinput.read(selection=chrom, fields=outfields),
                    window_size=wsize, step_size=wstep,
                    featurewise=featurewise), chrom=chrom)
        tout.close()
        self.new_file(output, 'smoothed_track')
        return self.display_time()
