from bsPlugins import *
from bbcflib.btrack import track
import tarfile, sys

class FilteringForm(BaseForm):
    class TrackMulti(twb.BsMultiple):
        label='Tracks: '
        tracks = twb.BsFileField(label=' ',
                                 help_text='Select files (e.g. bedgraph)',
                                 validator=twb.BsFileFieldValidator(required=True))
    minscore = twf.TextField(label='Minimum score: ',
                             validator=twb.FloatValidator(required=False))
    maxscore = twf.TextField(label='Maximum score: ',
                             validator=twb.FloatValidator(required=False))
    minlength = twf.TextField(label='Minimum length: ',
                              validator=twc.IntValidator(required=False))
    maxlength = twf.TextField(label='Maximum length: ',
                              validator=twc.IntValidator(required=False))
    chrom = twf.TextField(label='Chromosome: ',
                          help_text='Comma-separated list of chromosome names',
                          validator=twc.Validator(required=False))
    submit = twf.SubmitButton(id="submit", value="Filter")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [{'id': 'tracks', 'type': 'track', 'required': True, 'multiple':'TrackMulti'},
                 {'id': 'minscore', 'type': 'float'},
                 {'id': 'maxscore', 'type': 'float'},
                 {'id': 'minlength', 'type': 'int'},
                 {'id': 'maxlength', 'type': 'int'},
                 {'id': 'chrom', 'type': 'text'}]
out_parameters = [{'id': 'output', 'type': 'track'},
                  {'id': 'archive', 'type': 'file'}]

class FilteringPlugin(BasePlugin):
    """Select features from a track passing a filter."""
    info = {
        'title': 'Apply a filter to a track',
        'description': __doc__,
        'path': ['Files', 'Filtering'],
        'output': FilteringForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta}

    def __call__(self, **kw):
        tracks = kw['TrackMulti']['tracks']
        if not isinstance(tracks, list): tracks = [tracks]
        minscore = kw.get('minscore')
        maxscore = kw.get('maxscore')
        minlength = kw.get('minlength')
        maxlength = kw.get('maxlength')
        selection = [{'chr': c} for c in kw.get('chrom','').split(',')]
        if minscore or maxscore:
            if not minscore:
                minscore = -sys.maxint
            if not maxscore:
                maxscore = sys.maxint
            if minscore > maxscore: 
                raise ValueError("Empty range: %f:%f" %(minscore,maxscore))
            for s in selection:
                s['score'] = (float(minscore),float(maxscore))
        if minlength or maxlength:
            minlength = int(minlength or 0)
            maxlength = int(maxlength or sys.maxint)
            if minlength > maxlength: 
                raise ValueError("Empty range: %i:%i" %(minlength,maxlength))
            for s in selection:
                s['length'] = (minlength,maxlength)
        outtracks =[]
        for tin in [track(t) for t in tracks]:
            outname = self.temporary_path(tin.name+"_filtered."+tin.format)
            tout = track(outname)
            outtracks.append(outname)
            outstream = tin.read(selection=selection)
            tout.write(outstream)
            tout.close()

        if len(outtracks) > 1:
            tar_name = self.temporary_path('Filtered_tracks.tgz')
            tar = tarfile.open(tar_name, "w:gz")
            [tar.add(f) for f in outtracks]
            tar.close()
            self.new_file(tar_name, 'archive')
        else:
            self.new_file(outtracks[0], 'output')
        return self.display_time()
