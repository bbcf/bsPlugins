from bsPlugins import *
from bbcflib import genrep
from bbcflib.btrack import track
from bbcflib.bFlatMajor import common
import math
import tw2.forms as twf


class OperationForm(BaseForm):
    class SigMulti(Multi):
        label='Signals: '
        track = twf.FileField(label=' ',
        help_text='Select files (e.g. bedgraph)',
        validator=twf.FileValidator(required=True))
    function =  twf.SingleSelectField(label='Operation: ',
        options=["log2","log10","sqrt"],
        validator=twc.Validator(required=True),
        help_text='Select a method, by default: log2')
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bedgraph","bigwig","wig"],
        validator=twc.Validator(required=False),
        help_text='Output file(s) format, by default: same format as input file(s) format(s)')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [{'id': 'track', 'type': 'track', 'required': True},    # input parameters
                {'id': 'assembly', 'type': 'assembly', 'required': True},
                {'id': 'function', 'type': 'function', 'required': True},
                {'id': 'format', 'type': 'format'}]
out_parameters = [{'id': 'output', 'type': 'file'}]    # output parameters


class OperationPlugin(BasePlugin):
    description = """Apply a numeric transformation to the track scores - such as logarithm or square root."""
    info = {
        'title': 'Operation',
        'description': description,
        'path': ['Signal', 'Operation'],
        'output': OperationForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        def method(x):    # the function that is applied to the track score
            if kw.get('function')=="log2":
                return math.log( x,2) ;
            elif kw.get('function')=="log10":
                return math.log10( x) ;
            elif kw.get('function')=="sqrt":
                return math.sqrt(x) ;
            else:
                return math.log(x) ;

        assembly_id = kw.get('assembly') #or None
        assembly = genrep.Assembly(assembly_id)
        length=len( kw.get('track') )    # number of tracks
        for i in range( length  ) :
            tname = kw['track'][i]
            tinput = track(tname, chrmeta=kw.get('assembly'))
            (filepath, filename) = os.path.split(tname)    # ( path, name of one track)
            (shortname, extension) = os.path.splitext(filename)    # (name of one track (without path) , extension)
            modif= kw['function']
            if not isinstance(tinput, (tuple,list)):
                tinput = [tinput]
            for trackFILE in tinput:
                if  "score" in trackFILE.fields:
                    if kw.get('format')=="":
                        #create a temporary directory, name a new file in this directory
                        output1_name = self.temporary_path( shortname+'_'+kw.get('assembly')+'_'+modif+str(extension))
                        t1_name = track(output1_name,chrmeta=assembly.chrmeta) #define an output track file in the temporary directory
                        #apply method function to the score field trackFILE
                        #create and write in the new file
                        t1_name.write(common.apply(trackFILE.read(),'score',method ) ,mode='write')
                        t1_name.close()
                        trackFILE.close()
                    if kw.get('format')!="":
                        output2_name = self.temporary_path( shortname +'_'+kw.get('assembly')+'_'+ modif +'.'+ kw.get('format') )
                        t2_name = track(output2_name,chrmeta=assembly.chrmeta)
                        t2_name.write(common.apply(trackFILE.read(),'score',method ) ,mode='write')
                        t2_name.close()
                        trackFILE.close()
