import hashlib
import os
import wordlist
import time
import string
import random
import tempfile


random_name = lambda x: ''.join(random.choice(string.ascii_lowercase + string.digits) for i in xrange(x))
TMP_DIR = '.'


class BasePlugin(object):

    bs_plugin = 'bs-operation'

    def __init__(self):

        if not hasattr(self, 'info'): self.info = {}
        self.title = self.info.get('title', '')
        self.description = self.info.get('description', '')
        self.path = self.info.get('path', None)
        self.output = self.info.get('output', '')
        self.in_parameters = self.info.get('in', [])
        self.out_parameters = self.info.get('out', [])
        self.meta = self.info.get('meta', '')
        self.deprecated = self.info.get('deprecated', False)
        self.uid = None
        self.service = None
        self.output_files = []
        self.tmp_files = []
        self.start_time = 0
        self.end_time = 0
        self.is_debug = False
        self.debug_stack = []

        
    def _start_timer(self):
        self.start_time = time.time()

    def display_time(self):
        if self.end_time is None: self.end_time = time.time()
        t = max(0.0,self.end_time-self.start_time)
        return 'Time elapsed %0.3fs.' % t

    def html_doc_link(self):
        """
        The default link to the plugin documentation. You can override this function
        to reference the HTML page you want.
        """
        return 'http://bbcf.epfl.ch/bsplugins/content/bsPlugins.html#bsPlugins.%s.%s' % (self.__module__, self.__class__.__name__)

    def html_source_code_link(self):
        """
        The default link to the plugin source code. You can override this function
        to reference the HTML page you want.
        """
        return 'https://github.com/bbcf/bsPlugins/tree/master/bsPlugins/%s.py' % (self.__module__.split('.')[-1])

    def unique_id(self):
        '''
        It's a unique identifier for your plugin.
        Do not override
        '''
        if self.uid is None:
            tohash = str(self.__class__) + str(self.title) + str(self.in_parameters) + str(self.out_parameters)
            if 'version' in self.meta: tohash += self.meta['version']
            self.uid = hashlib.sha1(tohash).hexdigest()
        return self.uid

    def new_file(self, fpath, fparam):
        """
        Append a file to self.output_files .
        """
        added = False
        for p in self.out_parameters:
            if p.get('id') == fparam:
                added = True
                ftype = p.get('type')
                self.output_files.append([fpath, ftype])
                break
        if not added:
            raise Exception("You did not specify %s as a plugin output, only %s" % (fparam, self.out_parameters))
        print "%s (%s): %s" %(fparam, ftype, fpath)

    def in_params_typeof(self, typeof):
        return [param for param in self.in_parameters if wordlist.is_of_type(param.get('type'), typeof)]

    def temporary_path(self, fname=None, ext=None):
        """
        Create a temporary folder; generate and return a file name inside of this folder
        (the file is not created yet).
        The folder will be automatically deleted at the end of the Bioscript process.
        :param fname: the file name (a random string by default).
        :param ext: the file extension.
        :return: the absolute path to the file.
        """
        tmp_dir = tempfile.mkdtemp(dir=TMP_DIR)
        if fname is None or fname == '': fname = random_name(6)
        if ext is not None:
            if ext.startswith('.'): fname += ext
            else:                   fname = "%s.%s"%(fname, ext)
        if hasattr(self, 'outprefix') and self.outprefix: fname = self.outprefix+fname
        fpath = os.path.join(tmp_dir, fname)
        fpath = os.path.abspath(fpath)
        self.tmp_files.append(tmp_dir)
        return fpath


def test_form(clz, port=8000):
    """
    :param clz: Form class to test.
    Example ::

    >>> import bsPlugins
    >>> from bsPlugins.Smoothing import SmoothingForm
    >>> bsPlugins.base.test_form(SmoothingForm, port=8080)
    """
    import tw2.forms
    import tw2.devtools

    class Index(tw2.forms.FormPage):
        title = 'Test bs plugin form'
        child = clz

    tw2.devtools.dev_server(port=port)
