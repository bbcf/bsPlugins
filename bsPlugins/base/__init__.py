import hashlib
import os
import wordlist
import time
import string
import random
import tempfile


random_name = lambda x: ''.join(random.choice(string.ascii_lowercase + string.digits) for i in xrange(x))
TMP_DIR = '.'


class OperationPlugin(object):

    bs_plugin = 'bs-operation'

    def __init__(self):

        if not hasattr(self, 'info'):
            self.info = {}

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
        self.start_time = None
        self.end_time = None

    def _start_timer(self):
        self.start_time = time.time()

    def _end_timer(self):
        self.end_time = time.time()

    def time(self):
        if self.end_time is None:
            self._end_timer()
        return self.end_time - self.start_time

    def display_time(self):
        return 'Time elapsed %0.3fs.' % self.time()

    def unique_id(self):
        '''
        It's an unique identifier for your plugin.
        Do not override
        '''
        if self.uid is None:
            tohash = str(self.__class__) + str(self.title) + str(self.in_parameters) + str(self.out_parameters)
            if 'version' in self.meta:
                tohash += self.meta['version']
            self.uid = hashlib.sha1(tohash).hexdigest()
        return self.uid

    def new_file(self, fpath, fparam):
        """
        Append a file to the result
        """
        added = False
        for p in self.out_parameters:
            if p.get('id') == fparam:
                added = True
                ftype = p.get('type')
                self.output_files.append([fpath, ftype])
                break
        if not added:
            raise Exception("You're trying to add an output file but you did not specify it as an output in the plugin : you gave %s and the parameters allowed are %s" % (fparam, ', '.join([o['id'] for o in self.out_parameters])))

    def in_params_typeof(self, typeof):
        return [param for param in self.in_parameters if wordlist.is_of_type(param.get('type'), typeof)]

    def temporary_path(self, fname=None, ext=None):
        """
        Get a temporary path to write a file.
        The file will be automatically deleted at the end of the plugin process.
        :param fname: the file name
        :return: a path
        """
        tmp_dir = tempfile.mkdtemp(dir=TMP_DIR)
        if fname is None or fname == '':
            fname = random_name(6)
        if ext is not None:
            if ext.startswith('.'):
                fname += ext
            else:
                fname = '%s.%s' % (fname, ext)

        fpath = os.path.join(tmp_dir, fname)
        self.tmp_files.append(tmp_dir)
        return fpath
