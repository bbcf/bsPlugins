from tempfile import NamedTemporaryFile

class BaseForm():
    pass

class OperationPlugin():

    def __init__(self):
        self.files={}

    def temporary_path(self, fname='', ext=''):
        tmpf = NamedTemporaryFile(prefix=fname, suffix="."+ext, 
                                  dir="./", delete=False)
        tmpf.close()
        return tmpf.name

    def new_file(self, fname, descr):
        self.files[fname] = descr
        print "%s: %s" %(fname,descr)

    def display_time(self):
        return 1


class Multi():
    pass
