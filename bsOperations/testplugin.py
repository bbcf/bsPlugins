import os, sys


class TestPlugin():

    def __init__(self, plugin_name):
        self.plugin_name = plugin_name
        plugins = [f[:-3] for f in os.listdir('bs/plugins') if f.endswith('.py')]
        if plugin_name not in plugins:
            raise Exception("Plugin %s does not exist" % plugin_name)


    def process(self, **kw):
        mod    = __import__('bs.plugins', fromlist=[self.plugin_name])
        clz = getattr(getattr(mod, self.plugin_name), self.plugin_name)
        c = clz()
        c._pre_process('Testing service')
        c.process(**kw)
        print '[x] ok'




if __name__ == '__main__':
    t = TestPlugin('Test')
    t.process(**{})