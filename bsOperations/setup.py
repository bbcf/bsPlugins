from distutils.core import setup

install_requires=[
    "tw2.forms",
    "tw2.dynforms",
    "yapsy"
    ]

setup(
    name='plugins',
    version='0.1',
    description='Plugins for BioScript',
    author='Yohan Jarosz',
    author_email='yohan.jarosz@epfl.ch',
    url='http://bbcf.epfl.ch/bs',
    install_requires=install_requires,
    packages=['bs','bs.lib','bs.lib.plugins','bs.plugins']
)
