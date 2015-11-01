from distutils.core import setup
import glob

__version__ = '0.0.1'

setup_args = {
    'name': 'maps21cm',
    'author': 'Josh Dillon',
    'license': 'Creative Commons Attribution-Noncommercial-Share Alike license',
    'author_email': 'jsdillon@berkeley.edu',
    'package_dir': {'maps21cm':'Source'},
    'packages': ['maps21cm'],
    'scripts': glob.glob('Scripts/*'),
    'version': __version__
}

if __name__== '__main__':
    apply(setup, (), setup_args)