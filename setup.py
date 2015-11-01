import glob
from setuptools import setup, Command
import os

__version__ = '0.0.1'

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.egg-info')

setup_args = {
    'name': 'maps21cm',
    'author': 'Josh Dillon',
    'license': 'Creative Commons Attribution-Noncommercial-Share Alike license',
    'author_email': 'jsdillon@berkeley.edu',
    'package_dir': {'maps21cm':'Source'},
    'packages': ['maps21cm'],
#    'scripts': glob.glob('Scripts/*'),
    'version': __version__,
    'cmdclass':{'clean': CleanCommand}
}


if __name__== '__main__':
    apply(setup, (), setup_args)