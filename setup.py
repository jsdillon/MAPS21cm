from distutils.core import setup
import glob

__version__ = '0.0.1'

setup_args = {
    'name': 'Joint_Mapmaking_Power_Spectrum_Pipeline',
    'author': 'Josh Dillon',
    'license': 'Creative Commons Attribution-Noncommercial-Share Alike license',
    'author_email'='jsdillon@berkeley.edu',
    #'package_dir': {'uvdata':'src'},
    #'packages': ['uvdata'],
    #'scripts': glob.glob('scripts/*'),
    'version': __version__
}

if __name__== '__main__':
    apply(setup, (), setup_args)