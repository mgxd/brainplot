import os
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext

VERSION = None
with open (os.path.join('brainplot', '__init__.py'), 'r') as fp:
    for line in (line.strip() for line in fp):
        if line.startswith('__version__'):
            VERSION = line.split('=')[-1].strip().strip("''")
            break
if VERSION is None:
    raise RuntimeError('Version information cannot be found')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

DISTNAME = 'brainplot'
URL = 'https://github.com/mgxd/brainplot'
DOWNLOAD_URL = 'https://github.com/mgxd/brainplot/archive/master.zip'
MAINTAINER = 'Mathias Goncalves'
MAINTAINER_EMAIL = 'mathiasg@mit.edu'

# until nibabel releases
links = ['https://github.com/nipy/nibabel/tarball/master#egg=nibabel-2.2.0dev']

if __name__ == '__main__':
    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          url=URL,
          download_url=DOWNLOAD_URL,
          version=VERSION,
          packages=find_packages(),
          install_requires=requirements,
          dependency_links=links,
          entry_points={'console_scripts':
                        ['brainplot=brainplot.plotting:main']})
