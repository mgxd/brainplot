import os
from setuptools import setup

version = None
with open (os.path.join('brainplot', '__init__.py'), 'r') as fp:
    for line in (line.strip() for line in fp):
        if line.startswith('__version__'):
            version = line.split('=')[-1].strip().strip("''")
            break
if version is None:
    raise RuntimeError('Version information cannot be found')

DISTNAME = 'brainplot'
DOWNLOAD_URL = 'https://github.com/mgxd/brainplot'
MAINTAINER = 'Mathias Goncalves'
MAINTAINER_EMAIL = 'mathiasg@mit.edu'

# TODO: gather deps
if __name__ == '__main__':
    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          version=version)
