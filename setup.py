import os.path as op
from setuptools import findall, setup, find_packages

thispath = op.dirname(__file__)
ldict = locals()

# Get version and release info, which is all stored in heudiconv/info.py
info_file = op.join(thispath, 'brainplot', 'info.py')
with open(info_file) as infofile:
    exec(infofile.read(), globals(), ldict)

SETUP_REQUIRES = ['numpy']

setup(name=ldict['__packagename__'],
      maintainer=ldict['MAINTAINER'],
      maintainer_email=ldict['MAINTAINER_EMAIL'],
      url=ldict['__url__'],
      version=ldict['__version__'],
      packages=find_packages(),
      install_requires=ldict['REQUIRES'],
      setup_requires=SETUP_REQUIRES,
      tests_requires=ldict['TESTS_REQUIRES'],
      entry_points={'console_scripts': ['brainplot=brainplot.plotting:main']})
