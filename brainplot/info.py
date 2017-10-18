__version__ = "0.0.1"

__url__ = "https://github.com/mgxd/brainplot"
__packagename__ = "brainplot"
__description__ = "Pretty brain plots!"
__longdesc__ = ("Plot various statistics on the Conte69 atlas generated from "
                "69 healthy adults registered to the `fs_LR` surface mesh, and "
                "the HCP500 resting atlas. By default, plots are rendered off-"
                "screen using xvfbwrapper. Multiple views can be generated with"
                " the `--view` flag, and varying levels of brain mesh inflation"
                " can be generated with the `--inflation` flag.")

MAINTAINER = 'Mathias Goncalves'
MAINTAINER_EMAIL = 'mathiasg@mit.edu'


PROVIDES = ["brainplot"]

# currently, problems with latest VTK - TODO: set version limits
REQUIRES = [
    "mayavi",
    "matplotlib",
    "nibabel>=2.2.0",
    "xvfbwrapper",
]

# TODO: add testing
TESTS_REQUIRES = ['pytest-cov']
