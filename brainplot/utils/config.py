"""
Optional configuration file for automatically setting atlas paths.
Allows for different paths based on current host.
"""

import os
import socket
host = socket.gethostname()

# Set atlas paths
CONTE69_ATLAS = ''
REST_ATLAS  = ''

if host == 'openmind7':
    base = '/om/user/mathiasg/scripts/templates'
    CONTE69_ATLAS = os.path.join(base, 'Conte69_Atlas')
    REST_ATLAS = os.path.join(base, 'rfMRI_REST1_LR_Atlas.dtseries.nii')
elif host == 'mathias-N501VW':
    base = '/home/mathias/code/datasets'
    CONTE69_ATLAS = os.path.join(base, 'Conte69_Atlas')
    REST_ATLAS = os.path.join(base, 'rfMRI_REST1_LR_Atlas.dtseries.nii')
