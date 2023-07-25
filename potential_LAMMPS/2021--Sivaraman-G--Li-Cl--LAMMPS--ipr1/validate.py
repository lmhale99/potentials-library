"""
@author: gsivaraman@anl.gov
"""
import numpy as np
from ase.io import read, write
from sklearn.metrics import mean_absolute_error


inp = read('test.xyz',':')
print(len(inp))
inenergy = [ei.get_potential_energy() for ei  in inp ]
output = read('quip_train.xyz',':')
print(len(output))
outenergy = [eo.get_potential_energy() for eo in output ]
print("\nE MAE: {}".format( mean_absolute_error(np.asarray(outenergy), np.asarray(inenergy)) ) )



in_force, out_force = [], []


for num in range(len(inp)):
    incar =  inp[num]
    outcar = output[num]
    chemind = inp[num].get_chemical_symbols()
    for ind in range(len(chemind)):
        in_force.append(incar.get_forces()[ind])
        out_force.append(outcar.arrays['force'][ind])


def rms_dict(x_ref, x_pred):
    """ Takes two datasets of the same shape and returns a dictionary containing RMS error data"""

    x_ref = np.array(x_ref)
    x_pred = np.array(x_pred)

    if np.shape(x_pred) != np.shape(x_ref):
        raise ValueError('WARNING: not matching shapes in rms')

    error_2 = (x_ref - x_pred) ** 2

    average = np.sqrt(np.average(error_2))
    std_ = np.sqrt(np.var(error_2))

    return {'rmse': average, 'std': std_}

_rms = rms_dict(in_force, out_force)
print("\nF metrics",_rms)
