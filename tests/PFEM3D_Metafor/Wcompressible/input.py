import cupydo.interfaces.Cupydo as cupy
import os

# %% Input Parameters

def getFsiP():

    param = {}
    path = os.path.abspath(os.path.dirname(__file__))
    param['cfdFile'] = path+'\\'

    # Metafor and PFEM solvers
    
    param['fluidSolver'] = 'Pfem3D'
    param['solidSolver'] = 'Metafor'
    param['cfdFile'] += 'input_pfem'
    param['csdFile'] = 'input_meta'
    
    # FSI objects

    param['interpolator'] = 'Matching'
    param['criterion'] = 'Displacements'
    param['algorithm'] = 'AitkenBGS'
    
    # FSI parameters

    param['compType'] = 'unsteady'
    param['computation'] = 'direct'
    param['timeItTresh'] = 0
    param['nDim'] = 2
    param['dt'] = 1e-6
    param['tTot'] = 0.3
    param['tol'] = 1e-6
    param['maxIt'] = 25
    param['omega'] = 0.5
    
    return param

param = getFsiP()
cupydo = cupy.CUPyDO(param)
cupydo.run()