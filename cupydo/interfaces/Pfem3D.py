from ..genericSolvers import FluidSolver
from slpp import slpp as lua
import pfem3Dw as w
import numpy as np
import os

# %% Translate the Lua Table Into a Dict

def read(path):

    if(os.path.isfile(path)):
        with open(path,'r') as file:

            input = file.read().replace(' ','')
            input = input.replace('Problem=','')
            input = lua.decode(input)

    else: raise Exception('Cannot open lua file '+path)
    return input

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        print('\n***************************** Initializing PFEM3D *****************************')
        path = param['cfdFile']+'.lua'
        input = read(path)

        # Problem class initialization

        if input['id'] == 'IncompNewtonNoT':
            self.problem = w.ProbIncompNewton(path)
            self.typeBC = 'velocity'

        elif input['id'] == 'WCompNewtonNoT':
            self.problem = w.ProbWCompNewton(path)
            self.typeBC = 'acceleration'

        else: raise Exception('Problem type not supported '+input['id'])
        if not input['useCupydo']: raise Exception('UseCupydo must be True')

        # Gets some important objects and variables

        self.solver = self.problem.getSolver()
        self.mesh = self.problem.getMesh()
        self.load = w.VectorArrayDouble3()
        self.FSID = w.VectorInt()
        self.prevMesh = w.Mesh()

        # Number of nodes at the FSInterface

        self.mesh.getNodesIndexOfTag('FSInterface',self.FSID)
        self.nPhysicalNodes = self.FSID.size()
        self.nNodes = self.FSID.size()
        self.dim = self.mesh.getDim()

        # Initializes the tracking data

        self.pos = self.getNodalInitialPositions()
        self.problem.updateTime(-param['dt'])
        self.disp = np.zeros((self.nNodes,3))
        self.BC = np.zeros((self.nNodes,3))
        self.dt = param['dt']
        self.reload = False
        self.setNodalBC()

        # Initializes the fluid solver

        FluidSolver.__init__(self)
        self.problem.displayParams()
        self.problem.dump()
        
# %% Calculates One Increment From t1 to t2

    def run(self,t1,t2):

        time = self.problem.getCurrentSimTime()
        self.solver.setTimeStep(t2-t1)
        self.solver.computeNextDT()
        progress = True

        # Solves with remeshing until t2

        while progress:
            
            dt = self.solver.getTimeStep()
            if (dt+time-t2)/dt > -0.2:
                
                dt = t2-time
                self.solver.setTimeStep(dt)
                ok = self.solver.solveOneTimeStep()
                if ok: progress = False

            else: ok = self.solver.solveOneTimeStep()
            if ok: self.remesh()

            # Prints the current state

            self.solver.computeNextDT()
            time = self.problem.getCurrentSimTime()
            print(ok,': t = {:.6f} - dt = {:.2e}'.format(time,dt))

        # Computes nodal fluid loads

        self.setCurrentState()
        self.reload = True

# %% Get and Set Nodal Values

    def getNodalInitialPositions(self):

        coord = np.zeros((self.nNodes,3))
        for i in range(self.nNodes):

            idx = self.FSID[i]
            node = self.mesh.getNode(idx)
            coord[i] = [node.getCoordinate(j) for j in range(3)]

        return np.transpose(coord)

    # Sets current node states

    def setCurrentState(self):

        self.solver.computeLoads(self.load)
        for i in range(self.nNodes):

            idx = self.FSID[i]
            self.nodalLoad_X[i] = -self.load[idx][0]
            self.nodalLoad_Y[i] = -self.load[idx][1]
            self.nodalLoad_Z[i] = -self.load[idx][2]

    # Returns the index of the index-th interface node

    def getNodalIndex(self,index):
        return self.FSID[index]

    # Prescribes nodal positions from solid solver

    def applyNodalDisplacements(self,dx,dy,dz,dx_nM1,dy_nM1,dz_nM1,haloNodesDisplacements,time):

        if self.reload:
            self.problem.loadMesh(self.prevMesh,self.prevTime)
            self.FSID = w.VectorInt(self.prevFSID)

        vel = np.zeros((self.nNodes,3))
        deltaDisp = np.transpose([dx,dy,dz])-self.disp

        # BC for incompressible fluids

        if self.typeBC == 'velocity':
            self.BC = deltaDisp/self.dt
            self.setNodalBC()

        # BC for weakly compressible fluids

        elif self.typeBC == 'acceleration':
            for i in range(self.nNodes):
                for j in range(self.dim):
                    vel[i,j] = self.mesh.getNode(self.FSID[i]).getState(j)

            self.BC = 2*(deltaDisp-vel*self.dt)/(self.dt**2)
            self.setNodalBC()

# %% Update and Save Results

    def update(self,dt):

        self.problem.copyMesh(self.prevMesh)
        self.prevTime = self.problem.getCurrentSimTime()
        self.disp = np.transpose(self.getNodalInitialPositions()-self.pos)
        self.prevFSID = w.VectorInt(self.FSID)
        FluidSolver.update(self,dt)
        self.reload = False

    def save(self,nt):
        self.problem.writeExtractors()

    def setNodalBC(self):

        nodalBC = w.MapIntArrayDouble3()
        for i in range(self.nNodes): nodalBC[self.FSID[i]] = self.BC[i]
        self.solver.setNodalBC(nodalBC)

    def remesh(self):

        self.mesh.remesh(self.problem.isOutputVerbose())
        self.mesh.getNodesIndexOfTag('FSInterface',self.FSID)
        self.setNodalBC()

# %% Exits the Solver

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\n***************************** Exit PFEM3D *****************************')