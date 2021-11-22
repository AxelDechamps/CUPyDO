Problem = {
    id = "WCompNewtonNoT",
	simulationTime = math.huge,
	verboseOutput = false,
	useCupydo = true,
	
	Mesh = {
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		hchar = 0.01,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.01,-0.01,0.61,100},
		ignoreGroups = {"SolidBase","Solid"},
		mshFile = "geometry.msh",
		exclusionZones = {}
	},
	
	Extractors = {
		{
			kind = "GMSH",
			writeAs = "NodesElements",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p","velocity"},
			outputFile = "fluid.msh"
		}
	},
	
	Material = {
		mu = 1e-3,
		rhoStar = 1000,
        K0 = 2.2e+7,
        K0p = 7.6,
		gamma = 0
	},
	
	IC = {
		ReservoirFixed = true,
		FSInterfaceFixed = false
	},
	
	Solver = {
		adaptDT = true,
	    id = "CDS_dpdt",
		securityCoeff = 0.1,
		initialDT = math.huge,
		maxDT = math.huge,

        MomEq = {
            bodyForce = {0,-9.81},
            BC = {}
        },
        
        ContEq = {
            stabilization = "Meduri",
            BC = {}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0,0,0,Problem.Material.rhoStar,0,0}
end

function Problem.Solver.MomEq.BC:ReservoirV(pos,t)
	return {0,0}
end

function Problem.Solver.MomEq.BC:FSInterfaceV(pos,t)
	return nil
end