
pfset FileVersion 4

pfset Process.Topology.P        [lindex $argv 0]
pfset Process.Topology.Q        [lindex $argv 1]
pfset Process.Topology.R        [lindex $argv 2]

pfset FlowVR [if {[lindex $argv 3] == "--FlowVR"} {list True} {list False}]

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X                -10.0
pfset ComputationalGrid.Lower.Y                 10.0
pfset ComputationalGrid.Lower.Z                  1.0

pfset ComputationalGrid.DX	                 8.8888888888888893
pfset ComputationalGrid.DY                      10.666666666666666
pfset ComputationalGrid.DZ	                 1.0


#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names "domain_input background_input source_region_input \
		       concen_region_input"


#---------------------------------------------------------
# Domain Geometry Input
#---------------------------------------------------------
pfset GeomInput.domain_input.InputType            Box
pfset GeomInput.domain_input.GeomName             domain

#---------------------------------------------------------
# Domain Geometry
#---------------------------------------------------------
pfset Geom.domain.Lower.X                        -10.0
pfset Geom.domain.Lower.Y                         10.0
pfset Geom.domain.Lower.Z                          1.0

pfset Geom.domain.Upper.X                        150.0
pfset Geom.domain.Upper.Y                        170.0
pfset Geom.domain.Upper.Z                          9.0

pfset Geom.domain.Patches "left right front back bottom top"

#---------------------------------------------------------
# Background Geometry Input
#---------------------------------------------------------
pfset GeomInput.background_input.InputType         Box
pfset GeomInput.background_input.GeomName          background

#---------------------------------------------------------
# Background Geometry
#---------------------------------------------------------
pfset Geom.background.Lower.X -99999999.0
pfset Geom.background.Lower.Y -99999999.0
pfset Geom.background.Lower.Z -99999999.0

pfset Geom.background.Upper.X  99999999.0
pfset Geom.background.Upper.Y  99999999.0
pfset Geom.background.Upper.Z  99999999.0


#---------------------------------------------------------
# Source_Region Geometry Input
#---------------------------------------------------------
pfset GeomInput.source_region_input.InputType      Box
pfset GeomInput.source_region_input.GeomName       source_region

#---------------------------------------------------------
# Source_Region Geometry
#---------------------------------------------------------
pfset Geom.source_region.Lower.X    65.56
pfset Geom.source_region.Lower.Y    79.34
pfset Geom.source_region.Lower.Z     4.5

pfset Geom.source_region.Upper.X    74.44
pfset Geom.source_region.Upper.Y    89.99
pfset Geom.source_region.Upper.Z     5.5


#---------------------------------------------------------
# Concen_Region Geometry Input
#---------------------------------------------------------
pfset GeomInput.concen_region_input.InputType       Box
pfset GeomInput.concen_region_input.GeomName        concen_region

#---------------------------------------------------------
# Concen_Region Geometry
#---------------------------------------------------------
pfset Geom.concen_region.Lower.X   60.0
pfset Geom.concen_region.Lower.Y   80.0
pfset Geom.concen_region.Lower.Z    4.0

pfset Geom.concen_region.Upper.X   80.0
pfset Geom.concen_region.Upper.Y  100.0
pfset Geom.concen_region.Upper.Z    6.0

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names "background"

pfset Geom.background.Perm.Type     Constant
pfset Geom.background.Perm.Value    4.0

pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "background"

pfset Geom.background.Perm.TensorValX  1.0
pfset Geom.background.Perm.TensorValY  1.0
pfset Geom.background.Perm.TensorValZ  1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain"
pfset Geom.domain.SpecificStorage.Value 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	Constant
pfset Phase.water.Density.Value	1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names			""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------
pfset Geom.Retardation.GeomNames           ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

pfset TimingInfo.BaseUnit		0.0001
pfset TimingInfo.StartCount		0
pfset TimingInfo.StartTime		$::env(START_TIME)
pfset TimingInfo.StopTime               $::env(STOP_TIME)
pfset TimeStep.Value                    0.0001
pfset TimingInfo.DumpInterval	       -1
pfset TimeStep.Type                     Constant

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          background

pfset Geom.background.Porosity.Type    Constant
pfset Geom.background.Porosity.Value   1.0

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          domain
pfset Geom.domain.RelPerm.Alpha        0.005
pfset Geom.domain.RelPerm.N            2.0

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type            VanGenuchten
pfset Phase.Saturation.GeomNames       domain
pfset Geom.domain.Saturation.Alpha     0.005
pfset Geom.domain.Saturation.N         2.0
pfset Geom.domain.Saturation.SRes      0.2
pfset Geom.domain.Saturation.SSat      0.99

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names constant
pfset Cycle.constant.Names		"alltime"
pfset Cycle.constant.alltime.Length	 1
pfset Cycle.constant.Repeat		-1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames "left right front back bottom top"

pfset Patch.left.BCPressure.Type			DirEquilRefPatch
pfset Patch.left.BCPressure.Cycle			"constant"
pfset Patch.left.BCPressure.RefGeom			domain
pfset Patch.left.BCPressure.RefPatch			bottom
pfset Patch.left.BCPressure.alltime.Value		5.0

pfset Patch.right.BCPressure.Type			DirEquilRefPatch
pfset Patch.right.BCPressure.Cycle			"constant"
pfset Patch.right.BCPressure.RefGeom			domain
pfset Patch.right.BCPressure.RefPatch			bottom
pfset Patch.right.BCPressure.alltime.Value		3.0

pfset Patch.front.BCPressure.Type			FluxConst
pfset Patch.front.BCPressure.Cycle			"constant"
pfset Patch.front.BCPressure.alltime.Value		0.0

pfset Patch.back.BCPressure.Type			FluxConst
pfset Patch.back.BCPressure.Cycle			"constant"
pfset Patch.back.BCPressure.alltime.Value		0.0

pfset Patch.bottom.BCPressure.Type			FluxConst
pfset Patch.bottom.BCPressure.Cycle			"constant"
pfset Patch.bottom.BCPressure.alltime.Value		0.0

pfset Patch.top.BCPressure.Type			        FluxConst
pfset Patch.top.BCPressure.Cycle			"constant"
pfset Patch.top.BCPressure.alltime.Value		0.0

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames ""

pfset TopoSlopesX.Geom.domain.Value 0.0

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames ""

pfset TopoSlopesY.Geom.domain.Value 0.0

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames ""
pfset Mannings.Geom.domain.Value 0.

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      3.0
pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   bottom

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    background
pfset PhaseSources.water.Geom.background.Value        0.0

#----------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
#Less output to compare:
#Writing output (all pfb):
pfset Solver.PrintSubsurfData False
pfset Solver.PrintPressure False
pfset Solver.PrintSaturation False
pfset Solver.PrintMask False
pfset Solver.WriteCLMBinary False
pfset Solver.PrintCLM False
pfset Solver.WriteSiloSpecificStorage False
pfset Solver.WriteSiloMannings False
pfset Solver.WriteSiloMask False
pfset Solver.WriteSiloSlopes False
pfset Solver.WriteSiloSubsurfData False
pfset Solver.WriteSiloPressure False
pfset Solver.WriteSiloSaturation False
pfset Solver.WriteSiloEvapTrans False
pfset Solver.WriteSiloEvapTransSum False
pfset Solver.WriteSiloOverlandSum False
pfset Solver.WriteSiloCLM False


#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution


#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver                                             Richards
#pfset Solver.MaxIter                                     5

pfset Solver.Nonlinear.MaxIter                           10
pfset Solver.Nonlinear.ResidualTol                       1e-9
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          1e-5
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-2

pfset Solver.Linear.KrylovDimension                      10

pfset Solver.Linear.Preconditioner                       PFMGOctree
pfset Solver.Linear.Preconditioner.MGSemi.MaxIter        1
pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels      100

pfset NetCDF.NumStepsPerFile			1
