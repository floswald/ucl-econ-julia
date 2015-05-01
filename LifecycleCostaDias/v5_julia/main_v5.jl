
# Paul Rodriguez translation into Julia 0.3.5 of Monica Costa Dias and Corma O'Dea "Dynamic Economics in Practice"
# using Alex code for Getting Equiprobable Shock Points from a Normal Distributions
# March 2015, p.lesmes.11@ucl.ac.uk
# Please improve it!

# ------------------------------------------------------------------------
# DESCRIPTION
# This program solves and simulates a finite period consumption and saving
# problem. There is income that can be uncertain. The income program can be
# hardcoded by the user or can be set to follow a log normal autoregessive
# process

# IMPORTANT: Comment lines 188 and 189 if you want to change the number of 
# periods or use shocks generated in Julia

# ------------------------------------------------------------------------
# PREAMBLE
# Ensure that all storage spaces variables, matrices, memory, globals are
# clean from information stored in past runs


using Distributions   # For drwaing distributions
using Optim           # Basic optimization routines
using Grid            # Interpolation tools


# Set here the relevant directory. Please check 
cd("C://Users//Paul//Documents//GitHub//ucl-econ-julia//LifecycleCostaDias//v5_julia")

# ----------------------------------------------------------------------

mainDir=pwd() 

# Numeric tools
cd("code_numericTools")
require("getEquiprobNormalDeviates.jl")
require("truncate.jl")
require("interp1lin.jl")
require("interp2D.jl") 
cd(mainDir)

# Grids and other inputs
cd("code_gridsInputs")
require("checkInputs.jl")
require("getIncomeGrid.jl")
require("getMinandMaxAss.jl")
require("getGrid.jl")
cd(mainDir)

# Programs for the solution
cd("code_solution")
require("solveValueFunction.jl")
require("objectivefunc.jl")
require("utility.jl")
cd(mainDir)

# Programs for the simulation
cd("code_simulation")
require("simNoUncer.jl")
require("simWithUncer.jl")
cd(mainDir)

# Plot stuff
cd("code_plots")
require("plotPaths.jl")
require("plots.jl")
cd(mainDir)

tic()        # start the clock

# ------------------------------------------------------------------------
# DECLARE VARIABLES AND MATRICES THAT WILL BE 'GLOBAL'
# explicitly set variables and matrices to be shared throughout the routine
# as globals

#% ------------------------------------------------------------------------
# NUMERICAL METHODS
# select solution, interpolation and integration methods

#global interpMethod = 'pchip';      # interpolation method - IN JULIA, UP TO THE MOMENT, YOU CAN ONLY USE LINEAR WITH UNEVEN GRIDS
const linearise = 1;               # whether to linearise the slope of EV when using EE - set linearise=1 to do this, else = 0


# ------------------------------------------------------------------------
# NUMERICAL CONSTANTS
# set constants needed in numerical solution and simulation

# precision parameters
#--------------------------------------%
const tol = 1e-10;                 # max allowed error
const minCons = 1e-5;              # min allowed consumption

# where to truncate the normal distributions
#--------------------------------------%
const normBnd = 3;                 #Ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma


# information for simulations
#--------------------------------------%
const numSims = 2;                #How many individuals to simulate


# ------------------------------------------------------------------------
# THE ECONOMIC ENVIRONMENT
# Set values of structural economic parameters

const T = 40;                      # Number of time period (If you want to reproduce matlab results, be sure that it has the same number of periods!)
const r = 0.01;                    # Interest rate
const beta = 0.98;                 # Discount factor
const gamma = 1.5;                 # Coefficient of relative risk aversion
const mu = 0;                      # mean of initial log income
const sigma = 0.25;                   # variance of log income
const rho = 0.75;                     # persistency of log income
const Tretire = 41;                # age after which there is no income earned
const borrowingAllowed = 0;        # Is borrowing allowed
const isUncertainty = 1;           # Is income uncertain?
const startA = 0;                  # How much asset do people start life with

# ------------------------------------------------------------------------
# GRIDS
# choose dimensions, set matrices and select methods to construct grids

#The grid for assets
#--------------------------------------%
const numPointsA = 20;             # number of points in the discretised asset grid
gridMethod = "3logsteps";    # method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps

#The grid for income shocks
#--------------------------------------%
const numPointsY = 5;           #  points in grid for income (should be 2 if hard-coded)
const uncertaintyMethod = 1;    #  =0 if we enter shocks manually, =1 if put shocks on a grid using Tauchen (1986) method
hcIncome = [0.5941 ,   0.8179  ,  1.0000  ,  1.2227  ,  1.6832];  # hard-coded shocks, used if uncertaintyMethod = 1
hcIncPDF = [0.1016 ,   0.2492  ,  0.2983  ,  0.2492  ,  0.1017];  # and respective probabilities

# Check inputs
checkInputs()


## Get income grid
(Ygrid, incTransitionMrx, minInc, maxInc) = getIncomeGrid();


## ------------------------------------------------------------------------
# GET ASSET GRID
# populate grid for assets using 'gridMethod'
# populate matrix borrowingConstraints

(borrowCon, maxAss) = getMinAndMaxAss(borrowingAllowed, minInc, maxInc, startA);

Agrid = Array(Float64,T+1, numPointsA);
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod);
end

#% ------------------------------------------------------------------------
# SOLVE CONSUMER'S PROBLEM
# Get policy function and value function

(policyA1, policyC, val, exVal) = solveValueFunction();

	# Want to compare with Matlab version? In that way you can test how different are
	# simulations for the same policyFunc in Matlab and Julia
	# policyA1=reshape(readdlm("matlablObj\\policyA1matlab.csv", ',' ),40,20,5);

#% ------------------------------------------------------------------------
# SIMULATE CONSUMER'S PATHS
# start from initial level of assets and simulate optimal consumption and
# savings profiles over lifecycle

if isUncertainty == 0
    (ypath, cpath, apath, vpath) = simNoUncer(policyA1, exVal, startA);
else
	# Get the draws, and the feed the simulation with them (sligthly different than the original code)

	# For the discreate shocks of Income case... Julia makes our life way easier!
	srand(1223424) # Seed it!
	randYpath=rand( Categorical(vec(hcIncPDF))  ,(T,numSims));

	# For the normal distributed case
	srand(1223424) # Seed it!
 	e=rand(Normal(0, sigma),(T,numSims)); # normally distributed random draws for the innovation
 	srand(234636)  # Seed it!
 	sig_inc = sigma/ ((1-rho^2)^0.5);
 	logy1   =rand(Normal(mu, sig_inc),numSims);  # a random draw for the initial income

		# IMPORTANT: Want to compare with Matlab version? In that way you can test how different are
		# simulations using the same shocks. NOTICE THAT THIS OVERIDES THE GENERATED SHOCKS FROM
		# THE PREVIOUS LINES
		# Comment this lines if you want to change the number of periods
 		e      = readdlm("matlablObj\\ematlab.csv", ',' );
 		logy1  = readdlm("matlablObj\\logy1matlab.csv", ',' );

    (ypath, cpath, apath, vpath) = simWithUncer(policyA1,exVal, startA,e,logy1,randYpath);
end


toc();     # Stop the clock

#% ------------------------------------------------------------------------
# PLOTS
# Plots some features of the simulation and simulation
# Use a program that supports Gadfly, in these days, either JUNO, IJULIA

if 1==1  # In purspouse, Gadfly is just so slow!! But there are few good alternatives (March 2015)

	# Slow.........
	using Gadfly
	using DataFrames      # For making graphs or exporting data, is better as dataframes


	plotPaths()

	# Now plot value and policy functions
	whichyear = 20;
	plotNode1 = 3;
	plotNodeLast = numPointsA;

	plots();   #Plot Value Functions
end

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
