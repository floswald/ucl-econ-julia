%% ------------------------------------------------------------------------ 
% Dynamic Economics in Practice. 
% Monica Costa Dias and Cormac O'Dea
% Institute for Fiscal Studies
% 28-29th October 2013


%% ------------------------------------------------------------------------ 
% DESCRIPTION
% This program solves and simulates a finite period consumption and saving 
% problem. There is income that can be uncertain. The income program can be
% hardcoded by the user or can be set to follow a log normal autoregessive
% process

%% ------------------------------------------------------------------------ 
% PREAMBLE
% Ensure that all storage spaces variables, matrices, memory, globals are 
% clean from information stored in past runs

tic;        % start the clock
clear all;  % clear memory
close all;  % close any graphs


%% ------------------------------------------------------------------------ 
% DECLARE VARIABLES AND MATRICES THAT WILL BE 'GLOBAL'
% explicitly set variables and matrices to be shared throughout the routine
% as globals

global beta gamma r T sigma mu rho Tretire          % structural model parameters (1)
global isUncertainty                                % structural model parameters (2)
global numPointsA Agrid                             % assets grid and dimension
global numPointsY Ygrid incTransitionMrx            % income grid, dimension, and transition matrices (1)
global hcIncome hcIncPDF                            % income grid, dimension, and transition matrices (2)
global interpMethod linearise uncertaintyMethod     % numerical methods to be used
global tol minCons numSims normBnd                  % numerical constants
global plotNumber                                   % for graphing


%% ------------------------------------------------------------------------ 
% NUMERICAL METHODS
% select solution, interpolation and integration methods

interpMethod = 'linear';      % interpolation method - choose from 'linear', 'nearest', 'spline', 'pchip'
linearise = 1;               % whether to linearise the slope of EV when using EE - set linearise=1 to do this, else = 0


%% ------------------------------------------------------------------------ 
% NUMERICAL CONSTANTS
% set constants needed in numerical solution and simulation

% precision parameters
%--------------------------------------%
tol = 1e-10;                 % max allowed error
minCons = 1e-5;              % min allowed consumption

% where to truncate the normal distributions
%--------------------------------------%
normBnd = 3;                 %Ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma 

% information for simulations
%--------------------------------------%
numSims = 2;                % How many individuals to simulate


%% ------------------------------------------------------------------------ 
% THE ECONOMIC ENVIRONMENT
% Set values of structural economic parameters

T = 40;                      % Number of time period
r = 0.01;                   % Interest rate
beta = 0.98;                  % Discount factor
gamma = 1.5;                 % Coefficient of relative risk aversion
mu = 0;                      % mean of initial log income
sigma = 0.25;                   % variance of log income 
rho = 0.75;                     % persistency of log income
Tretire = 41;                % age after which there is no income earned
borrowingAllowed = 0;        % Is borrowing allowed
isUncertainty = 1;           % Is income uncertain?
startA = 0;                  % How much asset do people start life with

%% ------------------------------------------------------------------------ 
% GRIDS
% choose dimensions, set matrices and select methods to construct grids

%The grid for assets
%--------------------------------------%
numPointsA = 20;             % number of points in the discretised asset grid
gridMethod = '3logsteps';    % method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps

%The grid for income shocks
%--------------------------------------%
numPointsY = 5;           %  points in grid for income (should be 2 if hard-coded)
uncertaintyMethod =1;    %  =0 if we enter shocks manually, =1 if put shocks on a grid using Tauchen (1986) method
hcIncome = [0.5941    0.8179    1.0000    1.2227    1.6832]';   % hard-coded shocks, used if uncertaintyMethod = 1
hcIncPDF = [0.1016    0.2492    0.2983    0.2492    0.1016]';  % and respective probabilities

%% Check inputs
checkInputs;

%% Get income grid
[Ygrid, incTransitionMrx, minInc, maxInc] = getIncomeGrid;

%% ------------------------------------------------------------------------ 
% GET ASSET GRID
% populate grid for assets using 'gridMethod'
% populate matrix borrowingConstraints


[ borrowCon, maxAss ] = getMinAndMaxAss(borrowingAllowed, minInc, maxInc, startA);

Agrid = NaN(T+1, numPointsA);
for ixt = 1:1:T+1
    Agrid(ixt, :) = getGrid(borrowCon(ixt), maxAss(ixt), numPointsA, gridMethod);
end


%% ------------------------------------------------------------------------ 
% SOLVE CONSUMER'S PROBLEM
% Get policy function and value function 

[ policyA1, policyC, val, exVal ] = solveValueFunction;

% Export it in order to compare with Julia
dlmwrite('../v5_julia/matlablObj/policyA1matlab.csv',policyA1, ',');


%% ------------------------------------------------------------------------ 
% SIMULATE CONSUMER'S PATHS
% start from initial level of assets and simulate optimal consumption and
% savings profiles over lifecycle

if isUncertainty == 0
    [ ypath, cpath, apath, vpath ] = simNoUncer(policyA1, exVal, startA);
else 
    [ ypath, cpath, apath, vpath, e,logy1 ] = simWithUncer(policyA1,exVal, startA);
    
    dlmwrite('../v5_julia/matlablObj/ematlab.csv',e, ',');  % Export them so we can add them into Julia
    dlmwrite('../v5_julia/matlablObj/logy1matlab.csv',logy1, ',');
end

toc;     % Stop the clock
%% ------------------------------------------------------------------------ 
% PLOTS
% Plots some features of the simulation and simulation

% Plot paths of consumption, income and assets 
plotNumber = 0;
%plotCpath(cpath)
plotApath(apath, borrowCon)
plotYAndCpaths( ypath, cpath );
%plotYCAndApaths( ypath, cpath, apath );

% Now plot value and policy functions
whichyear = 20;
plotNode1 = 3;
plotNodeLast = numPointsA; 
plots; 


% ------------------------------------------------------------------------ 
% ------------------------------------------------------------------------   