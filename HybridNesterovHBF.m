%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNessterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Nonstrongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

clear all

set(0,'defaultTextInterpreter','tex');

% global variables
global delta M gamma lambda c_0 c_10 r tauMin tauMax tauMed c zeta cTilde_0 cTilde_10 d_0 d_10 alpha

%%%%%%%%%%%%%%%%%%%%%%%
% setting the variables
%%%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
M = 2;
zeta = 1;

% Heavy Ball constants
gamma = 2/3; % 
lambda = 40; % 

% Uniting parameters for \mathcal{U}_0 and \mathcal{T}_{1,0}:
c_0 = 320; % \mathcal{U}_0 
c_10 = 271.5835157; % \mathcal{T}_{1,0} 

% These are the same, since L = x^2
alpha = 1;

% eps_0 has to be bigger than eps_10
eps_0 = 10;
eps_10 = 5; % 5 10 12 15

cTilde_0 = eps_0*alpha
cTilde_10 = eps_10*alpha
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)

tauMin = 3;

c = 0.25;
deltaMed = 50000; %4010
r = 51; 

delta = 0.5;

deltaVecUniting = [0,0,0];
deltaVec = [0,0,0];
deltaVecHBF = [0,0,0];
deltaVecUniting = [0,0,0];

lDeltaNest = 0;
lDeltaHBF = 0;
lDeltaUniting = 0;
lDelta = 0;

% initial conditions
z1_0 = 50; 
z2_0 = 0;
z2_00 = 50;
q_0 = 1;
tau_0 = 0;
tauPN_0 = tauMin; 

tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN_HBF=[0 700];
TSPAN=[0 20];
JSPAN = [0 20000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

%% simulate

% Nesterov
[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest = zeros(1,length(tNest));
[deltaVecNest lNest lDeltaNest] = timeToConv(xNest,tNest);

lNestAvg = lNest.';
% This is the dotted line indicating the average value:
PO = polyfit(tNest,log10(lNestAvg(:,1)),1);
yO = polyval(PO,tNest);

% heavy ball
[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF = zeros(1,length(tHBF));
[deltaVecHBF lHBF lDeltaHBF] = timeToConv(xHBF,tHBF);

deltaVecHBF = timeToConv(xHBF,tHBF);

% Uniting
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lUniting = zeros(1,length(tUniting));
[deltaVecUniting lUniting lDeltaUniting] = timeToConv(xUniting,tUniting);

% HAND-1
[t,j,x] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lHAND = zeros(1,length(t));
[deltaVec lHAND lDelta] = timeToConv(x,t);

lHANDAvg = lHAND.';
% This is the dotted line indicating the average value:
PH = polyfit(t,log10(lHANDAvg(:,1)),1);
yH = polyval(PH,t);

minarc = min([length(x),length(xUniting)]);
ta = [tUniting(1:minarc),t(1:minarc)];
ja = [jUniting(1:minarc),j(1:minarc)];
xa = [xUniting(1:minarc,1),x(1:minarc,1)];
xb = [xUniting(1:minarc,2),x(1:minarc,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra simulations, for averaging percent improvement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulating the algorithms from z1(0,0) = 20:
% Retune r and deltaMed for HAND-1
deltaMed = 8090; 
r = 21;

c_0 = 200;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 74.9533626; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.2;

% initial conditions
z1_0 = 20; 
z2_00 = 20;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU1,jU1,xU1] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU1 = zeros(1,length(tU1));
[deltaVecU1 lU1 lDeltaU1] = timeToConv(xU1,tU1);

% HAND-1
[tH1,jH1,xH1] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH1= zeros(1,length(tH1));
[deltaVecH1 lH1 lDeltaH1] = timeToConv(xH1,tH1);

% Heavy ball
[tHBF1,jHBF1,xHBF1] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF1 = zeros(1,length(tHBF1));
[deltaVecHBF1 lHBF1 lDeltaHBF1] = timeToConv(xHBF1,tHBF1);

% Nesterov
[tNest1,jNest1,xNest1] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest1 = zeros(1,length(tNest1));
[deltaVecNest1 lNest1 lDeltaNest1] = timeToConv(xNest1,tNest1);

%% Simulating the algorithms from z1(0,0) = 30:
% Retune r and deltaMed for HAND-1
deltaMed = 18090; 
r = 31;

c_0 = 200;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 121.770066;
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.3;

% initial conditions
z1_0 = 30; 
z2_00 = 30;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU2,jU2,xU2] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU2 = zeros(1,length(tU2));
[deltaVecU2 lU2 lDeltaU2] = timeToConv(xU2,tU2);

% HAND-1
[tH2,jH2,xH2] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH2= zeros(1,length(tH2));
[deltaVecH2 lH2 lDeltaH2] = timeToConv(xH2,tH2);

% Heavy ball
[tHBF2,jHBF2,xHBF2] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF2 = zeros(1,length(tHBF2));
[deltaVecHBF2 lHBF2 lDeltaHBF2] = timeToConv(xHBF2,tHBF2);

% Nesterov
[tNest2,jNest2,xNest2] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest2 = zeros(1,length(tNest2));
[deltaVecNest2 lNest2 lDeltaNest2] = timeToConv(xNest2,tNest2);

%% Simulating the algorithms from z1(0,0) = 40:
% Retune r and deltaMed for HAND-1
deltaMed = 32060; 
r = 41;

c_0 = 300;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 187.3134505; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.4;

% initial conditions
z1_0 = 40; 
z2_00 = 40;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU3,jU3,xU3] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU3 = zeros(1,length(tU3));
[deltaVecU3 lU3 lDeltaU3] = timeToConv(xU3,tU3);

% HAND-1
[tH3,jH3,xH3] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH3 = zeros(1,length(tH3));
[deltaVecH3 lH3 lDeltaH3] = timeToConv(xH3,tH3);

% Heavy ball
[tHBF3,jHBF3,xHBF3] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF3 = zeros(1,length(tHBF3));
[deltaVecHBF3 lHBF3 lDeltaHBF3] = timeToConv(xHBF3,tHBF3);

% Nesterov
[tNest3,jNest3,xNest3] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest3 = zeros(1,length(tNest3));
[deltaVecNest3 lNest3 lDeltaNest3] = timeToConv(xNest3,tNest3);

%% Simulating the algorithms from z1(0,0) = 60:
% Retune r and deltaMed for HAND-1
deltaMed = 71910; 
r = 61;

c_0 = 500;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 374.580265; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.6;

% initial conditions
z1_0 = 60; 
z2_00 = 60;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU4,jU4,xU4] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU4 = zeros(1,length(tU4));
[deltaVecU4 lU4 lDeltaU4] = timeToConv(xU4,tU4);

% HAND-1
[tH4,jH4,xH4] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH4 = zeros(1,length(tH4));
[deltaVecH4 lH4 lDeltaH4] = timeToConv(xH4,tH4);

% Heavy ball
[tHBF4,jHBF4,xHBF4] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF4 = zeros(1,length(tHBF4));
[deltaVecHBF4 lHBF4 lDeltaHBF4] = timeToConv(xHBF4,tHBF4);

% Nesterov
[tNest4,jNest4,xNest4] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest4 = zeros(1,length(tNest4));
[deltaVecNest4 lNest4 lDeltaNest4] = timeToConv(xNest4,tNest4);

%% Simulating the algorithms from z1(0,0) = 70:
% Retune r and deltaMed for HAND-1
deltaMed = 97800; 
r = 71;

c_0 = 650;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 496.303695; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.7;

% initial conditions
z1_0 = 70; 
z2_00 = 70;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU5,jU5,xU5] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU5 = zeros(1,length(tU5));
[deltaVecU5 lU5 lDeltaU5] = timeToConv(xU5,tU5);

% HAND-1
[tH5,jH5,xH5] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH5 = zeros(1,length(tH5));
[deltaVecH5 lH5 lDeltaH5] = timeToConv(xH5,tH5);

% Heavy ball
[tHBF5,jHBF5,xHBF5] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF5 = zeros(1,length(tHBF5));
[deltaVecHBF5 lHBF5 lDeltaHBF5] = timeToConv(xHBF5,tHBF5);

% Nesterov
[tNest5,jNest5,xNest5] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest5 = zeros(1,length(tNest5));
[deltaVecNest5 lNest5 lDeltaNest5] = timeToConv(xNest5,tNest5);

%% Simulating the algorithms from z1(0,0) = 80:
% Retune r and deltaMed for HAND-1
deltaMed = 127650; 
r = 81;

c_0 = 800;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 636.753805; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.8;

% initial conditions
z1_0 = 80; 
z2_00 = 80;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU6,jU6,xU6] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU6 = zeros(1,length(tU6));
[deltaVecU6 lU6 lDeltaU6] = timeToConv(xU6,tU6);

% HAND-1
[tH6,jH6,xH6] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH6 = zeros(1,length(tH6));
[deltaVecH6 lH6 lDeltaH6] = timeToConv(xH6,tH6);

% Heavy ball
[tHBF6,jHBF6,xHBF6] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF6 = zeros(1,length(tHBF6));
[deltaVecHBF6 lHBF6 lDeltaHBF6] = timeToConv(xHBF6,tHBF6);

% Nesterov
[tNest6,jNest6,xNest6] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest6 = zeros(1,length(tNest6));
[deltaVecNest6 lNest6 lDeltaNest6] = timeToConv(xNest6,tNest6);

%% Simulating the algorithms from z1(0,0) = 90:
% Retune r and deltaMed for HAND-1
deltaMed = 161500;
r = 91;

c_0 = 1000;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 795.930591; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 0.9;

% initial conditions
z1_0 = 90; 
z2_00 = 90;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU7,jU7,xU7] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU7 = zeros(1,length(tU7));
[deltaVecU7 lU7 lDeltaU7] = timeToConv(xU7,tU7);

% HAND-1
[tH7,jH7,xH7] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH7 = zeros(1,length(tH7));
[deltaVecH7 lH7 lDeltaH7] = timeToConv(xH7,tH7);

% Heavy ball
[tHBF7,jHBF7,xHBF7] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF7 = zeros(1,length(tHBF7));
[deltaVecHBF7 lHBF7 lDeltaHBF7] = timeToConv(xHBF7,tHBF7);

% Nesterov
[tNest7,jNest7,xNest7] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest7 = zeros(1,length(tNest7));
[deltaVecNest7 lNest7 lDeltaNest7] = timeToConv(xNest7,tNest7);

%% Simulating the algorithms from z1(0,0) = 100:
% Retune r and deltaMed for HAND-1
deltaMed = 199300; 
r = 101;

c_0 = 1200;  
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 973.834065; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 1;

% initial conditions
z1_0 = 100; 
z2_00 = 100;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU8,jU8,xU8] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU8 = zeros(1,length(tU8));
[deltaVecU8 lU8 lDeltaU8] = timeToConv(xU8,tU8);

% HAND-1
[tH8,jH8,xH8] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH8 = zeros(1,length(tH8));
[deltaVecH8 lH8 lDeltaH8] = timeToConv(xH8,tH8);

% Heavy ball
[tHBF8,jHBF8,xHBF8] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF8 = zeros(1,length(tHBF8));
[deltaVecHBF8 lHBF8 lDeltaHBF8] = timeToConv(xHBF8,tHBF8);

% Nesterov
[tNest8,jNest8,xNest8] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest8 = zeros(1,length(tNest8));
[deltaVecNest8 lNest8 lDeltaNest8] = timeToConv(xNest8,tNest8);

%% Simulating the algorithms from z1(0,0) = 110:
% Retune r and deltaMed for HAND-1
deltaMed = 241050; 
r = 111;

c_0 = 1400; 
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)

% Retune c_10, d_10, and delta:
c_10 = 1170.46423; 
d_10 = c_10 - (cTilde_10/alpha)^2 - (zeta^2/(M))*((cTilde_10^2)/alpha)
delta = 1.1;

% initial conditions
z1_0 = 110; 
z2_00 = 110;

% Setting the rest of the HAND-1 parameters (since z1_0 now set):
tauMed = sqrt(((r^2)/(2*c) + (tauMin^2)*CalculateL(z1_0))/deltaMed) + tauMin
tauMax = tauMed + 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0;tau_0];
x00 = [z1_0;z2_00;tauPN_0];
x000 = [z1_0;z2_0;q_0];

% Uniting:
[tU9,jU9,xU9] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Uniting:
lU9 = zeros(1,length(tU9));
[deltaVecU9 lU9 lDeltaU9] = timeToConv(xU9,tU9);

% HAND-1
[tH9,jH9,xH9] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for HAND-1:
lH9 = zeros(1,length(tH9));
[deltaVecH9 lH9 lDeltaH9] = timeToConv(xH9,tH9);

% Heavy ball
[tHBF9,jHBF9,xHBF9] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x000,TSPAN_HBF,JSPAN,rule,options,'ode45');

% Find the L values for HBF:
lHBF9 = zeros(1,length(tHBF9));
[deltaVecHBF9 lHBF9 lDeltaHBF9] = timeToConv(xHBF9,tHBF9);

% Nesterov
[tNest9,jNest9,xNest9] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

% Find the L values for Nesterov:
lNest9 = zeros(1,length(tNest9));
[deltaVecNest9 lNest9 lDeltaNest9] = timeToConv(xNest9,tNest9);

%% Plots
figure(1) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
subplot(2,1,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(1),'k.','MarkerSize', 20)
strDelta = [num2str(deltaVec(3)), 's'];
text(deltaVec(3),deltaVec(1),strDelta,'HorizontalAlignment','left','VerticalAlignment','top');
plot(deltaVecUniting(3),deltaVecUniting(1),'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),deltaVecUniting(1),strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([0 20 -50 80])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
subplot(2,1,2), plotHarc(ta,ja,xb,[],modificatorF,modificatorJ);
hold on
plot(deltaVec(3),deltaVec(2),'k.','MarkerSize', 20)
plot(deltaVecUniting(3),deltaVecUniting(2),'k.','MarkerSize', 20)
axis([0 20 -100 70])
grid on
ylabel('z_2','Fontsize',16)
xlabel('t','Fontsize',16)
hold off
saveas(gcf,'Plots\ComparisonPlotsNSC','png')

figure(2)
clf
semilogy(tUniting,lUniting,'LineWidth',3);
hold on
semilogy(tHBF,lHBF,'Color',[0.4660 0.6740 0.1880],'LineWidth',3);
semilogy(tNest,lNest,'Color',[0.6350 0.0780 0.1840],'LineWidth',3,'LineStyle','--');
semilogy(tNest,10.^(yO(:,1)),'k--','LineWidth',3);
semilogy(t,lHAND,'Color',[0.8500 0.3250 0.0980],'LineWidth',3);
semilogy(t,10.^(yH(:,1)),'k:','LineWidth',3);

% Plotting times to convergence:
semilogy(deltaVec(3),lDelta,'k.','MarkerSize', 20)
strDelta = [num2str(deltaVec(3)), 's'];
text(deltaVec(3),lDelta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
semilogy(deltaVecUniting(3),lDeltaUniting,'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),lDeltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
semilogy(deltaVecNest(3),lDeltaNest,'k.','MarkerSize', 20)
strDeltaNest = [num2str(deltaVecNest(3)), 's'];
text(deltaVecNest(3),lDeltaNest,strDeltaNest,'HorizontalAlignment','left','VerticalAlignment','bottom');

hold off
axis([0 20 10^(-28) 10^(6)]);
ylabel('L(z_1)-L^*','FontSize',20)
xlabel('t','FontSize',20)
legend({'Hybrid','Heavy ball','Nesterov','Nesterov, average','HAND-1','HAND-1, average'},'Location','southwest')

axes('Position',[0.7 0.6 0.15 0.1])
box on
hold on
semilogy(tHBF,lHBF,'Color',[0.4660 0.6740 0.1880],'LineWidth',3);
semilogy(deltaVecHBF(3),lDeltaHBF,'k.','MarkerSize', 20)
strDeltaHBF = [num2str(deltaVecHBF(3)), 's'];
text(deltaVecHBF(3),lDeltaHBF,strDeltaHBF,'HorizontalAlignment','right','VerticalAlignment','bottom');
hold off
set(gca,'xtick',[0 100 200])
set(gca,'ytick',[10^(-2) 10^(2) 10^(6)])
axis([0 200 10^(-2) 10^(6)])
hold off

saveas(gcf,'Plots\Semilog','epsc')
saveas(gcf,'Plots\Semilog','png')

figure(3)
clf
semilogy(tU9,lU9,'LineWidth',3)
axis([0 5 10^(-20) 10^5])
grid on
ylabel('L(z_1)-L^*','Fontsize',16)
xlabel('t','Fontsize',16)
saveas(gcf,'Plots\TrajectoryPlotsTest','png')
