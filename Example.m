%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Example                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; 

% Set parameters 
n = 200;                     % image size 
savename = 'TestSol';        % savename 
epsilonStart = 0.1;         % initial epsilon 
epsilonEnd = 0.01;           % final epsilon 
maxIter = 1000;              % max. number of iterations
tol = 1e-6;                  % error tolerance 
a0 = 1e10;                   % a0
a = 0.05;                    % a
b = 1;                       % b
example = '3Points';         % Example (see getExample.m)
initSol = false;             % given initial solution
savenameInit = 'TestInit';   % filename/savename of given initial solution

%-------------------------------------------------------------------------
% Compute solution 
numPhaseFields = length(a);
if ( numPhaseFields == 1 && a0 > 1e7 ) 
    [sigma,phi,lambda] = SPFS(n,epsilonStart,epsilonEnd,maxIter,tol,a,b,example,savename);
elseif ( numPhaseFields > 1 && a0 > 1e7 )
    [sigma,phi,lambda] = MPFS(n,epsilonStart,epsilonEnd,maxIter,tol,a,b,example,savename,initSol,savenameInit);
else 
    [sigma,phi,lambda] = MPFSD(n,epsilonStart,epsilonEnd,maxIter,tol,a0,a,b,example,savename,initSol,savenameInit); 
end
%-------------------------------------------------------------------------