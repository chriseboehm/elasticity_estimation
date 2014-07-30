%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates elasticities
% Date: July 2014
% Authors: Nitya and Chris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

global K_tm1 L_tm1 VNAX_tm1 VmJ_tm1 VmnJ_tm1 pnJ_tm1 pJ_tm1 px_tm1 ...
    K_t L_t VNAX_t VmJ_t VmnJ_t pnJ_t pJ_t px_t

load('sample.mat')

% run estimation

% mu = 0.5;
% om = 0.2;
% psi = 0.1;

%mu0  = 0.5;
om0  = 0.3;
psi0 = 0.3;

startval = [om0 psi0]; %[mu0 om0 psi0];

theta = fmincon('est_obj',startval,[],[],[],[],[0 0 0],[1 10 10]);

disp(theta)