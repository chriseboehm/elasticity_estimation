%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sample for the simple version of the model
% Date: July 2014
% Authors: Nitya and Chris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

global mu alpha epsi om psi phi_i nu_i w_tm1 R_tm1 pnJ_tm1_i pJ_tm1_i x_tm1_i ...
    int_agg prod_fn

% Options

rng(12322)   % set seed for random number generator


% Parameters

mu = 0.5;
alpha = 1/3;
om = 0.2;
psi = 0.1;
epsi = 10;

phi_i = 1;          % later this may vary by firm
nu_i = 0.5;         % later this may vary by firm

% Functions

int_agg = @(mnJ,mJ) (nu_i.*(mnJ).^((om-1)/om) + ...
    (1-nu_i).*(mJ).^((om-1)/om)).^(om/(om-1));

prod_fn = @(K,L,mnJ,mJ) phi_i.*(mu*(K.^alpha.*L.^(1-alpha)).^((psi-1)/psi) + ...
    (1-mu).*int_agg(mnJ,mJ).^((psi-1)/psi) ).^(psi/(psi-1));

% prices are fixed

w_tm1 = 1;
R_tm1 = 1;


% First order conditions in period t-1
% Notice that we have constant returns to scale and perfect competition so
% the size of the firm will not be pinned down. For that reason, I will fix
% capital as an input.

N = 200;        % sample size

x_tm1     = 1+rand(N,1);
pnJ_tm1   = 1+rand(N,1);      
pJ_tm1    = 1+rand(N,1);  

% px_tm1    = (mu^psi*((R_tm1./alpha).^alpha.*(w_tm1./(1-alpha)).^(1-alpha)).^(1-psi) + ...
%     (1-mu)^psi.*((nu_i.^om.*pnJ_tm1.^(1-om)+(1-nu_i).^om.*pJ_tm1.^(1-om) ).^(1/(1-om))) ...
%     .^(1-psi) ).^(1/(1-psi));

K_tm1     = zeros(N,1);
L_tm1     = zeros(N,1);
mnJ_tm1   = zeros(N,1);
mJ_tm1    = zeros(N,1);
px_tm1    = zeros(N,1);  

startval = [1,1,1,1,1];


tic
for i=1:N
    
    x_tm1_i = x_tm1(i);
    pnJ_tm1_i = pnJ_tm1(i);
    pJ_tm1_i = pJ_tm1(i);
    sample = fsolve('foc_obj',startval);
    K_tm1(i) = sample(1);
    L_tm1(i) = sample(2);
    mnJ_tm1(i) = sample(3);
    mJ_tm1(i) = sample(4);
    px_tm1(i) = sample(5);

end
toc


% Our proxy with north american exports

kappa_i = 0.1 + 0.1.*rand(N,1);

VNAX_tm1    = kappa_i.*px_tm1.*x_tm1;
VmJ_tm1     = pJ_tm1.*mJ_tm1;
VmnJ_tm1    = pnJ_tm1.*mnJ_tm1;

% The following constitutes the sample for t-1:

% K_tm1
% L_tm1
% VNAX_tm1
% VmJ_tm1
% VmnJ_tm1
% pnJ_tm1
% pJ_tm1
% px_tm1

% Now set up sample for period t

K_t     = K_tm1;
L_t     = L_tm1;
px_t    = px_tm1;
pJ_t    = pJ_tm1;
pnJ_t   = pnJ_tm1;
mnJ_t   = mnJ_tm1;
mJ_t    = (0.5 + 0.5.*rand(N,1)).*mJ_tm1;
x_t     = prod_fn(K_t,L_t,mnJ_t,mJ_t);

VNAX_t  = kappa_i.*px_t.*x_t + 0.20.*randn(N,1);
VmJ_t   = pJ_t.*mJ_t;
VmnJ_t  = pnJ_t.*mnJ_t;

% The following constitutes the sample for t:

% K_t
% L_t
% VNAX_t
% VmJ_t
% VmnJ_t
% pnJ_t
% pJ_t
% px_t

save('sample.mat','K_tm1','L_tm1','VNAX_tm1','VmJ_tm1','VmnJ_tm1',...
    'pnJ_tm1','pJ_tm1','px_tm1', 'K_t', 'L_t', 'VNAX_t', 'VmJ_t', ...
    'VmnJ_t', 'pnJ_t', 'pJ_t', 'px_t');

%%

mu  = 0.5;
om  = 0.5;
psi = 0.5;


nu_i = (VmnJ_tm1./VmJ_tm1).^(1/om)./((pJ_tm1./pnJ_tm1).^((om-1)/om) + ...
    (VmnJ_tm1./VmJ_tm1).^(1/om));

temp1 = (mu*(px_tm1.*K_tm1.^alpha.*L_tm1.^(1-alpha)).^((psi-1)/psi) + ...
    (1-mu).*(px_tm1./pJ_tm1.*VmJ_tm1).^((psi-1)/psi).* ...
    ((nu_i.*(pJ_tm1./pnJ_tm1.*VmnJ_tm1./VmJ_tm1).^((om-1)/om) + ...
    (1-nu_i)).^(om/(om-1))).^((psi-1)/psi)).^(psi/(psi-1));

phi_i_tild = VNAX_tm1./temp1;

temp2 = (mu*(px_t.*K_t.^alpha.*L_t.^(1-alpha)).^((psi-1)/psi) + ...
    (1-mu).*(px_t./pJ_t.*VmJ_t).^((psi-1)/psi).* ...
    ((nu_i.*(pJ_t./pnJ_t.*VmnJ_t./VmJ_t).^((om-1)/om) + ...
    (1-nu_i)).^(om/(om-1))).^((psi-1)/psi)).^(psi/(psi-1));

eps_NA = (VNAX_t - phi_i_tild.*temp2);

f = -eps_NA'*eps_NA;

% Test whether we get the estimation right when nu is measured with error


 

