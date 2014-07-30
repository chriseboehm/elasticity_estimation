function f = est_obj(args)

global K_tm1 L_tm1 VNAX_tm1 VmJ_tm1 VmnJ_tm1 pnJ_tm1 pJ_tm1 px_tm1 ...
    K_t L_t VNAX_t VmJ_t VmnJ_t pnJ_t pJ_t px_t

% mu  = args(1);
mu  = 0.5; %+0.2
om  = args(1);
psi = args(2);

alpha = 1/3;

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

eps_NA = VNAX_t - phi_i_tild.*temp2;

f = eps_NA'*eps_NA;

end