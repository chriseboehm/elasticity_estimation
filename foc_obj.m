function f = foc_obj(args)
    
global mu alpha epsi om psi phi_i nu_i w_tm1 R_tm1 pnJ_tm1_i pJ_tm1_i x_tm1_i ...
    int_agg prod_fn

K   = args(1);
L   = args(2);     
mnJ = args(3);
mJ  = args(4);
px_tm1_i = args(5);
x = x_tm1_i;


f(1) = (1-1/epsi)*px_tm1_i.*phi_i^(1-1/psi)*x^(1/psi)*mu*...
    (K^alpha*L^(1-alpha))^((psi-1)/psi)*(1-alpha) - w_tm1*L;

f(2) = (1-1/epsi)*px_tm1_i.*phi_i^(1-1/psi)*x^(1/psi)*mu*...
    (K^alpha*L^(1-alpha))^((psi-1)/psi)*alpha - R_tm1*K;

M = int_agg(mnJ,mJ);

f(3) = (1-1/epsi)*px_tm1_i.*phi_i^(1-1/psi)*x^(1/psi)*(1-mu)*M^(1/om-1/psi)...
    *nu_i*mnJ^(-1/om) - pnJ_tm1_i;

f(4) = (1-1/epsi)*px_tm1_i.*phi_i^(1-1/psi)*x^(1/psi)*(1-mu)*M^(1/om-1/psi)...
    *(1-nu_i)*mJ^(-1/om) - pJ_tm1_i;

f(5) = prod_fn(K,L,mnJ,mJ) - x;

end