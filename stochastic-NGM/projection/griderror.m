function ssr = griderror(coeff,param, qnodes, qweights, kgrid, zgrid)

alfa =param(1);
beta =param(2);
gamma=param(3);
delta=param(4);
sigma=param(5);
rho  =param(6);


nk =  size(kgrid,1);
nz =  size(zgrid,1);
nq =  size(qnodes,1);

ssr      =  0;

for ik = 1:nk
    for iz = 1:nz
        k = kgrid(ik);
        z = zgrid(iz);
        c = consfun(k,z,coeff);
    knext = z*k^alfa+(1-delta)*k-c;

        expec  = 0;
            for iq = 1:nq
                znext = exp(rho*log(z)+sqrt(2)*sigma*qnodes(iq));
                cnext = consfun(knext,znext,coeff);

                expec = expec + qweights(iq)*beta* cnext^(-gamma)*(alfa*znext*knext^(alfa-1)+1-delta);            
            end
        expec = expec/sqrt(pi);  
        ssr = ssr+(expec-c^(-gamma))^2;
    end
end


