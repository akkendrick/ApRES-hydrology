function [Qsca,L_mie] = computeMie(r,d,phi)
% Compute mie scattering

    f = 300*10^6;

    % Looking up index of refractions at center ApRES wavelength (~100 cm) 
    nwater = 8.847727 + 1.0i * 6.727E-02;
    nice = 1.7861 + 1.0i*3.348E-004;
    nvac = 1; 

    c = 2.998 * 10^8; %m/s
    pi = 3.14159;
    lambda = c/f;

    dia = 2*r;

    musgp = getMieScatter(lambda,dia,phi,nwater,nice);
    Qsca = musgp(4);

    tau_s = (3*phi*d)/(4*r) * Qsca;
    L_mie_term = exp(-2*tau_s);
    L_mie = 10*log10(L_mie_term);

end

