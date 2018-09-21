function [Qsca,L_mie] = computeMie(r,d,phi)
% Define needed parameters and call Mie scattering code
%
%   Inputs
%       r - pore radius
%       d - porous layer thickness 
%       phi - porosity of porous layer
%
%   Outputs
%       Qsca - scattering efficiency 
%       L_mie - calculated scattering losses through porous layer
%
%   To download MATLAB Mie scattering code visit: 
%       https://omlc.org/software/mie/
%   Download "Maetzler's MATLAB code for Mie theory". 
%   Place codes in current working directory or add to path.
%

    % Center frequency of ApRES
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

