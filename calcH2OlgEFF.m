function kln2a_H2O_lg = calcH2OlgEFF(Tc)

% Calculate equilibrium water liquid-vapor hydrogen isotope fractionation
% based on experimental results of Horita and Wesolowski (GCA, 1994), for
% the range 0-374 oc. The that there are two available fits, the fit in the
% website AlphaDelta is slightly different, with up to 1.35 difference at
% 0oc.
%
% Input 
%             Tc - Temperature in degree celcius
%
% Output
%             kln2a_H2O_lg - EFF for given temperature

Tk = Tc + 273.15;

% Fit from the original paper
% kln2a_H2O_lg = 1158.8.*(Tk.^3./1e9) - 1620.1.*(Tk.^2./1e6) + 794.84.*(Tk./1e3) - ...
%     161.04 + 2.9992.*(1e9./Tk.^3);

% Fit from website AlphaDelta
kln2a_H2O_lg = -0.353e18./Tk.^6 + 42.17e12./Tk.^4 - 309.4e9./Tk.^3 + ...
    963.7e6./Tk.^2 - 1399e3./Tk + 766.2;
    
if Tc > 374
    kln2a_H2O_lg = NaN;
end

end