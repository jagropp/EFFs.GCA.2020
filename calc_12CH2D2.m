clear

a2eq_H2O_CH3SCoM  = 1.1498;
a2eq_H2O_HSCoB    = 2.284;
a13eq_CO2_CH3SCoM = 1.069;
a2kin_Mcr_s       = 0.85;
a2kin_Mcr_p       = 0.41;
a13kin_Mcr        = 0.96;
dDH2O             = 0;
d13CCO2           = -80;
R2H2O             = (dDH2O/1000+1)*1.5576e-4;
R13CO2            = (d13CCO2/1000+1)*0.0112;
gamma_22_s        = 1;
gamma_22_p        = 0.992;
gamma_132_s       = 0.998;
gamma_132_p       = 0.992;
a22kin_Mcr_s      = a2kin_Mcr_s*a2kin_Mcr_s*gamma_22_s;
a22kin_Mcr_p      = a2kin_Mcr_p*a2kin_Mcr_s*gamma_22_p;
a132kin_Mcr_s     = a13kin_Mcr*a2kin_Mcr_s*gamma_132_s;
a132kin_Mcr_p     = a13kin_Mcr*a2kin_Mcr_p*gamma_132_p;
Delta_22eq_CH3    = 15.6063; % at 25oC.
Delta_132eq_CH3   = 5.4911; % at 25oC.

a_CH4_H2O = 0.75*(a2kin_Mcr_s/a2eq_H2O_CH3SCoM) + ...
            0.25*(a2kin_Mcr_p/a2eq_H2O_HSCoB);
klna_CH4_H2O = 1000*log(a_CH4_H2O);
a_CO2_CH4 = a13eq_CO2_CH3SCoM*(1/a13kin_Mcr);

R13CH3       = R13CO2/a13eq_CO2_CH3SCoM;
R2CH3        = R2H2O/a2eq_H2O_CH3SCoM;
R13CH4        = R13CO2/a_CO2_CH4;
R2CH4         = a_CH4_H2O*R2H2O;

R_1 = a2kin_Mcr_p/a2eq_H2O_HSCoB;

D13CH3D  = 1000.*((3*gamma_132_s*a2kin_Mcr_s*...
    (Delta_132eq_CH3/1000+1) + gamma_132_p*R_1)/...
    (3*a2kin_Mcr_s + R_1) - 1);

D12CH2D2 = 1000*((8*a2kin_Mcr_s*...
    (gamma_22_s*a2kin_Mcr_s*(Delta_22eq_CH3/1000+1)...
    + gamma_22_p*R_1))/...
    (((3*a2kin_Mcr_s + R_1)^2)) - 1);

fprintf('13CH3D = %g\n\n',D13CH3D)
fprintf('12CH2D2 = %g\n\n',D12CH2D2)





