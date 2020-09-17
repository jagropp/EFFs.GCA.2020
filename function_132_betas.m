% 11.02.2020
%
% A function to caluclate beta values for EFFs of H, C and clumping from
% the coeficients of the polynomila fits for the beta values. Make sure
% that the .csv file (Table_main2_beta_fit.csv) is in the same folder as
% the function.

function beta_values = function_132_betas(Temp)

% INPUT
%        Temp           - Temperature in degree Celsius
%
% OUTPUT
%        compound_names - a list of the names of the metabolites
%        beta2_values   - a list of temperature dependent beta for H
%        beta13_values  - a list of temperature dependent beta for C
%        beta132_values - a list of temperature dependent beta for 13C-D clumps
%        beta22_values  - a list of temperature dependent beta for D-D clumps

file_name  = 'beta_fit_values.csv';
beta_coefs   = readtable(file_name);

Tc = Temp;
Tk = Tc + 273.15;
beta2_values =   beta_coefs.A_2b(2:19).*beta_coefs.A_2b(1).*1e12./Tk.^4 + ...
                 beta_coefs.B_2b(2:19).*beta_coefs.B_2b(1).*1e9./Tk.^3 + ...
                 beta_coefs.C_2b(2:19).*beta_coefs.C_2b(1).*1e6./Tk.^2 + ...
                 beta_coefs.D_2b(2:19).*beta_coefs.D_2b(1).*1e3./Tk + ...
                 beta_coefs.E_2b(2:19).*beta_coefs.E_2b(1);

beta13_values =  beta_coefs.A_13b(2:19).*beta_coefs.A_13b(1).*1e12./(Tk.^4) + ...
                 beta_coefs.B_13b(2:19).*beta_coefs.B_13b(1).*1e9./(Tk.^3) + ...
                 beta_coefs.C_13b(2:19).*beta_coefs.C_13b(1).*1e6./Tk.^2 + ...
                 beta_coefs.D_13b(2:19).*beta_coefs.D_13b(1).*1e3./Tk + ...
                 beta_coefs.E_13b(2:19).*beta_coefs.E_13b(1);
           
beta132_values = beta_coefs.A_132b(2:19).*beta_coefs.A_132b(1).*1e12./(Tk.^4) + ...
                 beta_coefs.B_132b(2:19).*beta_coefs.B_132b(1).*1e9./(Tk.^3) + ...
                 beta_coefs.C_132b(2:19).*beta_coefs.C_132b(1).*1e6./Tk.^2 + ...
                 beta_coefs.D_132b(2:19).*beta_coefs.D_132b(1).*1e3./Tk + ...
                 beta_coefs.E_132b(2:19).*beta_coefs.E_132b(1);

beta22_values  = beta_coefs.A_22b(2:19).*beta_coefs.A_22b(1).*1e12./(Tk.^4) + ...
                 beta_coefs.B_22b(2:19).*beta_coefs.B_22b(1).*1e9./(Tk.^3) + ...
                 beta_coefs.C_22b(2:19).*beta_coefs.C_22b(1).*1e6./Tk.^2 + ...
                 beta_coefs.D_22b(2:19).*beta_coefs.D_22b(1).*1e3./Tk + ...
                 beta_coefs.E_22b(2:19).*beta_coefs.E_22b(1);
             
compound_names = beta_coefs.Compound(2:end);

beta_values    = table(compound_names,beta2_values,...
                       beta13_values,beta132_values,beta22_values);
end




