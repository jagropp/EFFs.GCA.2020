

Tc = linspace(0,500,501);
Tk = Tc + 273.15;

EFF_CO2_CH4_1000lna = 0.2054.*1e12./Tk.^4 - 2.3655.*1e9./Tk.^3 + ...
                      11.754.*1e6./Tk.^2 + 0.160;


plot(Tc,EFF_CO2_CH4_1000lna,Tc,EFF_CO2_CH4_calc)
xlabel('Tc')
ylabel('CO_2-CH_4')




