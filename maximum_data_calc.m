function[maxH,minH] = maximum_data_calc(mu,L,H0)

%% set global variables to be used in other functions
global Ap a dt dx Tcl g f mode n D V H
%% calculate maximum data

% maximum pressure difference
maxH = max(max(H())) ;
output = sprintf('maximum pressure height at valve: %0.2f m', maxH) ;
disp(output)
minH = min(min(H())) ;
output = sprintf('minimum pressure height at valve: %0.2f m', minH) ;
disp(output)
excessH = maxH-H0 ;
output = sprintf('excess pressure height at valve: %0.2f m', excessH) ;
disp(output)

% comparison with maximum pressure equation
if mu/Tcl >= 1
    deltaH = a*V(n+1,1)/g ;
else
    deltaH = L/Tcl*V(n+1,1)/g ;
end
output = sprintf('excess pressure height at valve from equations in 4.2.1: %0.2f m'...
    , deltaH) ;
disp(output)
end