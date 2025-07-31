function pw = get_total_water_content(pres, T, q, Ps, Ts)

 
% this is to get total water content in cm
R = 8.314472/18.06e-3;

newp= [Ps; pres]; 
newT = [Ts; T];

q = [q(1); q];

% q is in g/kg, convert back to kg/kg
% p is in mb, convert to SI unit, Pasacal
pw = tpz_integral(newp(end:-1:1)*1e2, q(end:-1:1)*1e-3);
% -1 becase of the bottome up integrating (dp < 0)
pw =  pw/ 9.806; 

% convert kg/m^2 to m as if condensed, then to cm, for 1atm, 273K
rho_h2o = 1e3;

pw = 100 * pw / rho_h2o;

