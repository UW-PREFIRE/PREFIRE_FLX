function [sid, pft, sid6] = get_prefire_cld_scene_id(pw, fcld, sst, sctd, psi)

% use central-temperature to decide \psi
% if input \psi is out of range, use the closest one at this stage
if nargin ~= 5
	disp('wrong number of inputs');
	return;
end

if isnan(psi)
    sid = 0;
    pft = 0;
    sid6 = 0;
    return;
end

 
% cloud fraction
fmin = [0.001 0.5 0.75 0.99] ;
fmax = [0.5 0.75 0.99 1.0+1e-6];

  
 
% precipitable water (cm) 
wmin = [0 0.5 1 2 3 5];
wmax = [0.5 1 2 3 5 300];


% surface-cloud tyemperature difference 
dTmin = [-300 -15:5:85];
dTmax = [-15:5:85 300];

 % Ts (K)
Tsmin = [100 230:10:270 275:5:320];
Tsmax = [230:10:270 275:5:320 500];


ipw = find( wmin <= pw & wmax > pw);
ifcld = find( fmin <= fcld & fmax > fcld);
ist = find( Tsmin <= sst      & Tsmax > sst);
isctd = find( dTmin <= sctd & dTmax > sctd);

fr = [0.25 0.625 0.874 1];
dt = [-15 -12.5:5:82.5 85];
ts = [230 235:10:275 277.5:5:317.5 320 ];
es = 0.982;
sigma = 5.6696e-8;
sigma1 = sigma/pi;

f = fr(ifcld);
sctd = dt(isctd);
sfcT = ts(ist);
cldT = sfcT - sctd;
BTs = sigma1 .* sfcT .^ 4;
BTc = sigma1 .* cldT .^ 4;
ec = [1 0.05];
%[es  BTs  BTc]
psimm = (1-f) * es * BTs + f * (es * BTs .* (1 - ec) + ec * BTc);
if isctd <= 4
    psi_range = 4;
    if psi >= psimm(1)
       ipsi = 1;
    elseif psi <= psimm(2)
       ipsi = 4;
    else
       dpsi = (psimm(2) - psimm(1)) / psi_range;
       ipsi = floor((psi - psimm(1)) / dpsi) + 1;
    end
else
    icld = ifcld * (isctd - 4);
    if icld <= 20
	psi_range = 4;
    elseif icld <= 40
	psi_range = 12;
    else
	psi_range = 20;
    end

    if psi >= psimm(2)
       ipsi = psi_range;
    elseif psi <= psimm(1)
       ipsi = 1;
    else
       dpsi = (psimm(2) - psimm(1)) / psi_range;
       ipsi = floor((psi - psimm(1)) / dpsi) + 1;
    end
end
if ipsi > psi_range ipsi = psi_range; end
if ipsi < 1 ipsi = 1; end
 
sidstr = [int2str(ipw), int2str(ifcld), myint2str(ist,2), myint2str(isctd,2),int2str(ipsi)];
sid = str2num(sidstr);
pftstr = sidstr(1:4);
pft = str2num(pftstr);
sid6 = str2num(sidstr(1:6));

