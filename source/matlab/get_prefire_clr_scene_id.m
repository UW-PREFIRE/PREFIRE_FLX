function sid = get_prefire_clr_scene_id(pw, vtc, st)

% This code is to define scene subtype based on pw, vtc and st, over polar
% region only

if nargin < 3  | nargin >4 
	disp('wrong number of input variables');
	return;
end

 
% precipitable water (cm) 0-1, 1-3, 3-5, > 5
wmin = [0 0.5 1 2 3 5];
wmax = [0.5 1 2 3 5 300];

% deltaT (K)
dTmin = [-900 -20 -10 0 10 20 30 40];
dTmax = [-20 -10 0 10 20 30 40 900];

 % Ts (K)
Tsmin = [100, 230, 250, 270, 290, 310, 330];
Tsmax = [230, 250, 270, 290, 310, 330, 500];


if nargin ~= 3
	disp('wrong number of inputs');
	return;
end

id1 = find( wmin < pw & wmax > pw);
id2 = find( dTmin <= vtc & dTmax > vtc);
id3 = find( Tsmin <= st      & Tsmax > st);

sidstr = [int2str(id1), int2str(id2), int2str(id3)];
sid = str2num(sidstr);
