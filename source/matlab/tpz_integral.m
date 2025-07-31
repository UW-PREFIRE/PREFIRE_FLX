function f = tpz_integral(h, x)
if nargin ~= 2
	disp('wrong number of input');
	return;
elseif length(x) ~= length(h)	
	disp('wrong size of x or h');
	return;
elseif diff(h) <= 0
	disp('step can only be mono increasing');
	return;
end

len = length(x);

x = reshape(x, len, 1);
h = reshape(h, len, 1);

f = sum(0.5*( (x(1:len-1) + x(2:len)).*diff(h)));
