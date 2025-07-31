function str = myint2str(data, strlen)

tmp = 10^(strlen-1);

str = [];

if data > 0
flag = 0;
for i = 1:strlen-1
	if data < 10^i & data >= 10^(i-1)
		for k = 1:strlen-i
		str = [str, '0'];
		end
		str = [str, int2str(data)];
		flag = 1;
	end
end

if flag == 0 
	str = int2str(data);
end

elseif data == 0
	str = repmat('0', 1, strlen);

else
	disp('input must be non-negative integer');
	str =' ';
	return
end
