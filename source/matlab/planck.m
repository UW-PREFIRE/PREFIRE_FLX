function  radiance=planck(freq,temp)
% freq in wavenumber, temp in Kelvin degree and 
% radiance in 1e-3 W per square meter per sr per wavenumber

%      ca    = 3.741832e-05 / 3.14159;
      ca = 1.191043934e-05;

%      cb    = 1.438837896;
      cb = 1.438769911;

      cof   = ca * (freq.^3);
      arg   = cb * freq;
      zeroind =	find(temp == 0);
      radiance    = cof ./ ( exp( (arg./ temp)) - 1.0 );
