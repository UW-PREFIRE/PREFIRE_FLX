function [sctd,pseudo,ec] = get_psi_data(fcld,cloudt,st,...
  tvis,wi, re_in,ri_in,lwp,es) 

% following values are from Table 8, Minnis et al. 1998
% for the transformation from Vis Tau to IR Tau
re = [0 2 4 6 8 12 16 32 3200];  % add end*100 point (ensure re_in is in range)
qvis = [2.277 2.277 2.191 2.143 2.118 2.090 2.074 2.029 2.029];
qir = [0.335 0.335 0.749 1.156 1.499 1.973 2.208 2.271 2.271];  % 10.8 um values
w0 = [0.1556 0.1556 0.3166 0.4109 0.4659 0.5178 0.5323 0.5092 0.5092];

ri = [0 5.83 18.15 23.86 30.36 41.20 45.30 67.6 75.2 104.9 123.1 134.9 13490];  % add end*100 point (ensure ri_in is in range)
qvis2 = [2.101 2.101 2.049 2.027 2.023 2.016 2.018 2.012 2.006 2.004 2.002 2.001 2.001];
qir2 = [0.846 0.846 1.410 1.752 1.815 1.935 1.892 1.947 2.021 2.014 2.012 2.022 2.022];
w02 = [0.1777 0.1777 0.3556 0.4110 0.4315 0.4667 0.4566 0.4763 0.4997 0.5028 0.5054 0.5102 0.5102];

% Stefan-Boltzmann constant
sigma = 5.6696e-8;
sigma1 = sigma/pi;

       
        ct = cloudt;
        sctd = st - ct;

        

        tvislogt = log(tvis);
        tvislogt(tvislogt > 1e5) = 100;
        tvislog = tvislogt;
        tvist = exp(tvislogt);
        tvist(tvist > 1e5) = 0;         % to avoid NaN
        
        
       
        rwater = wi;

        ret = re_in;
        ret(ret > 1e5 & wi >0.5) = 10; % adapted from CERES ADM
        rit = ri_in;
        rit(rit > 1e5 & wi <0.5) = 50;

        qire = interp1(re, qir, re_in);
        qiri = interp1(ri, qir2, ri_in);

        qvise = interp1(re, qvis, re_in); 
        qvisi = interp1(ri, qvis2, ri_in); 

        w0w = interp1(re, w0, re_in);
        w0i = interp1(ri, w02, ri_in);

    
        
        % tau vis water for mix cloud
        tvisw = 0.75 .* qvise .* lwp ./ re_in;
        
        % tau vis ice for mix cloud
        tvisi = exp((tvislog - rwater .* log(tvisw)) ./ ...
            (1 - rwater));
       
        % considering pure and mix phase for water and ice
        indmix = find(rwater > 0.0001 & rwater < 0.9999);
        indw = find(rwater >= 0.9999);
        indi = find(rwater <= 0.0001 & rwater >= 0);
        
        % convert to absorbed tau_ir in consideration with phase, Q, and ref tau_vis       
        tir = zeros(size(rwater));

        tir(indw) = (1 - w0w(indw)) .* tvis(indw) .* qire(indw) ./ qvise(indw);
        
        tir(indi) = (1 - w0i(indi)) .* tvis(indi) .* qiri(indi) ./ qvisi(indi);
        
        tirw = (1 - w0w) .* tvisw .* qire ./ qvise;
        
        tiri = (1 - w0i) .* tvisi .* qiri ./ qvisi;
        
        indileg = find(tirw <= 0 | tirw >= 128 | ...
            tiri <= 0 | tiri >= 128);
        indileg = find(tirw <= 0 | tirw >= 128 | ...
            tiri <= 0 | tiri >= 128);
        indileg = intersect(indileg,indmix);
        
        tir(indmix) = exp(log(tirw(indmix)) .* rwater(indmix) + ...
            log(tiri(indmix)) .* (1 - rwater(indmix)));

        ec = 1 - exp(-1 .* tir);
        fclr = 1-fcld;

        % Contribution from earth
        pseudo = fclr .* es .* sigma1 .* st .^ 4;
       
        
        
        pseudo = pseudo + fcld .* (es .* (1 - ec) .* sigma1 .* st .^ 4 + ...
           ec .* sigma1 .* cloudt .^ 4);
       
        
    
        pseudo = pseudo ./ 100;
