function [flux_est_mean flux_est_unc] = rad2flux_PREFIRE_clrsky(lat, lon, ...
         TCWV, Ts, lapse, ice, snow, rad_clrsky, rad_clrsky_unc, ADM_info, ...
         wv1, wv2, sfctype, sensor_num, merged_surface_type)

% This code is to derive spectral flux from PREFIRE radiance, for clear-sky
% only and for PREFIRE channels 6-63. Flux at channels 1-5 are set to NaN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat    latitude in degree, range from -90 to 90
% lon    longitude in degree, range from 0 to 360
% TCWV   total column water vapor in cm
% Ts     surface skin temperature in K
% lapse  Ts minus air temperature at 300hPa above surface
% ice    ice fraction, unit is 1
% snow   snow fraction, unit is 1
% rad_clrsky   clear-sky radiance at PREFIRE 63 channels for each record, 
%              dimension is (channel, record), unit is (Wm^-2/sr/um)
% rad_clrsky_unc the uncertainty of rad_clrsky, unit is Wm^-2/sr/um
% wv1,wv2, lower and upper limits of each PREFIRE channel
% sfctype     IGBP surface type
% ADM_info  structure array containing ADM info
% idw_valid  valid PREFIRE channels for flux retrieval 
% sensor_num       1: TIRS1 or 2: TIRS2
% merged_surface_type   merged surface/supersurface type from AUX-MET or AUX-SAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flux_est_mean  derived clear-sky spectral flux (Wm^-2/um), same dimension as rad_clrsky
% flux_est_unc   uncertainty of derived clear-sky spectral flux (Wm^-2/um), same dimension as rad_clrsky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin < 15  | nargin >16 % count input number
	disp('wrong number of input variables');
	return;
end

using_aux = true;  % default mode

for i=1:length(lat)
        
        sid(i) = get_prefire_clr_scene_id(TCWV(i), lapse(i), Ts(i));

end


flux_est_mean = NaN(size(rad_clrsky));
flux_est_unc = NaN(size(rad_clrsky));

idw2 = 6:63;   % we only derive flux for PREFIRE channels 6-63
% idx_16 = find(idw_valid==16);  % see if channel 16 is valid
% idx_20 = find(idw_valid==20);  % see if channel 20 is valid

% Now we add channel 64 that covers 50-183 cm-1
wvdiff = wv2 - wv1;
wvdiff(64) = 1;

% {"sea ice", "partial sea ice", "open water",
%  ["permanent land ice" "Antarctic ice shelf" "snow covered land"],
%  "partial snow covered land", "snow free land"}
iset_to_aux_index = {2, 3, 1, [4 5 6], 7, 8};

 for iset =1:6

    if using_aux
       this_index = iset_to_aux_index{iset};
       if length(this_index) == 3
          type_idx = find(merged_surface_type == this_index(1) | ...
                          merged_surface_type == this_index(2) | ...
                          merged_surface_type == this_index(3));
       else
          type_idx = find(merged_surface_type == iset_to_aux_index{iset});
       end
    else

       if iset ==1  % Sea-ice
         type_idx=find(ice>=0.95 & snow<0.001);
       elseif iset==2  % Melted-ice
         type_idx=find(ice>=0.05 & ice<0.95  & snow<0.001);
       elseif iset==3  % Ocean
         type_idx=find(sfctype==17 & ice<0.05  & snow<0.001);
       elseif iset==4  % Permanent-snow
         type_idx=find(snow>=0.5);
       elseif iset==5  % Fresh-snow
         type_idx=find(snow>=0.001 & snow<0.5);
       else  % Non-snow-land
         type_idx=find( ice<0.05 & snow<0.001 & sfctype~=17);
       end
    end

    subtype = ADM_info{iset}.subtype;
    sid_ansR = ADM_info{iset}.sid_ansR;
    pca = ADM_info{iset}.pca;
    ch17_b = ADM_info{iset}.ch17_b;
    ch18_b = ADM_info{iset}.ch18_b;

   for ii=1:length(type_idx)
   
    idx=find(subtype==sid(type_idx(ii)));
    
    if isempty(idx)
        % borrow the "neighboring" scene type
       for is = 1: length(subtype)
         sidstr = num2str(sid(type_idx(ii)));
         sid1 = str2num(sidstr(1:1));
         sid2 = str2num(sidstr(2:2));
         if length(sidstr) > 2
            sid3 = str2num(sidstr(3:3));
         else
            sid3 = 0;
         end
         sidstr = num2str(subtype(is));
         sid10 = str2num(sidstr(1:1));
         sid20 = str2num(sidstr(2:2));
         if length(sidstr) > 2
            sid30 = str2num(sidstr(3:3));
         else
            sid30 = 0;
         end
         mindiff(is) = 0.5*(sid1-sid10)^2+0.2*(sid2-sid20)^2+0.3*(sid3-sid30)^2;
       end
     
       [val idx] = min(mindiff);
    end
    
    if ~isempty(idx)   
       scene_id=idx;
       if scene_id > size(sid_ansR, 1)
          continue  % scene_id out of bounds (why???), should quality flag this
       end
       
       idw_valid = find(~isnan(rad_clrsky(6:end,type_idx(ii))))+5;

      if length(idw_valid)>10  % Avoid all channels that have NaN values

       for ic=1:2
           if  ic==2
               rad_clrsky(:,type_idx(ii)) = rad_clrsky(:,type_idx(ii)) + rad_clrsky_unc(:,type_idx(ii)); 
           end

           flux_est(idw2,1)=pi.*rad_clrsky(idw2,type_idx(ii))./sid_ansR(scene_id,1:end-1)';
        
           % To get full spectral flux  %%%%%%%%%%%%%%%%
           rad2 = flux_est(idw_valid,1);
           new_mrad = pca{scene_id, 1}(idw_valid-5);   % -5 here because ADM stored are from channels 6-63
       
           new_pc = pca{scene_id, 2}(idw_valid-5, :);
 
           coeff = inv(new_pc' * new_pc) * new_pc' * (rad2' - new_mrad)';
           coeff = reshape(coeff, 1, length(coeff));
   
           % Project back to get the full spectral OLR
           %Now we add channel 64: 50-183 cm-1
           flux_est(6:64,1) =  pca{scene_id, 1}' + pca{scene_id, 2} * coeff';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
           % Flux at two CO2 channels (17 and 18) are derived using regression coefficients
           % obtained by fitting radiances at neighbour channels
            
           if sensor_num == 1
               if ~isnan(rad_clrsky(20,type_idx(ii)))
                  % don't do fitting
               else
                  flux_est(17,1) = ch17_b(1) + ch17_b(2) .* rad_clrsky(19,type_idx(ii))' + ch17_b(3) .* rad_clrsky(20,type_idx(ii))';
                  flux_est(18,1) = ch18_b(1) + ch18_b(2) .* rad_clrsky(19,type_idx(ii))' + ch18_b(3) .* rad_clrsky(20,type_idx(ii))';
               end
           elseif sensor_num == 2
               if ~isnan(rad_clrsky(16,type_idx(ii)))
                  %  don't do fitting
               else
                  flux_est(17,1) = ch17_b(1) + ch17_b(2) .* rad_clrsky(16,type_idx(ii))' + ch17_b(3) .* rad_clrsky(19,type_idx(ii))';
                  flux_est(18,1) = ch18_b(1) + ch18_b(2) .* rad_clrsky(16,type_idx(ii))' + ch18_b(3) .* rad_clrsky(19,type_idx(ii))';
               end
           end
           if ic==1
               flux_est_mean(1:64, type_idx(ii))= flux_est./wvdiff';
           elseif ic==2
               flux_est_unc(1:64, type_idx(ii))  = flux_est./wvdiff' - flux_est_mean(:, type_idx(ii));
           end

%            idn = find(flux_est_mean(:,type_idx(ii))<0 | flux_est_mean(:,type_idx(ii))>100);
%            if ~isempty(idn)
%               'sss'
%            end
       end
      else

          flux_est_mean(1:64, type_idx(ii)) = NaN;
          flux_est_unc(1:64, type_idx(ii)) = NaN;
      end
    else
       
       disp('no such scene found'); 
        
    end
   
   
   end
   
 end

