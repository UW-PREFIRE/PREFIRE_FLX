function [flux_est_mean flux_est_unc] = rad2flux_PREFIRE_overcast(lat, lon, ...
         TCWV, Ts, ice, snow, fcld, cloudt, tvis, wi, re_in, ri_in, lwp, ...
         es, rad_overcast, rad_overcast_unc, ADM_info, wv1, wv2, sensor_num)

% This code is to derive spectral flux from PREFIRE radiance, for overcast
% only and for PREFIRE channels 6-63. Flux at channels 1-5 are set to NaN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat    latitude in degree, range from -90 to 90
% lon    longitude in degree, range from 0 to 360
% TCWV   total column water vapor in cm
% Ts     surface skin temperature in K
% lapse  Ts minus air temperature at 300hPa above surface
% ice    ice fraction, 0-1
% snow   snow fraction, 0-1
% fcld   cloud fraction, 0-1
% cloudt  cloud top temperature in K
% tvis   cloud optical depth
% wi     the ratio of water cloud path to total cloud path
% re_in  the radius of water cloud
% ri_in  the radius of ice cloud
% lwp    the cloud water path
% es     the broadband  surface emissivity
% rad_overcast   overcast radiance at PREFIRE 63 channels for each record, 
%              dimension is (channel, record), unit is (Wm^-2/sr/um)
% rad_overcast_unc  the uncertainty of rad_overcast, unit is (Wm^-2/sr/um)
% wv1,wv2, lower and upper limits of each PREFIRE channel
% ADM_info  structure array containing ADM info
% sensor_num       1: TIRS1 or 2: TIRS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flux_est_mean  derived overcast spectral flux (Wm^-2/um), same dimension as rad_overcast
% flux_est_unc   uncertainty of derived overcast spectral flux (Wm^-2/um), same dimension as rad_overcast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 20  | nargin >21 % count input number
	disp('wrong number of input variables');
	return;
end


% mask=load('land_mask_1degree.mat'); % 1 is land, other is not land
sid =[];
for i=1:length(lat)
%         [val idlat] = min(abs(lat(i) - mask.lat));
%         [val idlon] = min(abs(lon(i) - mask.lat));
%         land(i) = mask.mask(idlat,idlon);
        [sctd(i),pseudo(i),ec(i)] = get_psi_data(fcld(i),cloudt(i),Ts(i),...
                                                   tvis(i),wi(i), re_in(i),ri_in(i),lwp(i),es(i)); 
        [sid(i,:), pft(i,:), sid6(i,:)] = get_prefire_cld_scene_id(TCWV(i), fcld(i), Ts(i), sctd(i), pseudo(i));  
end


flux_est_mean=zeros(size(rad_overcast)) + NaN;
flux_est_unc=zeros(size(rad_overcast)) + NaN;   

idw2 = 6:63;   % we only derive flux for PREFIRE channels 6-63

% Now we add channel 64 that covers 0-183 cm-1
wvdiff = wv2 - wv1;
wvdiff(64) = 1;

 for iset =1 %[1:2 4:6] 
 
%    if iset ==1  % Sea-ice
%         type_idx=find(ice>=0.95 & snow<0.001);  
%    elseif iset==2  % Melted-ice
%         type_idx=find(ice>=0.05 & ice<0.95  & snow<0.001); 
%    elseif iset==3  % Ocean
%         type_idx=find(sfctype==17 & ice<0.05  & snow<0.001); 
%    elseif iset==4  % Permanent-snow
%         type_idx=find(snow>=0.5);
%    elseif iset==5  % Fresh-snow
%         type_idx=find(snow>=0.001 & snow<0.5);
%    else  % Non-snow-land
%         type_idx=find( ice<0.05 & snow<0.001 & sfctype~=17); 
%    end

    type_idx = 1:length(sid);

    subtype = ADM_info{iset}.subtype;
    sid_ansR = ADM_info{iset}.sid_ansR;
    pca = ADM_info{iset}.pca;
    ch17_b = ADM_info{iset}.ch17_b;
    ch18_b = ADM_info{iset}.ch18_b;
   
   for ii=1:length(type_idx)
   
    idx=find(subtype==sid(type_idx(ii)));
    
    if isempty(idx)
        % borrow the "neighboring" scene type
        % the preference is [Psi_radiance, Diff_cloud_Surf_T]
        % .I.e search Psi_radiance first, if a solution, then stop and use it;
        % if no solution, then search Diff_cloud_Surf_T.
        strsid = num2str(sid(type_idx(ii)));
        sidtmp(1) = str2num(strsid(1:4));
        sidtmp(2) = str2num(strsid(5:6));
        sidtmp(3) = str2num(strsid(7:end));
        for iscan = 3:-1:2
           for ip=1:6
                   sid2 = sidtmp;
                   sid2(iscan) = sidtmp(iscan) -ip;
                   if sid2(iscan)>0
                   if sid2(3)>=10
                     sidout = sid2(1)*10000 + sid2(2)*100 + sid2(3);
                   else
                     sidout = sid2(1)*1000 + sid2(2)*10 + sid2(3);
                   end
                   idx = find(subtype == sidout);
                   if ~isempty(idx)
                       break;
                   end
                   end
                  sid2 = sidtmp;
                  sid2(iscan) = sidtmp(iscan) +ip;
                    if sid2(3)>=10
                     sidout = sid2(1)*10000 + sid2(2)*100 + sid2(3);
                   else
                     sidout = sid2(1)*1000 + sid2(2)*10 + sid2(3);
                   end
                   idx = find(subtype == sidout);
                   if ~isempty(idx)
                       break;
                   end

          end
          if ~isempty(idx)
                       break;
          end
        end 
       
    end
    
    
    if ~isempty(idx)
    
       scene_id=idx;
       if scene_id > size(sid_ansR, 1)
          continue  % scene_id out of bounds (why???), should quality flag this
       end
       
       idw_valid = find(~isnan(rad_overcast(6:end,type_idx(ii))))+5;

      if length(idw_valid)>10  % Avoid all channels that have NaN values

       for ic=1:2
           if  ic==2
               rad_overcast(:,type_idx(ii)) = rad_overcast(:,type_idx(ii)) + rad_overcast_unc(:,type_idx(ii)); 
           end
            flux_est(idw2,1)=pi.*rad_overcast(idw2,type_idx(ii))./sid_ansR(scene_id,1:end-1)';
        
           % To get full spectral flux  %%%%%%%%%%%%%%%%
           rad2 = flux_est(idw_valid,1);
           new_mrad = pca{scene_id, 1}(idw_valid-5);   % -5 here because ADM stored are from channels 6-63
       
           new_pc = pca{scene_id, 2}(idw_valid-5, :);
 
           coeff = inv(new_pc' * new_pc) * new_pc' * (rad2' - new_mrad)';
           coeff = reshape(coeff, 1, length(coeff));
   
           % Project back to get the full spectral OLR
           %Now we add channel 64: 0-183 cm-1
           flux_est(6:64,1) =  pca{scene_id, 1}' + pca{scene_id, 2} * coeff';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            % Flux at two CO2 channels (17 and 18) are derived using regression coefficients
           % obtained by fitting radiances at neighbour channels
           if sensor_num ==1
               if ~isnan(rad_overcast(20,type_idx(ii)))
                  % don't do fitting
               else
                  flux_est(17,1) = ch17_b(1) + ch17_b(2) .* rad_overcast(19,type_idx(ii))' + ch17_b(3) .* rad_overcast(20,type_idx(ii))';
                  flux_est(18,1) = ch18_b(1) + ch18_b(2) .* rad_overcast(19,type_idx(ii))' + ch18_b(3) .* rad_overcast(20,type_idx(ii))';
               end
           elseif sensor_num == 2
               if ~isnan(rad_overcast(16,type_idx(ii)))
                  % don't do fitting
               else
                  flux_est(17,1) = ch17_b(1) + ch17_b(2) .* rad_overcast(16,type_idx(ii))' + ch17_b(3) .* rad_overcast(19,type_idx(ii))';
                  flux_est(18,1) = ch18_b(1) + ch18_b(2) .* rad_overcast(16,type_idx(ii))' + ch18_b(3) .* rad_overcast(19,type_idx(ii))';
               end
           end
           if ic==1
               flux_est_mean(1:64, type_idx(ii))= flux_est./wvdiff';
           elseif ic==2
               flux_est_unc(1:64, type_idx(ii))  = flux_est./wvdiff' - flux_est_mean(:, type_idx(ii));
           end
           
       end

      else

           flux_est_mean(1:64, type_idx(ii))= NaN;
           flux_est_unc(1:64, type_idx(ii))= NaN;
      end
    else
     
       disp('no such scene found');  
        
    end
   
   
   end
   
 
 end
