function [successful_completion, error_ca] = produce_FLX_granulepart(varargin)
% Produce part of a 2B-FLX granule, based on the input parameters
%  specified by a 'pd' struct returned from configure_toplevel_IO.
%
% Returns a successful completion boolean flag and a cell array error code.

successful_completion = false; 

error_ca = {'#NONE#', 0};  % Default

pd = configure_toplevel_IO;  % Get/set shared filepaths, dirs, control params

% Read in product file data specifications (done here to have access to
%  fill_value parameters for the calculations below):
JSON_fspecs_FLX_fpath = fullfile(pd.ancillary_data_dir, ...
                                 'Flx_product_filespecs.json');
[FLX_schema, FLX_jdat] = init_nc_schema_from_JSON(JSON_fspecs_FLX_fpath, false);

% Define/initialize output data structure array:
odat = struct;
odat.global_atts = struct;
odat.Geometry = struct;
odat.Geometry.group_atts = struct;
odat.Geometry.group_dims = struct;
odat.Flx = struct;
%odat.Flx.group_atts = struct;
odat.Flx.group_dims = struct;

% Open input L1B granule data file, read in some useful parameters:
gg_att_C = netcdf.getConstant('NC_GLOBAL');
L1B.ncid = netcdf.open(pd.L1B_rad_fpath, 'NC_NOWRITE');
[L1B.n_dims, ~, L1B.n_globalatts, ~] = netcdf.inq(L1B.ncid);

L1B.G_gid = netcdf.inqNcid(L1B.ncid, 'Geometry');
L1B.R_gid = netcdf.inqNcid(L1B.ncid, 'Radiance');
dimid = netcdf.inqDimID(L1B.ncid, 'atrack');
[~, dims.max_nframes] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'xtrack');
[~, dims.nxtrack] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'UTC_parts');
[~, dims.nUTCparts] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'FOV_vertices');
[~, dims.nvertices] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'spectral');
[~, dims.nspectral] = netcdf.inqDim(L1B.ncid, dimid);

% Copy some L1B global and group attribute values to the output structure array:
g_atts_to_copy = {'granule_ID', 'spacecraft_ID', 'sensor_ID', ...
                  'ctime_coverage_start_s', 'ctime_coverage_end_s', ...
                  'UTC_coverage_start', 'UTC_coverage_end', ...
                  'orbit_sim_version'};
for ia=1:length(g_atts_to_copy)
   g_att_name = g_atts_to_copy{ia};
   g_att_value = netcdf.getAtt(L1B.ncid, gg_att_C, g_att_name);
   odat.global_atts.(g_att_name) = g_att_value;
end

gp_atts_to_copy = {'image_integration_duration_ms', 'solar_beta_angle_deg', ...
                   'TAI_minus_ctime_at_epoch_s', ...
                   'start_granule_edge_ctime_s', 'end_granule_edge_ctime_s'};
for ia=1:length(gp_atts_to_copy)
   gp_att_name = gp_atts_to_copy{ia};
   odat.Geometry.group_atts.(gp_att_name) = ...
                               netcdf.getAtt(L1B.G_gid, gg_att_C, gp_att_name);
end

% Determine instrument/sensor ID:
sensor_num = str2num(odat.global_atts.sensor_ID(5:6));

% (pd.idxbeg_atrack:pd.idxend_atrack) is the inclusive range/subset to process.
% NOTE that when pd.idxend_atrack == -1, the actual atrack dimension length
%  (determined from the L1B file, since it varies across files) should be used
%  for the end index of the subset.
%
% This means that to process the first N frames, we expect idxbeg_atrack to be 1
% and idxend_atrack to be N.
frame_beg = pd.idxbeg_atrack;
frame_end = pd.idxend_atrack;
if (frame_end == -1)
   % if idxend is -1, that means do the whole file.
   % the inclusive end frame is then the max N frame in the L1B file.
   frame_end = dims.max_nframes;
end
dims.nframes = frame_end-frame_beg+1;

% Read some Geometry group data from the L1B file:
fields_to_read = {'obs_ID', 'latitude', 'longitude'};
%@ TO-DO: make sure viewing_zenith_angle = 0 for each scene?
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.G_gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some Geometry group data from the L1B file:
fields_to_read = {'time_UTC_values'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.G_gid, fields_to_read, ...
                                 [1,frame_beg], [dims.nUTCparts,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some Radiance group data from the L1B file:

fields_to_read = {'wavelength', 'idealized_wavelength'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.R_gid, fields_to_read, ...
                  [1,1], [dims.nspectral,dims.nxtrack]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_read = {'spectral_radiance', 'spectral_radiance_unc', ...
                  'radiance_quality_flag'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.R_gid, fields_to_read, ...
                  [1,1,frame_beg], [dims.nspectral,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some meteorology data from the AUX-MET file:

AUXM.ncid = netcdf.open(pd.AUX_MET_fpath, 'NC_NOWRITE');
[AUXM.n_dims, ~, AUXM.n_globalatts, ~] = netcdf.inq(AUXM.ncid);
AUXM.gid = netcdf.inqNcid(AUXM.ncid, 'Aux-Met');
dimid = netcdf.inqDimID(AUXM.ncid, 'zlevels');
[~, dims.nzlevels] = netcdf.inqDim(AUXM.ncid, dimid);
dimid = netcdf.inqDimID(AUXM.ncid, 'n_igbp_classes');
[~, dims.nIGBPc] = netcdf.inqDim(AUXM.ncid, dimid);

fields_to_read = {'surface_pressure', 'skin_temp', 'seaice_concentration', ...
                  'snow_cover', 'total_column_wv'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

if pd.proc_mode == 1
   fields_to_read = {'merged_surface_type_prelim'};
   [error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);

   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
   merged_surface_type = AUXM.merged_surface_type_prelim;
end

fields_to_read = {'VIIRS_surface_type'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                     [1,1,frame_beg], [dims.nIGBPc,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end
[~, surface_type_xat] = max(AUXM.VIIRS_surface_type, [], 1);

fields_to_read = {'temp_profile', 'wv_profile'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                   [1,1,frame_beg], [dims.nzlevels,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_read = {'pressure_profile'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                                        [1], [dims.nzlevels]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(AUXM.ncid)  % Done with input AUX-MET file

% Possibly read some surface type data from the AUX-SAT file:
if pd.proc_mode == 2
   AUXS.ncid = netcdf.open(pd.AUX_SAT_fpath, 'NC_NOWRITE');
   [AUXS.n_dims, ~, AUXS.n_globalatts, ~] = netcdf.inq(AUXS.ncid);
   AUXS.gid = netcdf.inqNcid(AUXS.ncid, 'Aux-Sat');

   fields_to_read = {'merged_surface_type_final'};
   [error_ca, AUXS] = read_PREFIRE_granule(AUXS, AUXS.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
   merged_surface_type = AUXS.merged_surface_type_final;

   netcdf.close(AUXS.ncid)  % Done with input AUX-SAT file
end

% Read cloud product from the 2B-CLD file:

CLD.ncid = netcdf.open(pd.L2B_cld_fpath, 'NC_NOWRITE');
[CLD.n_dims, ~, CLD.n_globalatts, ~] = netcdf.inq(CLD.ncid);
CLD.gid = netcdf.inqNcid(CLD.ncid, 'Cld'); 

fields_to_read = {'cloudtop_pressure', 'cloud_tau', 'cloud_d_eff', ...
                  'cld_quality_flag'};
[error_ca, CLD] = read_PREFIRE_granule(CLD, CLD.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(CLD.ncid)  % Done with input 2B-CLD file

% Read some cloud mask data from the 2B-MSK file:

MSK.ncid = netcdf.open(pd.L2B_msk_fpath, 'NC_NOWRITE');
[MSK.n_dims, ~, MSK.n_globalatts, ~] = netcdf.inq(MSK.ncid);
MSK.gid = netcdf.inqNcid(MSK.ncid, 'Msk');

fields_to_read = {'cloud_mask'};
[error_ca, MSK] = read_PREFIRE_granule(MSK, MSK.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(MSK.ncid)  % Done with input 2B-MSK file

% Read useful data from SRF file:
SRF_fpath = fullfile(pd.ancillary_data_dir, ...
            sprintf('PREFIRE_TIRS%1d_%s.nc', sensor_num, pd.SRF_disambig_str));
ncid = netcdf.open(SRF_fpath, 'NC_NOWRITE');
fields_to_read = {'channel_wavelen1', 'channel_wavelen2'};
for iv=1:length(fields_to_read)
   v_name = fields_to_read{iv};
   varid = netcdf.inqVarID(ncid, v_name);
   SRF.(v_name) = netcdf.getVar(ncid, varid);  % C-indexing
end
SRF.im_version = netcdf.getAtt(ncid, gg_att_C, 'instrument_model_version');
netcdf.close(ncid)  % Done with input SRF file

% Read some reference emissivity fields:
fid = fopen(fullfile(pd.ancillary_data_dir, 'EMISSIVITY_id1'), 'r');
for i=1:18
    fgetl(fid);  % Read (and ignore) header line
end
for i=1:964
    str = fgetl(fid);
    data = str2num(str);
    wn(i) = data(1);
    emis(i,:) = data(2:end);
end
fclose(fid);
wn2 = 2:1999;
for i=1:18
    emis2(:,i) = interp1(wn,emis(:,i), wn2);
end

% Read clear-sky ADM info:
ADMstr = {'Sea-ice', 'Melted-ice', 'Ocean', 'Permanent-snow', 'Fresh-snow', ...
          'Non-snow-land'};
if sensor_num == 1
   ADM_path = fullfile(pd.ancillary_data_dir, ...
                       'ADM_clear/ADM_clear_SRF_v12_TIRS1/');
elseif sensor_num == 2
   ADM_path = fullfile(pd.ancillary_data_dir, ...
                       'ADM_clear/ADM_clear_SRF_v12_TIRS2/');
end
for iADMt=1:length(ADMstr)
   % ADMs stored are from PREFIRE channels 6 to 63
   ADM_fpfx = fullfile(ADM_path, ADMstr{iADMt});
   ivars = {'subtype', 'sid_ansR'};
   tmp_sa = load([ADM_fpfx '_ansR'], ivars{:});
   for ifn=1:length(ivars)
      ADM_clear{iADMt}.(ivars{ifn}) = tmp_sa.(ivars{ifn});
   end
   tmp_sa = load([ADM_fpfx '_PREFIRE_full_spectrum_PCA_normalized'], ...
                 'PREFIRE_full_spectrum_PCA');
   ADM_clear{iADMt}.pca = tmp_sa.PREFIRE_full_spectrum_PCA;
     % Regression coefficients for channel 17 and 18
   ivars = {'ch17_b', 'ch18_b'};
   tmp_sa = load([ADM_fpfx '_co2_channel_coef'], ivars{:});
   for ifn=1:length(ivars)
      ADM_clear{iADMt}.(ivars{ifn}) = tmp_sa.(ivars{ifn});
   end
end

% Read "cloudy" ADM info:
ADMstr = {'cloudy'};
if sensor_num == 1 
   ADM_path = fullfile(pd.ancillary_data_dir, ...
                       'ADM_overcast/ADM_cloudy_SRF_v12_TIRS1/');
elseif sensor_num == 2 
   ADM_path = fullfile(pd.ancillary_data_dir, ...
                       'ADM_overcast/ADM_cloudy_SRF_v12_TIRS2/');
end
for iADMt=1:length(ADMstr)
   % ADMs stored are from PREFIRE channels 6 to 63
   ADM_fpfx = fullfile(ADM_path, ADMstr{iADMt});
   ivars = {'subtype', 'sid_ansR'};
   tmp_sa = load([ADM_fpfx '_ansR'], ivars{:});
   for ifn=1:length(ivars)
      ADM_cloudy{iADMt}.(ivars{ifn}) = tmp_sa.(ivars{ifn});
   end
   tmp_sa = load([ADM_fpfx '_PREFIRE_full_spectrum_PCA_normalized'], ...
                 'PREFIRE_full_spectrum_PCA');
   ADM_cloudy{iADMt}.pca = tmp_sa.PREFIRE_full_spectrum_PCA;
     % Regression coefficients for channel 17 and 18
   ivars = {'ch17_b', 'ch18_b'};
   tmp_sa = load([ADM_fpfx '_co2_channel_coef'], ivars{:});
   for ifn=1:length(ivars)
      ADM_cloudy{iADMt}.(ivars{ifn}) = tmp_sa.(ivars{ifn});
   end
end

c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;
R = 287.05;
g = 9.806;
mm_to_cm = 0.1;  % mm -> cm unit conversion factor

%-- Some constraints about which FOVs a retrieval will be attempted for:

% Ch 14 is good on both SAT1 and SAT2; for certain quality-check operations,
%  this channel is selected in order to yield a boolean array of the correct
%  shape that reflects the FOV (i.e., cross-track + along-track) dependence of
%  a field.
qc_ch = 14;

limit_by_latitude = true;  % Limit processing by latitude?
latlim = 60.;  % latitude range check is okay if (abs(L1B.latitude) > latlim)

rqflim = 1;  % Any obs with a radiance quality flag value > this is not used

% These are per FOV (dimensioned {dims.nxtrack,dims.nframes}):
% (note that AUX-MET info is created for every FOV present in a 1B-RAD granule,
%  and so does not have one of these boolean masks)
if limit_by_latitude
   lat_okay = abs(L1B.latitude) > latlim;
else
   lat_okay = L1B.latitude > -999;  % Always evaluates to an array of 'true'
end
cmsk_okay = ( MSK.cloud_mask >= 0 );  % Simply a check for missing values
L1B_FOV_okay = squeeze(L1B.radiance_quality_flag(qc_ch,:,:) <= rqflim);
cld_okay = ( (CLD.cld_quality_flag <= 3 & MSK.cloud_mask > 1) | ...
             (MSK.cloud_mask <= 1) );  % Should also be true for clear sky
cld_q_better = ( (CLD.cld_quality_flag <= 1 & MSK.cloud_mask > 1) | ...
             (MSK.cloud_mask <= 1) );  % Should also be true for clear sky
calc_for_this_FOV = ( lat_okay & cmsk_okay & cld_okay );

% Change longitude to the range of 0-360 deg:
lon_0to360 = L1B.longitude;
lonind = ( L1B.longitude < 0 );
lon_0to360(lonind) = L1B.longitude(lonind)+360;

% Mask unwanted spectral_radiance* field elements:
b_idx = ( L1B.radiance_quality_flag > rqflim );
L1B.spectral_radiance(b_idx) = NaN;
L1B.spectral_radiance_unc(b_idx) = NaN;

flux = NaN(dims.nspectral, dims.nxtrack, dims.nframes);
flux_unc = NaN(dims.nspectral, dims.nxtrack, dims.nframes);
flux_tmp = NaN(dims.nspectral+1, dims.nxtrack, dims.nframes);

flux_flag = zeros(dims.nxtrack, dims.nframes, 'int8');
flux_flag(:,:) = FLX_jdat.Flx.flx_quality_flag.fill_value;

atmpropv_okay = false(dims.nxtrack, dims.nframes);

for j=1:dims.nxtrack  
   
   rad_flag = reshape(L1B.radiance_quality_flag(:,j,:), dims.nspectral, ...
                      dims.nframes);

   use_this_frame = calc_for_this_FOV(j,:);

   cldmask = MSK.cloud_mask(j,:);
   lat = L1B.latitude(j,:)';
   lon = L1B.longitude(j,:)';

   rad_obs = reshape(L1B.spectral_radiance(:,j,:), dims.nspectral, ...
                     dims.nframes);
   rad_obs_unc = reshape(L1B.spectral_radiance_unc(:,j,:), dims.nspectral, ...
                         dims.nframes);

   Ps = AUXM.surface_pressure(j,:);
   Ts = AUXM.skin_temp(j,:);
   sfc = surface_type_xat(1,j,:);
   merged_sfctype = merged_surface_type(j,:);

   ice = AUXM.seaice_concentration(j,:);  % fraction with range 0-1
   snow = AUXM.snow_cover(j,:);  % fraction with range 0-1
   TCWV = AUXM.total_column_wv(j,:)*mm_to_cm;  % convert to [cm]

   % Replace/clip some values, if present:
   ice(isnan(ice)) = 0;
   snow(isnan(snow)) = 0;
   snow(snow > 1) = 1;

   for k=1:dims.nframes
      T = AUXM.temp_profile(:,j,k)';
      q = AUXM.wv_profile(:,j,k)';  % g/kg 
      pres = AUXM.pressure_profile(:); 
      Ttmp = interp1(log(pres), T, log(Ps(k) - 300));
      lapse(k) = Ts(k) - Ttmp;
   end

   % for clear-sky footprints: 0, clear; 1, likely clear
   idx = find(cldmask' <= 1 & use_this_frame');
   if ~isempty(idx)
      [flux0 flux_unc0] = rad2flux_PREFIRE_clrsky(lat(idx), lon(idx), ...
        TCWV(idx), Ts(idx), lapse(idx), ice(idx), snow(idx), rad_obs(:,idx), ...
        rad_obs_unc(:,idx), ADM_clear, SRF.channel_wavelen1(j,:), ...
        SRF.channel_wavelen2(j,:), sfc(idx), sensor_num, merged_sfctype(idx));
      % Add 50-183 cm-1 band to flux_tmp; flux and flux_unc still have
      %  63 channels
      n_idx = length(idx);
      flux_tmp(:,j,idx) = reshape(flux0,dims.nspectral + 1, 1, n_idx);
      flux(:,j,idx) = reshape(flux0(1:dims.nspectral,:), ...
                              dims.nspectral, 1, n_idx);
      flux_unc(:,j,idx) = reshape(flux_unc0(1:dims.nspectral,:), ...
                                  dims.nspectral, 1, n_idx);

      flux_flag(j,idx) = 0;
      atmpropv_okay(j,idx) = true;
   end

   % for overcast footprints: 2, uncertain; 3, likely cloud; 4, cloud
   for k=1:dims.nframes

      if cldmask(k) > 1 & use_this_frame(k)
          T = AUXM.temp_profile(:,j,k)';
          pres = AUXM.pressure_profile(:)';
          CTP(k) = CLD.cloudtop_pressure(j,k);
          cldopt(k) = CLD.cloud_tau(j,k);
          re_in(k) = CLD.cloud_d_eff(j,k)*0.5;  % convert diameter to radius

          ri_in(k) = re_in(k);
          sfctmp =  sfc(k);  % IGBP surface type from AUX-MET
          mgd_sfctype = merged_sfctype(k);  % merged surface/supersurface type
                                            % from AUX-MET or AUX-SAT
                    
          if CTP(k) > 440  % [hPa]
             wratio(k) = 0.0;  % all ice cloud
             lwp(k) = cldopt(k)/(0.003448 + 2.431/ri_in(k)); % g/m2
          else 
             wratio(k) = 1;    % all water cloud
             lwp(k) = cldopt(k)/3*re_in(k)*2; % g/m2
          end
                     
          if (CTP(k) >= pres(1) & CTP(k) <= pres(end)) & ...
               (cldopt(k) > 0 & cldopt(k) < 1e5) & ...
               (re_in(k) > 0 & re_in(k) < 1e3)
             atmpropv_okay(j,k) = true;
             CTT(k) = interp1(log(pres), T, log(CTP(k)));

             if mgd_sfctype <= 3  % Water (excludes perennial ice shelves)
                sfctmp2(k) = 11;  % Default
                if ice(k) < 0.2
                        sfctmp2(k) = 13;
                elseif ice(k) < 0.4
                        sfctmp2(k) = 14;
                elseif ice(k) < 0.6
                        sfctmp2(k) = 15;
                elseif ice(k) < 0.8
                        sfctmp2(k) = 16;
                elseif ice(k) < 0.95
                        sfctmp2(k) = 17;
                elseif ice(k) <= 1
                        sfctmp2(k) = 18;
                end

             else  % Land or perennial ice shelf
                if sfctmp <= 2  % Evergreen needleleaf and broadleaf forest
                   sfctmp2(k) = 1;
                elseif sfctmp <= 5  % deciduous or mixed forest
                   sfctmp2(k) = 2;
                elseif sfctmp == 6  % closed shrublands
                   sfctmp2(k) = 3;
                elseif sfctmp == 7  % open shrublands
                   sfctmp2(k) = 4;
                elseif sfctmp <= 9  % nominal and woody savannas
                   sfctmp2(k) = 5;
                elseif sfctmp <= 12  % grassland, permanent wetland, cropland
                   sfctmp2(k) = 6;
                elseif sfctmp == 13  % urban and built-up land
                   sfctmp2(k) = 7;
                elseif sfctmp == 14  % mosiacs of variously-vegetated land
                   sfctmp2(k) = 8;
                elseif sfctmp == 15  % land with perennial snow/ice cover
                   sfctmp2(k) = 9;
                elseif sfctmp == 16  % barren land
                   sfctmp2(k) = 10;
                elseif sfctmp==18  % ???? Not a VIIRS IGBP surface type
                   sfctmp2(k) = 12;
                else
                   sfctmp2(k) = 6;  % Poorly-classified as 17=water?
                                    % Treat as if a common coastal land-class
                end
                  
                % surface type in snow region
                if snow(k) > 0.5
                   sfctmp2(k) = 8; % use median snow
                end
             end

             prad=planck(wn2,Ts(k));
             ec(k) = sum(prad'.*emis2(:,sfctmp2(k)))/sum(prad);
          end
      end
   end

   % Uses "atmpropv_okay" because some 2B-CLD fields may be invalid for use in
   %  this application
   idx = find(cldmask' > 1 & use_this_frame' & atmpropv_okay(j,:)');
   if ~isempty(idx)
      cldfrac = ones(size(idx));  % No better way to set this currently
      [flux0 flux_unc0] = rad2flux_PREFIRE_overcast(lat(idx), ...
                         lon(idx), TCWV(idx), Ts(idx), ice(idx), snow(idx), ...
                         cldfrac, CTT(idx), cldopt(idx), ...
                         wratio(idx), re_in(idx), ri_in(idx), lwp(idx), ...
                         ec(idx),  rad_obs(:,idx), rad_obs_unc(:,idx), ...
                         ADM_cloudy, SRF.channel_wavelen1(j,:), ...
                         SRF.channel_wavelen2(j,:), sensor_num);
      % Add 50-183 cm-1 band to flux_tmp; flux and flux_unc still have
      %  63 channels
      n_idx = length(idx);
      flux_tmp(:,j,idx) = reshape(flux0, dims.nspectral+1, 1, n_idx);
      flux(:,j,idx) = reshape(flux0(1:dims.nspectral,:), ...
                              dims.nspectral, 1, n_idx);
      flux_unc(:,j,idx) = reshape(flux_unc0(1:dims.nspectral,:), ...
                                  dims.nspectral, 1, n_idx);

      flux_flag(j,idx) = 1;
      clear cldfrac;
   end

end

% adjust to 20-km reference level because GCMs assume parallel instead of
%  spherical
% r should be adjusted for 100-0km, 0km-20km, the latter is done for correction
%  factor
r = (6371+100)^2/(6391)^2;
flux = flux .* r;
flux_tmp = flux_tmp .* r;

olr = reshape(sum(flux_tmp(1:end-1,:,:) .* repmat( ...
        (SRF.channel_wavelen2-SRF.channel_wavelen1)', [1,1,dims.nframes]), ...
                   1, "omitnan")+flux_tmp(end,:,:), ...
                              size(flux,2), size(flux,3));
clear flux_tmp;

%=== NOTE: The order of the following quality checks/assignments is important:

qc_bitflags = zeros(dims.nxtrack, dims.nframes, 'uint16');

% If any OLR is basically zero in the calculation above, that means all that
%  FOV's spectral fluxes were NaN, so that FOV should be flagged with FillValue
%  (i.e., invalid).
idx = olr < 1.e-20;
flux_flag(idx) = FLX_jdat.Flx.flx_quality_flag.fill_value;
olr(idx) = FLX_jdat.Flx.olr.fill_value;

% Set any values that are NaN or otherwise invalid to the proper output
%  missing value:
idx = find(isnan(olr));
olr(idx) = FLX_jdat.Flx.olr.fill_value;

idx = find(isnan(flux));
flux(idx) = FLX_jdat.Flx.spectral_flux.fill_value;
flux_unc(idx) = FLX_jdat.Flx.spectral_flux_unc.fill_value;

% Set some "why retrieval was not attempted" bitflags:
idx = find(~lat_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 0+1, 'uint16');
idx = find(~L1B_FOV_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 1+1, 'uint16');
idx = find(~cmsk_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 2+1, 'uint16');
idx = find(~cld_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 3+1, 'uint16');
idx = find(~atmpropv_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 4+1, 'uint16');
idx = find(~cld_q_better);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 5+1, 'uint16');

% Fill flx_quality_flag and flx_qc_bitflags output fields;
odat.Flx.flx_quality_flag = int8(flux_flag);
odat.Flx.flx_qc_bitflags = qc_bitflags;

% Load any additional data (including dimension lengths) into output structure
%  array(s):

odat.global_atts.full_versionID = pd.product_fullversion;

line_parts = strsplit(pd.product_fullversion, '_');
odat.global_atts.archival_versionID = strrep(line_parts{1}, 'R', '');

tmp_fpath = fullfile(pd.top_path, 'dist', ...
                     sprintf('prdgit_version_m%d.txt', pd.proc_mode));
fid = fopen(tmp_fpath, 'rt');
tmp_line = fgetl(fid);
fclose(fid);
idx = strfind(tmp_line, '(');
odat.global_atts.provenance = sprintf('%s%s (%s', ...
                   extractBefore(tmp_line, idx(1)), pd.product_fullversion, ...
                   extractAfter(tmp_line, idx(1)));

iftmp = {pd.L1B_rad_fpath, pd.AUX_MET_fpath, pd.L2B_msk_fpath, ...
         pd.L2B_cld_fpath};
if (pd.proc_mode == 2)
   iftmp{5} = pd.AUX_SAT_fpath;
end
for ia=1:length(iftmp)
   [~, name, ext] = fileparts(iftmp{ia});
   intmp{ia} = [name ext];
end
odat.global_atts.input_product_files = strjoin(intmp, ', ');

odat.global_atts.SRF_NEdR_version = SRF.im_version;

tmp_fpath = fullfile(pd.top_path, 'VERSION.txt');
fid = fopen(tmp_fpath, 'rt');
odat.global_atts.processing_algorithmID = fgetl(fid);
fclose(fid);

odat.Geometry.group_dims.atrack = dims.nframes;
odat.Geometry.group_dims.xtrack = dims.nxtrack;
odat.Geometry.group_dims.UTC_parts = dims.nUTCparts;
odat.Geometry.group_dims.FOV_vertices = dims.nvertices;

fields_to_copy = {'obs_ID', 'latitude', 'longitude', 'time_UTC_values'};
for i=1:length(fields_to_copy)
    odat.Geometry.(fields_to_copy{i}) = L1B.(fields_to_copy{i});
end

fields_to_readcopy = {'ctime', 'ctime_minus_UTC', 'subsat_latitude', ...
                      'subsat_longitude', 'sat_altitude', ...
                      'sat_solar_illumination_flag', 'orbit_phase_metric', ...
                      'satellite_pass_type'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                                 [frame_beg], [dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_readcopy = {'vertex_latitude', 'vertex_longitude', ...
                      'maxintgz_verts_lat', 'maxintgz_verts_lon'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                  [1,1,frame_beg], [dims.nvertices,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_readcopy = {'land_fraction', ...
                   'elevation', 'elevation_stdev', 'viewing_zenith_angle', ...
                   'viewing_azimuth_angle', 'solar_zenith_angle', ...
                   'solar_azimuth_angle', 'solar_distance', ...
                   'geoloc_quality_bitflags'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                                 [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(L1B.ncid);  % Done with input L1B file

odat.Flx.group_dims.atrack = dims.nframes;
odat.Flx.group_dims.xtrack = dims.nxtrack;
odat.Flx.group_dims.spectral = dims.nspectral;

odat.Flx.idealized_wavelength = L1B.idealized_wavelength;
odat.Flx.wavelength = L1B.wavelength;

odat.Flx.spectral_flux = flux;
odat.Flx.spectral_flux_unc = flux_unc;

odat.Flx.olr = olr;

now_UTCn_DT = datetime('now', 'TimeZone', 'UTC', 'Format', ...
                      'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');  % faux-UTC (no leap-s)
odat.global_atts.UTC_of_file_creation = string(now_UTCn_DT);

% UTC shape is (7, n), where the 7th row is millisec.
% change to a UTC array shaped (6, n) where the 6th row is
% fractional seconds.
UTC_granule_start = L1B.time_UTC_values(1:6,1);  % int16 array

% Construct output filepath, inserting a suffix to describe the
% 'granule part' (as an inclusive atrack range.)
fn_fmtstr = 'raw-PREFIRE_SAT%1d_2B-FLX_%s_%04d%02d%02d%02d%02d%02d_%s';
output_fn = sprintf(fn_fmtstr, sensor_num, pd.product_fullversion, ...
                    UTC_granule_start(1:6,1), odat.global_atts.granule_ID);
suffix_fmtstr = '-%s_%05d_%05d_of_%05df.nc';

% Output frame range (in filename) should be the original inclusive
%  0-index values.
if pd.idxend_atrack == -1
   iend = frame_end-1;  % 0-based index
else
   iend = pd.idxend_atrack-1;  % 0-based index
end
part_suffix = sprintf(suffix_fmtstr, 'atrack', pd.idxbeg_atrack-1, iend, ...
                      dims.max_nframes);
output_fpath = fullfile(pd.output_dir, [output_fn part_suffix]);

[~, name, ~] = fileparts(output_fpath);
odat.global_atts.file_name = name;

% Prepare for NetCDF-format file write, then define and write the file:
[m_ncs, vars_to_write] = modify_nc_schema(FLX_schema, odat);
if isfile(output_fpath)
   delete(output_fpath);
end
ncwriteschema(output_fpath, m_ncs);
write_predef_nc_vars(output_fpath, vars_to_write);

if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

successful_completion = true;

end
