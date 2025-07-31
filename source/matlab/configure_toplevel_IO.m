function [pd] = configure_toplevel_IO
% Assigns filepaths, directories/paths, and other control parameters to the
%  fields of a new structure array 'pd'

[func_path, ~, ~] = fileparts(which('configure_toplevel_IO'));
top_path = fullfile(func_path, '..', '..');  % has 'dist','source','test'
pd.top_path = top_path;

% These values are the current operationally-nominal ones, and will be used
%  unless the associated environment variables have been explicitly set:
pd.SRF_disambig_str = 'SRF_v13_2024-09-15';

if isunix

   tmp = getenv('PROC_MODE');
   if strlength(tmp) < 1  % relevant environment variables are not set
      pd.proc_mode = 1;

      pd.ancillary_data_dir = fullfile(top_path, 'dist/ancillary');

      pd.output_dir = fullfile(top_path, 'test/outputs');

      pd.product_fullversion = 'R01_P00';

      pd.L1B_rad_fpath = '../../test/inputs/PREFIRE_SAT1_1B-RAD_R01_P00_20241007075724_01877.nc';
      pd.AUX_MET_fpath = '../../test/inputs/PREFIRE_SAT1_AUX-MET_R01_P00_20241007075724_01877.nc';
      pd.L2B_msk_fpath = '../../test/inputs/PREFIRE_SAT1_2B-MSK_R01_P00_20241007075724_01877.nc';
      pd.L2B_cld_fpath = '../../test/inputs/PREFIRE_SAT1_2B-CLD_R01_P00_20241007075724_01877.nc';

      if (pd.proc_mode == 2)
          pd.AUX_SAT_fpath = 'test/inputs/PREFIRE_SAT1_AUX-SAT_R01_P00_20241007075724_01877.nc';
      end

      atrack_idx_range_0bi = 'ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END';
      tmp = strsplit(atrack_idx_range_0bi, ':');
      pd.idxbeg_atrack = str2num(tmp{2})+1;  % 1-based index
      if strcmp(tmp{3}, 'END')
          pd.idxend_atrack = -1;  % Sentinel value
      else
          pd.idxend_atrack = str2num(tmp{3})+1;  % 1-based index
      end

   else  % Use values obtained from environment variables
      pd = obtain_from_env_vars(pd, tmp);
   end

end


 function pd = obtain_from_env_vars(pd, procmode)
    pd.top_path = getenv('PACKAGE_TOP_DIR');
    pd.proc_mode = str2double(procmode);

    pd.ancillary_data_dir = getenv('ANCILLARY_DATA_DIR');
    pd.output_dir = getenv('OUTPUT_DIR');

    pd.product_fullversion = getenv('PRODUCT_FULLVER');

    tmp = getenv('SRF_DISAMBIG_STR');
    if strlength(tmp) > 0
       pd.SRF_disambig_str = tmp;
    end

    pd.L1B_rad_fpath = getenv('L1B_RAD_FILE');
    pd.AUX_MET_fpath = getenv('AUX_MET_FILE');
    pd.L2B_msk_fpath = getenv('L2B_MSK_FILE');
    pd.L2B_cld_fpath = getenv('L2B_CLD_FILE');

    if (pd.proc_mode == 2)
        pd.AUX_SAT_fpath = getenv('AUX_SAT_FILE');
    end

    tmp = strsplit(getenv('ATRACK_IDX_RANGE_0BI'), ':');
    pd.idxbeg_atrack = str2num(tmp{2})+1;  % 1-based index
    if strcmp(tmp{3}, 'END')
       pd.idxend_atrack = -1;  % Sentinel value
    else
       pd.idxend_atrack = str2num(tmp{3})+1;  % 1-based index
    end

 end

end
