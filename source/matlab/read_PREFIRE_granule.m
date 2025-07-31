function [error_ca, IO_sa] = read_PREFIRE_granule(IO_sa, gid, ...
                                                  fields_to_read, start, count)
% Reads in one or more variables from an already-open PREFIRE granule data file.
%
%-- Needed input parameters:
%   IO_sa :  (structure array) input, possibly altered, and then output
%   gid :  NetCDF-file group ID value (or file ID, if no group is relevant) of
%           already-open granule data file
%   fields_to_read :  cell array of names of variables to read from the granule
%                     (NOTE: per call of this routine, all requested variables
%                      should have the same dimensions)
%   start :  scalar start index, or array of start indices (the number of values
%            provided should match the number of dimensions of the requested
%            variables)
%   count :  scalar element count, or array of element counts (the number of
%            values provided should match the number of dimensions of the
%            requested variables)
%
%-- Output:
%   Various new or updated fields of the structure array IO_sa, with those
%    fields named the same as given in 'fields_to_read'

error_ca = {'#NONE#', 0};  % Default

enotvar_C = netcdf.getConstant('NC_ENOTVAR');
for iv=1:length(fields_to_read)
   v_name = fields_to_read{iv};

   % Determine varid:
   varid = netcdf.inqVarID(gid, v_name);
   if (varid == enotvar_C)
      msg = sprintf('ERROR: Variable %s does not exist in L1A file.', v_name);
      error_ca = {msg, 80};
      return
   end

   [noFillMode, fillValue] = netcdf.inqVarFill(gid, varid);
   if noFillMode
      IO_sa.(v_name) = netcdf.getVar(gid, varid, start-1, count);  % C-indexing
   else
      tmp_var = netcdf.getVar(gid, varid, start-1, count);  % C-indexing
      IO_sa.(v_name) = standardizeMissing(tmp_var, [fillValue NaN]);
      clear tmp_var;
   end
end

end
