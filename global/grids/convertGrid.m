% convert a NetCDF file into a CF complaint GRIDSPEC format for ESMF
% 
% Usage:  convertGrid(filename, latname, lonname, cornerlatname, cornerlonname)
%   where filename is the input grid file name, the new bound variables will be written back to the same
%   file.  latname and lonname are the coordinate variable names for the cell center vertices. 
%   cornerlatname and cornerlonname are the coordinate varilabe names for the corner vertices. We will 
%   construct the bound variables using the corner vertex coordinates and add bounds attribute into the
%   cell center vertex variables.
%
% Example: convertGrid('grid_spec_sub.nc','gridlat_t','gridlon_t','gridlat_vert_t','gridlon_vert_t')
%   The bound variable names are hard coded to 'gridlat_bound' and 'gridlon_bound', it can be passed in
%   as arguments as well
%   Note:  the new variables will be written into the same file.  So, make sure you save a copy of the original
%   file first.
%
% History:
%  10/22/2013: created by Peggy Li
%
function convertGrid(filename, latname, lonname, cornerlatname, cornerlonname)

%filename = 'grid_spec_sub.nc';
%latname = 'gridlat_t';
%lonname = 'gridlon_t';
%cornerlatname = 'gridlat_vert_t';
%cornerlonname = 'gridlon_vert_t';

boundlatname = 'gridlat_bound';
boundlonname = 'gridlon_bound';

nc = netcdf.open(filename, 'NC_WRITE');
clatvid = netcdf.inqVarID(nc, cornerlatname);
lats = netcdf.getVar(nc, clatvid);
clonvid = netcdf.inqVarID(nc, cornerlonname);
lons = netcdf.getVar(nc, clonvid);

netcdf.reDef(nc)
% add bounds attribute for latitude
vid = netcdf.inqVarID(nc, latname);
[name, xtype, dimids, natts]=netcdf.inqVar(nc, vid);
netcdf.putAtt(nc, vid, 'bounds', boundlatname);
% define bound variable for latitude
% the code only support 1D latitude variable for now
if (size(dimids,1)==1);
   % define a bound dimension first
   bounddim=netcdf.defDim(nc,'bound',2);
   dims2d(2)=dimids(1);
   [~,ydim] = netcdf.inqDim(nc, dimids(1));
   dims2d(1)=bounddim
   blatvid=netcdf.defVar(nc,boundlatname,xtype,dims2d);
   netcdf.copyAtt(nc,clatvid,'units',nc,blatvid);
else
   display('This program does not support 2D lat/lon coordinates yet');
   exit;
end

% add bounds attribute for longitude
vid = netcdf.inqVarID(nc, lonname);
[name, xtype, dimids, natts]=netcdf.inqVar(nc, vid);
netcdf.putAtt(nc, vid, 'bounds', boundlonname);
% define bound variable for longitude
% the code only support 1D longitude variable for now
if (size(dimids,1)==1) 
   % define a new dimension
   dims2d(2)=dimids(1);
   [~,xdim] = netcdf.inqDim(nc, dimids(1));
   dims2d(1)=bounddim;
   blonvid=netcdf.defVar(nc,boundlonname,xtype,dims2d);
   netcdf.copyAtt(nc,clonvid,'units',nc,blonvid);
else
   display('This program does not support 2D lat/lon coordinates yet');
   exit;
end
netcdf.endDef(nc);

% write out the bound data 
data=zeros(2,ydim);
data(1,:)=lats(1:ydim);
data(2,:)=lats(2:ydim+1);
netcdf.putVar(nc, blatvid, data);
clear data;
data=zeros(2,xdim);
data(1,:)=lons(1:xdim);
data(2,:)=lons(2:xdim+1);
netcdf.putVar(nc, blonvid, data);
%clear data;

netcdf.close(nc);
