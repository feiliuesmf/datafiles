% Adding bound variables to a NetCDF grid file to make it a CF complaint GRIDSPEC format for ESMF
% 
% Usage:  addBound(filename, latname, lonname, bdlatname, bdlonname)
%   where filename is the input grid file name, the new bound variables will be written back to the same
%   file.  latname and lonname are the coordinate variable names for the cell center vertices. 
%   bdlatname and bdlonname are the bound variable names for latitude and longitude, respectively.
%   The corner vertices will be calculated using the mid point of the two center grid points. To garantee
%   the continity across the longitude, the coordinates are converted into 3D Cartisian coordinates.
%
    % Example: convertGrid('global_gx3.grid.nc','ulat','ulon','ulat_bound','ulon_bound')
%
%   Note:  the new variables will be written into the same file.  So, make sure you save a copy of the original
%   file first.
%
% History:
%  4/10/2014: created by Peggy Li
%
% function convertGrid(filename, latname, lonname, boundlatname, boundlonname)

filename = 'fei_test.nc';
latname = 'ulat';
lonname = 'ulon';
boundlatname = 'ulat_bound';
boundlonname = 'ulon_bound';

deg2rad = pi/180.0;
rad2deg = 180.0/pi;

nc = netcdf.open(filename, 'NC_WRITE');
latvid = netcdf.inqVarID(nc, latname);
lats = netcdf.getVar(nc, latvid);
lonvid = netcdf.inqVarID(nc, lonname);
lons = netcdf.getVar(nc, lonvid);
units = netcdf.getAtt(nc, latvid, 'units');
if (strcmp(units,'radians'));
    % convert to degrees
    deglats = lats * rad2deg;
    deglons = lons * rad2deg;
    netcdf.putVar(nc, latvid, deglats);
    netcdf.putVar(nc, lonvid, deglons);
end

netcdf.reDef(nc)
% change units to degrees_east and degrees_north
%
if (strcmp(units,'radians'));
  netcdf.putAtt(nc, latvid,'units','degrees_north');
  netcdf.putAtt(nc, lonvid,'units','degrees_east');
end
% add bounds attribute for latitude
%[name, xtype, dimids, natts]=netcdf.inqVar(nc, latvid);
netcdf.putAtt(nc, latvid, 'bounds', boundlatname);
% define bound variable for latitude
% the code only support 1D latitude variable for now
dimids=zeros(2,1);
dimids(1)=netcdf.inqDimID(nc,'x');
dimids(2)=netcdf.inqDimID(nc,'y');
xtype = 'double';
if (size(dimids,1)==1) 
   % define a bound dimension first
   bounddim=netcdf.defDim(nc,'bound',2);
   dims2d(2)=dimids(1);
   [~,ydim] = netcdf.inqDim(nc, dimids(1));
   dims2d(1)=bounddim
   blatvid=netcdf.defVar(nc,boundlatname,xtype,dims2d);
else
   bounddim=netcdf.defDim(nc,'bound',4);
   dims3d(2)=dimids(1);
   dims3d(3)=dimids(2);
   [~,xdim] = netcdf.inqDim(nc, dimids(1));
   [~,ydim] = netcdf.inqDim(nc, dimids(2));
   dims3d(1)=bounddim
   blatvid=netcdf.defVar(nc,boundlatname,xtype,dims3d);
end

% add bounds attribute for longitude
%[name, xtype, dimids, natts]=netcdf.inqVar(nc, lonvid);
netcdf.putAtt(nc, lonvid, 'bounds', boundlonname);
% define bound variable for longitude
% the code only support 1D longitude variable for now
if (size(dimids,1)==1) 
   % define a new dimension
   dims2d(2)=dimids(1);
   [~,xdim] = netcdf.inqDim(nc, dimids(1));
   dims2d(1)=bounddim;
   blonvid=netcdf.defVar(nc,boundlonname,xtype,dims2d);
else
   dims3d(2)=dimids(1);
   dims3d(3)=dimids(2);
   [~,xdim] = netcdf.inqDim(nc, dimids(1));
   [~,ydim] = netcdf.inqDim(nc, dimids(2));
   dims3d(1)=bounddim
   blonvid=netcdf.defVar(nc,boundlonname,xtype,dims3d);
end
netcdf.endDef(nc);

% generate bound variables
if (~strcmp(units,'radians'));
   lats = lats * deg2rad;
   lons = lons * deg2rad;
end
% convert the coordinates into 3D cartesian
% first duplicate the first column at the end of x dimension and last column at the beginning of x dim
lon1 = cat(1, lons(end,:), lons, lons(1,:));
lat1 = cat(1, lats(end,:), lats, lats(1,:));
% add one row at the bottom and one row at the top

lastlat = lats(1,1)+(lats(1,1)-lats(1,2));
lat2 = zeros(xdim+2, 1);
lat2(:,1)=lastlat;
lat2 = cat(2, lat2(:,1), lat1);
lon2 = cat(2, lon1(:,1), lon1);
% the top layer is just the diplaced north pole
polelon = mean(lons(:,end));
polelat = mean(lats(:,end));
lat2(:,end+1) = polelat;
lon2(:,end+1) = polelon;

[X Y Z] = sph2cart(lon2, lat2, 1.0);
  
% average in X direcion
X1=(X(1:end-1,:)+X(2:end,:))/2;
Y1=(Y(1:end-1,:)+Y(2:end,:))/2;
Z1=(Z(1:end-1,:)+Z(2:end,:))/2;

% average in Y direction
X2 = (X1(:,1:end-1)+X1(:,2:end))/2;
Y2 = (Y1(:,1:end-1)+Y1(:,2:end))/2;
Z2 = (Z1(:,1:end-1)+Z1(:,2:end))/2;

[lon3 lat3 z3] = cart2sph(X2, Y2, Z2);

% change radian to deg
lon3 = lon3 * rad2deg;
lat3 = lat3 * rad2deg;
ind = find(lon3 < 0);
lon3(ind) = lon3(ind)+360;

% write out the bound data 
data=zeros(4,xdim,ydim);
data(1,:,:)=lat3(1:end-1,1:end-1);
data(2,:,:)=lat3(2:end,1:end-1);
data(4,:,:)=lat3(1:end-1,2:end);
data(3,:,:)=lat3(2:end,2:end);
netcdf.putVar(nc, blatvid, data);
data(1,:,:)=lon3(1:end-1,1:end-1);
data(2,:,:)=lon3(2:end,1:end-1);
data(4,:,:)=lon3(1:end-1,2:end);
data(3,:,:)=lon3(2:end,2:end);
netcdf.putVar(nc, blonvid, data);
%clear data;

netcdf.close(nc);
