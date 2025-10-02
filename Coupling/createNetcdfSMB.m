function createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)
%%%%%%%%%%%%%%%%%
% Write a .nc-file with mean SMB from a period of choice
%%%%%%%%%%%%%%%%

%Meshgrid
[Xdata,Ydata]=meshgrid(xdata,ydata);
ny=length(Xdata(:,1));
nx=length(Ydata(1,:));
nt=length(tdata);
% smbGridAll=NaN(ny,nx,nt);

% disp(['-----------   Interpolating smb from mesh to grid']);
% i=0;
% for t=tdata(1):tdata(end)
% 	i=i+1;
% 	if ismeancalc,
% 		smbMesh=md.results.TransientSolution(t).SmbMassBalance;
% 	else
% 		smbMesh=md.smb.mass_balance(1:end-1,i);
% 	end
% 	if isHO,
% 		smbMesh=project2d(md,smbMesh,1);
% 		smbGrid = InterpFromMeshToGrid(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,smbMesh,xdata,ydata,NaN);
% 	else
% 		smbGrid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,smbMesh,xdata,ydata,NaN);
% 	end
% 
% 	smbGridAll(:,:,i)=smbGrid;
% end

%% Calculate mean SMB over this period
% smbMean=NaN(ny,nx);
% for i=1:ny
% 	for j=1:nx
% 		smbMean(i,j)=mean(smbGridAll(i,j,:));
% 	end
% end
if ismeancalc,
	smbMean=mean(smbGridAll,3);
	contourf(Xdata,Ydata,smbMean)
end


%% Write to netcdf
disp(['-----------   Writing smb to netcdf']);
filename=[pathsmblong];
out=exist(filename)
if out==2
	disp('--Deleting existing file, and writing new one')
	delete(filename)
end

%Create file and set format
nccreate(filename,'X','dimensions',{'X',1,nx},'format','netcdf4')
nccreate(filename,'Y','dimensions',{'Y',1,ny},'format','netcdf4')

ncwrite(filename,'X',xdata) % add data to variables
ncwrite(filename,'Y',ydata) % add data to variables

S = ncinfo(filename); %Read variable, dimension, group definitions. This information defines the file's schema.
S.Format = 'netcdf4';

%Add time variable if needed
if ismeancalc,
	nccreate(filename,'mb_monthly','dimensions',{'Y','X'})
	ncwrite(filename,'mb_monthly',smbMean) % add data to variables
else

	nccreate(filename,'year_month','dimensions',{'year_month',nt},'datatype','string','format','netcdf4')
	nccreate(filename,'mb_monthly','dimensions',{'x','y','year_month'})

	%add data to netcdf variables
	ncwrite(filename,'year_month',tstrings) % add data to variables
	ncwrite(filename,'mb_monthly',smb) % add data to variables
end

S = ncinfo(filename); %Read variable, dimension, group definitions. This information defines the file's schema.
S.Format = 'netcdf4';
