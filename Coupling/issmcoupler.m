function md=issmcoupler(md,nyrs,iscoupledfuture,coupling_interval,ncpus,isHO,modelstring,org,pathsmb,isconstantsmb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ISSM and provide output for SMB
%% coupled run
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
issmbyearly=1;
issmbsouthanomaly=0;
issmbnorthanomaly=1;

%% Time
if iscoupledfuture,
	if isconstantsmb,
		first_year=2100;
	else
		first_year=2020;
	end
else
	first_year=1959;
end



%%% COUPLING %%%

%% Loop over all years: simulate SMB and dynamics over time
for tt = 1:coupling_interval:nyrs
	disp(['***** tt = ' num2str(tt)])

	startyear=first_year+tt; 
	endyear=startyear+coupling_interval;

	%If we're hitting end of loop, run only the few final years remaining
	if nyrs-tt<coupling_interval
		coupling_interval=nyrs-tt;
		endyear=startyear+coupling_interval;
	end

	%If we have no years left to run, don't continue
	if startyear==endyear
		break
	end

	disp(['-----------   Coupling SMB-dynamics, year ' num2str(startyear) '-' num2str(endyear)]);

   %%% SMB: downscale reference SMB from 1 km to 100 m using new simulated surface topo
   % input: modelled gridded surface topo. at time tt (netcdf or mat, 32N); reference SMB at year tt to tt + coupling_interval
   % output: downscaled monthly SMB from year tt to tt + coupling_interval (netcdf 32N)
	disp(['-----------   Downscaling SMB, year ' num2str(startyear) '-' num2str(endyear)]);
	if iscoupledfuture,
		pyrunfile("./Coupling/Downscaling/downscale_smb_spline_future.py", "done",...
			start_year=startyear, end_year=endyear, path_smb=pathsmb, model_string=md.miscellaneous.name);
	else
		pyrunfile("./Coupling/Downscaling/downscale_smb_spline.py", "done",...
			start_year=startyear, end_year=endyear, path_smb=pathsmb, model_string=md.miscellaneous.name);
	end

	%Load ISSM model results from previous step in loop
	if tt==1
		%do nothing if we're starting from 1960/2020
	else
		if iscoupledfuture,
			%load results from previous years modelled
			md = loadmodel([modelstring 'TransientCoupledFuture']);
		else
			%load results from previous years modelled
			md = loadmodel([modelstring 'TransientCoupledHist']);
		end
		
		%Initialize geometry from previous transient results
		md=transientrestart(md);
	end

	%%%%% ISSM: run ice cap evolution for coupling_interval years %%%%
   % input: downscaled monthly SMB at year tt

	%% Setup updated SMB forcing from year tt to tt + coupling_interval
	%Set up forcing times
	smbForcingTimes=[1:coupling_interval];
	smb_nyrs=length(smbForcingTimes);



	disp(['-----------   Loading SMB rates at year ' num2str(startyear) '-' num2str(endyear)]);
		%%% Read downscaled 100m SMB model data %%%
		ncdata	=['./Coupling/Output/' md.miscellaneous.name 'smboutput.nc'];
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

			nmonths=coupling_interval*12;
			smb=NaN(md.mesh.numberofvertices,nmonths);

			size(smbdata)

			%Fill in SMB data for each month
			for i=1:nmonths
				smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
			end

			%Insert SMB as model forcing
			%setup arrays with model forcing times, with monthly forcing for each year
			forcingtimes=[];
			for i=1:coupling_interval %yrs
				for j=1:12
					forcingtime=i-1+j/12;
					forcingtimes=[forcingtimes startyear+forcingtime];
				end
			end

			%%YEARLY SMB FORCING
			if issmbyearly,
				forcingtimes=linspace(1,coupling_interval,coupling_interval);
				forcingtimes=startyear+forcingtimes;
				nyears=length(forcingtimes);

				%Extract smb for all months for all years
				%Calculate yearly mean for each year, for each vertex
				smb_sum_yearly=NaN(md.mesh.numberofvertices,nyears);
				%Go through all years
				for i=1:coupling_interval %yrs
					 %Sum mass balance for each year
					 for j=1:md.mesh.numberofvertices
						 smb_sum_yearly(j,i)= sum(smb(j,(i-1)*12+1:(i-1)*12+12));
			% 			 smb_sum_yearly(j,i)= sum(smb(j,year*12+1:year*12+12));
					 end
				end

				smb=smb_sum_yearly;
			end

			smb=smb/1000; %convert from mm to m w.e.
			smb=smb*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

% 		end

		%% If running commit, regrow or reverse, we're using a constant SMB. 
		if isconstantsmb,
			disp(['--- Computing mean SMB for all years'])
			%Take mean of all years, for each vertex
			smb_mean=NaN(md.mesh.numberofvertices,1);
			for i=1:md.mesh.numberofvertices
				 smb_mean(i,1)=mean(smb(i,:));
			end

			%insert new mean smb as the smb forcing
			smb=smb_mean;

			clear smb_mean
		end



	%%% Add spatial anomalies (manual bias-corrections)
	%SOUTH
	if issmbsouthanomaly,
		disp('--- Adding SMB anomaly in the south ---')
      in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostSouth.exp',1);
		smb_anomaly=-0.5;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);
	end

	%NORTH
	if issmbnorthanomaly,
		disp('--- Adding SMB anomaly in the North ---')
		in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostNorth.exp',1);
		smb_anomaly=0.75;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);
	end



	%Reset and fill in model smb
	if issmbyearly==1 %&& isconstantsmb==0
		md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,nyears);
	else %if issmbyearly==0
		md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,nmonths);
	end

% 	if isconstantsmb, %temporally constant SMB
% 		md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,1); %reset SMB
% 		md.smb.mass_balance=smb_mean; %only need one smb since it's the mean
% 		
% 		clear smb_mean

% 	else %temporally variable SMB
		for i=1:md.mesh.numberofvertices
			md.smb.mass_balance(i,:)=smb(i,:);
		end
	% 	md.smb.mass_balance(1:end-1,:)=smb; %insert SMB values as model forcing
		md.smb.mass_balance(end,:)=forcingtimes; %insert model forcing times (months)

		clear smb
% 	end


	%Update model time settings
	md.timestepping.start_time=startyear;
	md.timestepping.final_time=endyear;

% 	pos=find(isnan(md.smb.mass_balance(1:end-1,1)));
% 	md.smb.mass_balance(pos,1)=100;
% 	plotmodel(md,'data',md.smb.mass_balance(1:end-1,1),'layer',1,'expdisp','./Exp/Others/HenningJostDomainMedium202405d.exp')

	disp(['-----------   Running ISSM, year ' num2str(startyear) '-' num2str(endyear)]);
	% Run ISSM from year tt to tt + coupling_interval
	md=solve(md,'Transient');

	%Save results
	disp(['-----------   Saving ISSM model results at year ' num2str(endyear)]);
	savemodel(org,md);

	%%% ISSM: interpolate simulated ice geometry from mesh onto 100 m regular grid, write to netcdf4
	%% Interpolate simulated geometry to grid
	disp(['-----------   Interpolating modelled surface to grid at year ' num2str(endyear)]);
% 	srfMesh=md.results.TransientSolution3(end).Surface; %save modelled surface
	srfMesh=md.results.TransientSolution(end).Surface; %save modelled surface

	%% Make ice/ice-free mask
	iceMask=NaN(md.mesh.numberofvertices,1);
	pos=md.results.TransientSolution(end).Thickness>10;
	iceMask(pos)=1;
	pos=md.results.TransientSolution(end).Thickness<=10;
	iceMask(pos)=0;

	if isHO,
		srfMesh=project2d(md,srfMesh,1);
		iceMask=project2d(md,iceMask,1);
		srfGrid = InterpFromMeshToGrid(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,srfMesh,xdata,ydata,NaN);
		iceMaskGrid = InterpFromMeshToGrid(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,iceMask,xdata,ydata,NaN);
	else
% 		srfGrid = InterpFromMeshToGrid(md.mesh.x,md.mesh.y,srfMesh,xdata,ydata,NaN);
		srfGrid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,srfMesh,xdata,ydata,NaN);
	end

	%Meshgrid
	[Xdata,Ydata]=meshgrid(xdata,ydata);
	ny=length(Xdata(:,1));
	nx=length(Ydata(1,:));

	% contourf(X,Y,srfGrid)


	%% Write to netcdf
	disp(['-----------   Write modelled surface to netcdf at year ' num2str(endyear)]);
	filename=['./Coupling/Output/' md.miscellaneous.name 'issmoutput.nc'];
	out=exist(filename)
	if out==2
		disp('--Deleting existing modelled surface, and writing new one')
		delete(filename)
	end
	
	%Create file and set format
	nccreate(filename,'X','dimensions',{'X',1,nx})
	nccreate(filename,'Y','dimensions',{'Y',1,ny})
	nccreate(filename,'surface_modelled','dimensions',{'Y','X'})
	nccreate(filename,'ice_mask','dimensions',{'Y','X'})

	ncwrite(filename,'X',xdata) % add data to variables
	ncwrite(filename,'Y',ydata) % add data to variables
	ncwrite(filename,'surface_modelled',srfGrid) % add data to variables
	ncwrite(filename,'ice_mask',iceMaskGrid) % add data to variables

	S = ncinfo(filename); %Read variable, dimension, group definitions. This information defines the file's schema.
	S.Format = 'netcdf4';

end %coupling loop end



