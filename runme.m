%%%%%%%%%%%%%%%%%%%%%%
% Script to perform simulations of Jostedalsbreen ice cap
% written by Henning Åkesson (henning.akesson at geo.uio.no),
% 2022-2025
%%%%%%%%%%%%%%%%%%%%%%

%Define model steps
steps=[1 2 3 4 5 7];

% Steps:
% 1 - Mesh
% 2 - Parameterization
% 3 - Stressbalance
% 4 - Transient spinup
% 5 - Transient Perturb (historic 1960-2020)
% 6 - Transient coupled historic
% 7 - Transient coupled future
% 8 - Transient uncoupled future

%Naming & cpus
modelstring='Jost20251002a_';
ncpus=30;

%Future experiments
isparis1=1; %RCP4.5
isparis2=0;
isparis3=0;
isparis4=0;
iswarm1=0; %RCP8.5
iswarm2=0;
iswarm3=0;
iswarm4=0;
iscommitnow=0;
iscommitrcp45=0;
iscommitrcp85=0;
isreversercp45=0;
isreversercp85=0;
isregrownow=0;
isregrowrcp45=0;
isregrowrcp85=0;

%Coupling options
iscoupled=1;
iscoupledfuture=1;
coupling_interval=20; %yrs

%Time
nyrs=200; %simulation length for Transient step 4, spinup
nyrsPerturb=1; %not used
nyrsHistoric=61; %default 61 for historic
nyrsFuture=80;  %starts from 2020 (80yrs), or 2100 (200yrs for commit), 1000 for regrow
istimesteppingadaptive=0;

%meshing
ishighres=0; %special high-resolution mesh
isboya=0; %special domain not used here
isboyajost=0; %special domain not used here

%Historic or present-day run (always use = 1)
ishistoric=1;

%Spinup from no ice
isspinupnoice=0;

%Friction
isfrictioncomposite=1;
issecondinversion=0;
if isfrictioncomposite,
	isinversiontunsberg=1;
	isinversionnigard=1;
else
	disp('Using inversion everywhere!')
end
isfrictionshift=0; %manual shift of the friction field
isfrictionmosaic=1; %use multiple velocity datasets as observations

%Choose friction law
isbuddfriction=1;
ispismfriction=0;
isschooffriction=0;
isweertmanfriction=0;


%Ice approx.
isHO=1;
isMLHO=0; %force vertical shear function to be a polynomial of degree 4; the rest is standard HO
iscollapse=0; %collapse a HO 3d model into a SSA 2d model in TransientPerturb step, and run with SSA
if isHO==1; nlayers=4; end

modelpath='/uio/hypatia/geofag-felles/projects/jostice/issm/jost/Models/'; % modify to your own path

%%%% Model name %%%%
org = organizer('repository',modelpath,'prefix',[modelstring],'steps',steps);
%%%%%%%%%%%%%%%%%%%%

%Cluster
isbcgsolver=1; %1 means use solver biconjugate gradient with block Jacobi preconditioner (likely faster)
ismallocdebug=1; %1 means flag -malloc_debug activated


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if perform(org,'Mesh'),% {{{
    %Isotropic or anisotropic mesh?
    isisomesh=0;
    isrefinemeshbed=1; %create anisotropic mesh based on bed gradients only
    isrefinemeshvel=0; %create isomesh based on stress balance velocities
    isrefinevel=0;
    isrefinevelgradient=0; %refine mesh based on observed velocities
    isrefinenigard=1; %manually refine mesh in Nigardsbreen area
    isrefinetunsberg=1; %manually refine mesh in Tunsbergdalsbreen
    isrefineoldedalen=0; %manually refine mesh for Oldedalen outlets
    isrefineoldedalenausterdalen=1; %manually refine mesh for Oldedalen and Austerdalen outlet glaciers

    if isisomesh==1
        md=model;
        md=bamg(md,'domain','./Exp/Others/HenningNorthDomain.exp','hmax',1000);
    elseif isrefinemeshbed==1
        md=model;
		  md.miscellaneous.name=modelstring;
			if ishighres,
				md=bamg(md,'domain','./Exp/Others/HenningJostDomainMedium202408.exp','hmax',50); 
         else
				md=bamg(md,'domain','./Exp/Others/HenningJostDomainMedium202408.exp','hmax',150); % default
         end

        %Load bed
		  DEMname = './Data/bed/JOSTICE_bed_20243001mh_dtm10_merge.tif';
        currentDEM=DEMname;
        info = geotiffinfo(currentDEM);
        [Z,R] = geotiffread(currentDEM);
        Z=double(Z);

		R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
		[DEMx,DEMy] = worldGrid(R,'gridvectors');
        [xUTM,yUTM] = meshgrid(DEMx,DEMy);

        %Interpolate bed onto mesh
        md.geometry.bed=InterpFromGrid(DEMx,DEMy,Z,double(md.mesh.x),double(md.mesh.y));

        %Load thickness
		  DEMname = './Data/thickness/JOSTICE_thk_20243001_32N.tif';
        currentDEM=DEMname;
        info = geotiffinfo(currentDEM);

        [Z,R] = geotiffread(currentDEM);
        Z=double(Z);

		R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
		[DEMx,DEMy] = worldGrid(R,'gridvectors');
        [xUTM,yUTM] = meshgrid(DEMx,DEMy);

        %Interpolate onto mesh
        md.geometry.thickness=InterpFromGrid(DEMx,DEMy,Z,double(md.mesh.x),double(md.mesh.y));

        % Correct negative thicknesses
        pos=find(md.geometry.thickness<10);
        md.geometry.thickness(pos)=1;

        %Calculate surface
        md.geometry.surface=md.geometry.bed+md.geometry.thickness;

        %Set min and max mesh resolution
		  if ishighres,
            hmin=50;
				hmax=150;
        else
            hmin=100;
				hmax=300;
			end
				disp(['--- Creating mesh with resolution ' num2str(hmin) '-' num2str(hmax) ' m' ])

            hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);

		 if isrefinenigard,
				in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/HenningNigardFineMesh.exp',1);
            hmaxVertices(find(in))=100;
        end

        if isrefinetunsberg,
            in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OutlineTunsberg.exp',1);
            hmaxVertices(find(in))=100;
        end

        if isrefineoldedalen,
            in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OldedalenBasins.exp',1);
            hmaxVertices(find(in))=100;
        end

        if isrefineoldedalenausterdalen,
            in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OldedalenAusterdalenBasins.exp',1);
            hmaxVertices(find(in))=100;
        end

        %Save original mesh
        mdO=md;

        %Refine using gradient in bed topography
        md=bamg(md,'hmin',hmin,'hmax',hmax,'field',md.geometry.bed,...
            'gradation',1.5,'err',75,'anisomax',1.,'hmaxVertices',hmaxVertices);
        % 		'gradation',1.1,'err',75,'anisomax',1.,'hmaxVertices',hmaxVertices);

        %Reinterpolate observations onto the new mesh
        md.geometry.bed=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
            mdO.mesh.y,mdO.geometry.bed,md.mesh.x,md.mesh.y);
        md.geometry.surface=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
            mdO.mesh.y,mdO.geometry.surface,md.mesh.x,md.mesh.y);
        md.geometry.thickness=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
            mdO.mesh.y,mdO.geometry.thickness,md.mesh.x,md.mesh.y);

        %Adjust surface
        md.geometry.surface=max(md.geometry.surface,1);
        md.geometry.base=md.geometry.bed;
        md.geometry.thickness = md.geometry.surface-md.geometry.bed;
        disp('      -- Adjusting ice thickness');
        pos=find(md.geometry.thickness<=1);
        md.geometry.thickness(pos)=1;
        md.geometry.bed=md.geometry.surface-md.geometry.thickness;
        pos=find(md.geometry.base<md.geometry.bed);
        md.geometry.bed(pos) = md.geometry.base(pos);



        if isrefinevel,
			  %NOT IMPLEMENTED FOR JOSTEDALSBREEN, A GREENLAND EXAMPLE IS PROVIDED HERE
            disp('Refining mesh')
            %Get thickness from Bedmachine
            md.geometry.thickness = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'thickness');

            %Get velocities
            [md.inversion.vx_obs md.inversion.vy_obs]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
            pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));
            md.inversion.vx_obs(pos)=0;
            md.inversion.vy_obs(pos)=0;
            md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);

            %Define refinement everywhere where vel is fast (> vel_threshold)
            hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
            in=find(md.inversion.vel_obs>500); %find elemens where observed vel is fast
            hmaxVertices(in)=500;

            %Remesh again, adding refined region
            md=bamg(md,'hmin',hmin,'hmax',hmax,'field',...
                md.geometry.bed,...
                'gradation',1.1,'err',50,'anisomax',1.,'hmaxVertices',hmaxVertices);

            %Interpolate bed from netcdf observations onto the new mesh
            md.geometry.bed = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'bed');

            %  		plotmodel(md,'data',md.geometry.bed,'edgecolor','w',...
            %  			'colormap',demmap(220,-700,1500),'caxis',[-700 1500]);

        elseif isrefinemeshvel==1
				%NOT IMPLEMENTED

        elseif isrefinevelgradient, %refine mesh based on observed velocities
            disp('hey')
            %Set min and max mesh resolution
            hmin=50;
            hmax=500;
            hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
            in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'./Exp/Others/HenningJostDomainLarge.exp','node',1);
            hmaxVertices(find(in))=hmax;

            if isrefinenigard,
                in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/HenningNigardFineMesh.exp',1);
                hmaxVertices(find(in))=100;
            end

            %Save original mesh
            mdO=md;

            disp('   Reading velocities ');
            tiffname = './Data/velocity/millan_2022_velx_32N.tif';
            info = geotiffinfo(tiffname);
            [velObs,R] = geotiffread(tiffname);
            velObs=double(velObs);

				R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
				[TIFFx,TIFFy] = worldGrid(R,'gridvectors');
            [xUTM,yUTM] = meshgrid(TIFFx,TIFFy);

            %Interpolate velocity data onto mesh
            md.inversion.vx_obs=InterpFromGrid(TIFFx,TIFFy,velObs,double(md.mesh.x),double(md.mesh.y));

            tiffname = './Data/velocity/millan_2022_vely_32N.tif';
            info = geotiffinfo(tiffname);
            [velObs,R] = geotiffread(tiffname);
            velObs=double(velObs);

				R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
				[TIFFx,TIFFy] = worldGrid(R,'gridvectors');
            [xUTM,yUTM] = meshgrid(TIFFx,TIFFy);

            %Interpolate velocity data onto mesh
            md.inversion.vy_obs=InterpFromGrid(TIFFx,TIFFy,velObs,double(md.mesh.x),double(md.mesh.y));

            clear velObs TIFFx TIFFy xUTM yUTM

            pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));
            md.inversion.vx_obs(pos)=0;
            md.inversion.vy_obs(pos)=0;
            md.inversion.vel_obs  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
            md.initialization.vx  = md.inversion.vx_obs;
            md.initialization.vy  = md.inversion.vy_obs;
            md.initialization.vz  = zeros(md.mesh.numberofvertices,1);
            md.initialization.vel = md.inversion.vel_obs;

            %Refine using gradient in vel_obs
            md=bamg(md,'hmin',hmin,'hmax',hmax,'field',md.initialization.vel,...
                'gradation',1.4,'err',20,'anisomax',1.,'hmaxVertices',hmaxVertices);
            % 			'gradation',1.4,'err',20,'anisomax',1.,'hmaxVertices',hmaxVertices);

            %Reinterpolate observations onto the new mesh
            md.geometry.bed=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.bed,md.mesh.x,md.mesh.y);
            md.geometry.surface=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.surface,md.mesh.x,md.mesh.y);
            md.geometry.thickness=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.thickness,md.mesh.x,md.mesh.y);

            %Adjust surface
            md.geometry.surface=max(md.geometry.surface,1);
            md.geometry.base=md.geometry.bed;
            md.geometry.thickness = md.geometry.surface-md.geometry.bed;
            disp('      -- Adjusting ice thickness');
            pos=find(md.geometry.thickness<=1);
            md.geometry.thickness(pos)=1;
            md.geometry.bed=md.geometry.surface-md.geometry.thickness;
            pos=find(md.geometry.base<md.geometry.bed);
            md.geometry.bed(pos) = md.geometry.base(pos);

        elseif isboyajost, %refine area around Boyabreen, S. Jostedalsbreen
            %Set min and max mesh resolution
            hmin=50;
            hmax=5000;
            hVertices=NaN*ones(md.mesh.numberofvertices,1);
            in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/HenningJostBoya.exp',1);
            hVertices(find(in))=100;

            %Save original mesh
            mdO=md;

            %Refine mesh in Boyabreen area
            % 		md=bamg(md,'hmin',hmin,'hmax',hmax,'hVertices',hVertices);
            md=bamg(md,'hmin',hmin,'hmax',hmax,'field',md.geometry.bed,...
                'gradation',2,'err',75,'anisomax',1.,'hVertices',hVertices);

            %Reinterpolate observations onto the new mesh
            md.geometry.bed=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.bed,md.mesh.x,md.mesh.y);
            md.geometry.surface=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.surface,md.mesh.x,md.mesh.y);
            md.geometry.thickness=InterpFromMeshToMesh2d(mdO.mesh.elements,mdO.mesh.x,...
                mdO.mesh.y,mdO.geometry.thickness,md.mesh.x,md.mesh.y);

            %Adjust surface
            md.geometry.surface=max(md.geometry.surface,1);
            md.geometry.base=md.geometry.bed;
            md.geometry.thickness = md.geometry.surface-md.geometry.bed;
            disp('      -- Adjusting ice thickness');
            pos=find(md.geometry.thickness<=1);
            md.geometry.thickness(pos)=1;
            md.geometry.bed=md.geometry.surface-md.geometry.thickness;
            pos=find(md.geometry.base<md.geometry.bed);
            md.geometry.bed(pos) = md.geometry.base(pos);

        end

    end

    savemodel(org,md);
    % end

end%}}}
if perform(org,'Parameterization'),% {{{
    md = loadmodel(org,'Mesh');

    md = parameterize(md,'./Scripts/Jost.par');

    savemodel(org,md);
end%}}}
if perform(org,'StressBalance'),% {{{
    md = loadmodel(org,'Parameterization');

    %Set levelset options
    md.levelset.stabilization=1;
    md.mask.ice_levelset(pos)=+1; %set levelset for no-ice vertices

    %remove 0 in ice_levelset (advection will fail if used)
    md.mask.ice_levelset(find(md.mask.ice_levelset==0))=-1;

    %make it a signed distance
    md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    % Reset levelset boundary conditions
    md.levelset.spclevelset(find(md.mesh.vertexonboundary)) = md.mask.ice_levelset(find(md.mesh.vertexonboundary));
    % 	plotmodel(md,'data',md.mask.ice_levelset)

    	disp('Spc levelset on domain boundaries only')
    	md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1); %reset spclevelset to no constraints

    	% Reset levelset boundary conditions on domain boundary
    	pos = find(md.mesh.vertexonboundary==1);
    	md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);

    %Set BCs along domain boundaries
    idxB=find(md.mesh.vertexonboundary==1); %vertices on boundary
    md.masstransport.spcthickness=NaN(md.mesh.numberofvertices,1);

    if isHO,
    	%Extrude mesh vertically
    	md=extrude(md,nlayers,1);

    	%Make sure bed is below base
    	pos1=find(md.geometry.bed>md.geometry.base);
    	md.geometry.base(pos1)=md.geometry.bed(pos1);

    	%Recalculate surface
    	md.geometry.surface=md.geometry.base+md.geometry.thickness;

    	%Set flowequation
    	md = setflowequation(md,'HO','all');

    elseif isMLHO,
    	%Set flowequation
    	md = setflowequation(md,'MOLHO','all');

    	%Set BCs
    	md.stressbalance.spcvx_shear = md.stressbalance.spcvx;
    	md.stressbalance.spcvy_shear = md.stressbalance.spcvx;
    	md.stressbalance.spcvx_base = md.stressbalance.spcvx;
    	md.stressbalance.spcvy_base = md.stressbalance.spcvx;

    else
    	md = setflowequation(md,'SSA','all');
    end


    %% FRICTION
    if ispismfriction==0 && isschooffriction==0 && isweertmanfriction==0 %use Budd
    	%PERFORM INVERSION

    	%Activate m1qn3-type inversion
    	md.inversion=m1qn3inversion();

    	%Get velocities for the inversion
    	md.inversion.vx_obs=md.initialization.vx;
    	md.inversion.vy_obs=md.initialization.vy;
    	md.inversion.vz_obs=md.initialization.vz;
    	md.inversion.vel_obs=md.initialization.vel;

		if isfrictionmosaic,

			% Load Schellenberger
			tiffname = './Data/velocity/schellenberger_1996_merge_SN_x.tif';
			 info = geotiffinfo(tiffname);
			 [velObs,R] = geotiffread(tiffname);
			 velObs=double(velObs);
			 velObs=velObs*365.25; %convert m/day to m/year
			R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
			[TIFFx,TIFFy] = worldGrid(R,'gridvectors');
			 [xUTM,yUTM] = meshgrid(TIFFx,TIFFy);
			vx_obs_schell=InterpFromGrid(TIFFx,TIFFy,velObs,double(md.mesh.x),double(md.mesh.y));

			tiffname = './Data/velocity/schellenberger_1996_merge_SN_y.tif';
			 info = geotiffinfo(tiffname);
			 [velObs,R] = geotiffread(tiffname);
			 velObs=double(velObs);
			 velObs=velObs*365.25; %convert m/day to m/year
			R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
			[TIFFx,TIFFy] = worldGrid(R,'gridvectors');
			 [xUTM,yUTM] = meshgrid(TIFFx,TIFFy);
			vy_obs_schell=InterpFromGrid(TIFFx,TIFFy,velObs,double(md.mesh.x),double(md.mesh.y));

			%Fix no-data areas
			pos=find(vx_obs_schell<-5000);
			vx_obs_schell(pos)=md.inversion.vx_obs(pos);
			pos=find(vy_obs_schell<-5000);
			vy_obs_schell(pos)=md.inversion.vy_obs(pos);


	% 		%Save millan vel
	% 		vel_obs_millan=md.inversion.vel_obs;

			%Get velocities in high elevations and low surface slopes from Schell
			[sx,sy,s]=slope(md);
			sslope=averaging(md,s,10);
			pos=find(sslope<0.08);
			md.inversion.vx_obs(pos)=vx_obs_schell(pos);
			md.inversion.vy_obs(pos)=vy_obs_schell(pos);
			pos=find(md.geometry.surface>1775);
			md.inversion.vx_obs(pos)=vx_obs_schell(pos);
			md.inversion.vy_obs(pos)=vy_obs_schell(pos);
			md.initialization.vel  = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
			md.inversion.vel_obs=md.initialization.vel;

	% 		plotmodel(md,'data',vel_obs_millan-md.initialization.vel,'layer',1)
		end

    	%Control general
    	md.inversion.iscontrol=1;
    	md.inversion.maxsteps=40;
    	md.inversion.maxiter=3; %default 15

    	%Cost functions
    	md.inversion.cost_functions=[101 103 501];
    	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);

    	%Jost values
    	md.inversion.cost_functions_coefficients(:,1)=2e3; 
    	md.inversion.cost_functions_coefficients(:,2)=1;  
    	md.inversion.cost_functions_coefficients(:,3)=3.2e-6;

    	%Controls
    	md.inversion.control_parameters={'FrictionCoefficient'};
    	md.inversion.min_parameters=10*ones(md.mesh.numberofvertices,1); %default 10
    	md.inversion.max_parameters=400*ones(md.mesh.numberofvertices,1);


    	% If running inversion only for selected areas, activate those areas and exclude the rest
		if isfrictioncomposite,
    		md.inversion.cost_functions_coefficients(:,1:3) = 0; %reset cost function to zero everywhere

			if isinversiontunsberg,
				%Tunsberg
				in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OutlineTunsberg.exp',1);
				md.inversion.cost_functions_coefficients(find(in),1)=2e3; %Josh/Andreas P: 20, Henning: 2000
				md.inversion.cost_functions_coefficients(find(in),2)=1;  %Josh/Andreas P: 1
				md.inversion.cost_functions_coefficients(find(in),3)=3.2e-6; %Josh/Andreas P: 3.2e-6
			end

			if isinversionnigard,
				%Nigard
				in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OutlineNigard.exp',1);
				md.inversion.cost_functions_coefficients(find(in),1)=2e3; %Josh/Andreas P: 20, Henning: 2000
				md.inversion.cost_functions_coefficients(find(in),2)=1;  %Josh/Andreas P: 1
				md.inversion.cost_functions_coefficients(find(in),3)=3.2e-6; %Josh/Andreas P: 3.2e-6
			end
			if isinversiontunsberg==0 && isinversionnigard==0
				md.inversion.iscontrol=0; %deactivate inversion
			end

    	else

    	end

    	%Exclude areas with low surface slopes (typically at ice divides), where vel_obs can be high due to errors in obs data
		disp('---- Excluding areas with low surface slopes from inversion')
    	[sx,sy,s]=slope(md);
    	sslope=averaging(md,s,10);
    	pos=find(sslope<0.08); %default 0.1
    	md.inversion.cost_functions_coefficients(pos,1:3) = 0; %set cost function to zero for those areas, for all costfunctions

		disp('---- Excluding high-elevation areas from inversion')
    	pos=find(md.geometry.surface>1775);
    	md.inversion.cost_functions_coefficients(pos,1:3) = 0; %set cost function to zero for those areas, for all costfunctions

    	%Exlude areas with vel_obs < 1 (large relative errors) and thickness < 10 (no ice)
    	pos = find(md.geometry.thickness<10);
    	md.inversion.cost_functions_coefficients(pos,1:3) = 0; %set cost function to zero for those areas, for all costfunctions
    	pos = find(md.inversion.vel_obs<10);
    	md.inversion.cost_functions_coefficients(pos,1:3) = 0; %set cost function to zero for those areas, for all costfunctions

    	%Additional parameters
    	md.stressbalance.restol=0.01;
    	md.stressbalance.reltol=0.1;
    	md.stressbalance.abstol=NaN;

    	%%% INVERSION END %%%
    end

    %Dont use thermal model
    md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

    md.toolkits = toolkits();
    md.verbose=verbose('solution',false,'module',false,'convergence',false);


    %Settings for PISM friction
    if ispismfriction,
		 %NOT USED HERE, SEE https://doi.org/10.1038/s41467-022-29529-5

    elseif isschooffriction,
		 %NOT USED HERE, SEE https://doi.org/10.1038/s41467-022-29529-5
    	%Deactivate inversion
    	md.inversion.iscontrol=0;

    	%Solve
    	disp('Stress balance with inverted rheology_B plugged, using Schoof friction')
    	md=solve(md,'sb','np',4);

    else
        %Budd
        %if we have setup for inversion, not using pism or schoof or weertman friction, solve stressbalance inverting for friction
    	if isMLHO==1
    		disp('	Deactivating inversion, MLHO inversion not implemented yet')
    		md.inversion.iscontrol=0; %inversion active?
    	end

    	if md.inversion.iscontrol==1
    		disp('	Basal friction inversion active - hey')
			%First guess
%     		%%Set friction in areas not included in inversion
% 			%Find these areas:
			pos=find(md.inversion.cost_functions_coefficients(:,1)==0);
    		z_max=2000; %m asl
    		beta_max=200; %max friction coeff %default 400
    		md.friction.coefficient(pos) = (min(max(0,md.geometry.base(pos)+300),z_max)/z_max).*beta_max; %default 800

    	else %no inversion
    		%Adjust friction coeff depending on flow equation
    		if isHO==1
    			disp('	HO: Deactivating inversion, constructing friction coeffs dependent on bed elevation')
    			z_max=2000; %m asl
    			beta_max=300; %max friction coeff
    			md.friction.coefficient = (min(max(0,md.geometry.base+500),z_max)/z_max).*beta_max;
%     			md.friction.coefficient = (min(max(0,md.geometry.base+500),z_max)/z_max).*beta_max;
    		elseif isMLHO==1
    			disp('	MLHO: Deactivating inversion, constructing friction coeffs dependent on bed elevation')
    			z_max=2000; %m asl %default 2000
    			beta_max=800; %max friction coeff %default 200
				md.friction.coefficient = (min(max(0,md.geometry.base+500),z_max)/z_max).*beta_max;
                % 			md.friction.coefficient = (min(max(0,md.geometry.base+800),2000)/2000).*120; %BoknaHard parameters
% 					 md.friction.coefficient=300*ones(md.mesh.numberofelements,1);
    		else
    			%dont change, we're using SSA, friction coeff already set in .par file
    		end
    	end

% 		%Plotting friction field before inversion
% 		disp('---Plotting friction coefficient')
% 		if isHO,
% 			plotmodel(md,'data',md.friction.coefficient,'layer',1)
% 		else
% 			plotmodel(md,'data',md.friction.coefficient)
% 		end

    	%Solve
    	md.verbose=verbose('solution',false,'module',false,'convergence',false);

		disp(['	Solving on mimi using ' num2str(ncpus) ' cpus'])
		md.cluster=generic('name',oshostname,'np',ncpus);

		%Solve
		md=solve(md,'sb');


    	%Update model friction fields according to stress balance solution incl the inversion
    	if md.inversion.iscontrol,
    		md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

            % 		plotmodel(md,'data',md.friction.coefficient,'layer',1);
    		if issecondinversion, %run a second inversion, starting from the first one, doing only absolute vel misfit
    			md.inversion.cost_functions=[101];
    			md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
    			md.inversion.cost_functions_coefficients(:,1)=2000; %Josh/Andreas P: 20

    			%Exlude areas with vel_obs < 1 (large relative errors) and thickness < 10 (no ice)
    			pos = find(md.geometry.thickness<10);
    			md.inversion.cost_functions_coefficients(pos,:) = 0; %set cost function to zero for those areas, for all costfunctions
    			pos = find(md.inversion.vel_obs<1);
    			md.inversion.cost_functions_coefficients(pos,:) = 0; %set cost function to zero for those areas, for all costfunctions
    			%Solve
    			disp('--- Running second inversion with only absolute misfit')
    			md=solve(md,'sb');

    			%Update model friction fields from second inversion
    			md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

    		end
    	end

        % 	%Set friction for floating ice
    	pos=find(md.mask.ocean_levelset<0);
    	md.friction.coefficient(pos)=30;%md.inversion.min_parameters(1,1);
    	md.friction.p=ones(md.mesh.numberofelements,1);
    	md.friction.q=ones(md.mesh.numberofelements,1);
    end

		%Shift friction field manually
		if isfrictionshift,
			disp('---Shifting friction field manually')
			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OutlineTunsberg.exp',1);
			md.friction.coefficient(find(in))=md.friction.coefficient(find(in))+25;
			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OutlineNigard.exp',1);
			md.friction.coefficient(find(in))=md.friction.coefficient(find(in))+25;
% 			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OldedalenBasins.exp',1);
% 			md.friction.coefficient(find(in))=md.friction.coefficient(find(in))-100;
% 			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostSouth.exp',1);
% 			md.friction.coefficient(find(in))=md.friction.coefficient(find(in))-100;

			%Make sure no areas have low or negative friction
			pos=find(md.friction.coefficient<10);
			md.friction.coefficient(pos)=10;
		end


    %Make sure base is not below bed
    pos1=find(md.geometry.bed>md.geometry.base);
    md.geometry.base(pos1)=md.geometry.bed(pos1);

    %Recalculate surface
    md.geometry.surface=md.geometry.base+md.geometry.thickness;

    %Make sure bed = base for grounded ice
    pos1=find(md.mask.ocean_levelset>0);
    md.geometry.base(pos1)=md.geometry.bed(pos1);

    %Recalculate surface
    md.geometry.surface=md.geometry.base+md.geometry.thickness;

    % 	%Test plot
    % 	plotmodel(md,...
    % 		'data',md.inversion.vel_obs,'title','Observed velocity',...
    % 		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
    % 		'data',md.geometry.base,'title','Bed elevation',...
    % 		'data',md.results.StressbalanceSolution.FrictionCoefficient,'title','Friction Coefficient',...
    % 		'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
    % 		'caxis#1-2',([1.5,1500]),...
    % 		'colorbartitle#3','(m)', 'log#1-2',10);

    savemodel(org,md);

    % md=loadmodel([modelpath num2str(org.prefix) 'StressBalance']);


end%}}}
if perform(org,'Transient'),% {{{
    md = loadmodel(org,'StressBalance');

    %Load historic ice surface DEM and replace surface height for areas within the 1966 domain
    if ishistoric==1 %start from 1966 ice-surface DEM

        if isHO,
            %Collapse current model 3d -> 2d
            md=collapse(md);
        end

        %Load 1966 ice-surface DEM
        DEMname = './Data/surface/DTM1966_TerrHexN50Mosaic_v2_32N.tif';
        currentDEM=DEMname;
        info = geotiffinfo(currentDEM);
        [Z,R] = geotiffread(currentDEM);
        Z=double(Z);

		  R = refmatToMapRasterReference(info.RefMatrix, [info.Height info.Width]);
		  [DEMx,DEMy] = worldGrid(R,'gridvectors');
        [xUTM,yUTM] = meshgrid(DEMx,DEMy);

        mdO=md; %save current mesh
        %Interpolate 1966 ice surface DEM onto mesh
        md.geometry.surface=InterpFromGrid(DEMx,DEMy,Z,double(md.mesh.x),double(md.mesh.y));

		  %Find areas with negative surface in 1966 new DEM and fill those with pr-day values
        pos=find(md.geometry.surface<0);
        md.geometry.surface(pos)=mdO.geometry.surface(pos);

        %Fix areas with high thickness SW of ice cap (local adjacent small ice caps that we're not interested in)
        pos=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'./Exp/Others/Henning1966_SW_ThkSpikes.exp','node',1);
        md.geometry.thickness(find(pos))=1;
        md.geometry.thickness=md.geometry.surface-md.geometry.bed; %recalculate thickness

        % Correct negative thicknesses
        pos=find(md.geometry.thickness<1);
		  md.geometry.thickness(pos)=1;
		  md.geometry.surface(pos)=md.geometry.bed(pos)+md.geometry.thickness(pos); %recalculate surface

        if isHO,
            %Reextrude mesh vertically
            md=extrude(md,nlayers,1);

            %Make sure bed is below base
            pos1=find(md.geometry.bed>md.geometry.base);
            md.geometry.base(pos1)=md.geometry.bed(pos1);

            %Recalculate surface
            md.geometry.surface=md.geometry.base+md.geometry.thickness;

            %Set flowequation
            md = setflowequation(md,'HO','all');
        end

    end


	%Spinup from no ice?
	if isspinupnoice,
		md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
		md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface
	end

    %%%% Set SMB
    issmbspinup=1; %constant forcing
    issmbhistoric=0; %variable forcing 1960-2020
    ismonthlysmb=0;
    issmbanomaly=1;
	 istempanomaly=1; %use SMB forcing with perturbed temperature, e.g. +1C
	 issmbkill=0; %add strongly negative smb in initial ice-free regions

	 if issmbanomaly,
		 issmblinearanomaly=0; %anomaly linearly declining with surface elevation
		 issmbnorthanomaly=1; %anomaly only in .exp-defined area in the north
		 issmboldeanomaly=0; %anomaly for outlets feeding the Oldedalen catchment 
		 issmbsouthanomaly=0; 
	 end

	%Get SMB from Kamilla
	disp('   Loading SMB rates');

	if ismonthlysmb,
		ncdata	='./Data/SMB/monthly_ref_smb_JOB_100m_1960_2020_UTM32.nc';
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		%Interp SMB onto model mesh
		nmonths=length(smbdata); %length = 732 (monthly data: 1960-2020)
		nyears=nmonths/12; %calculate number of years = 61
		smb=NaN(md.mesh.numberofvertices,nmonths);

		%Fill in SMB data for each year or month
		for i=1:nmonths
			smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
		end

	%Insert SMB as model forcing
		%% monthly %%
		%setup arrays with model forcing times, with monthly forcing for each year
		forcingtimes=[];
		for i=1:nyears
			for j=1:12
				forcingtime=i-1+j/12;
				forcingtimes=[forcingtimes forcingtime];
			end
		end
		smb=smb/1000; %convert from mm to m w.e.
		smb=smb*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
		md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,nmonths);
		md.smb.mass_balance(1:end-1,:)=smb;
		md.smb.mass_balance(end,:)=forcingtimes;
	else
% 		%% yearly %%
% 		startyear=years(1)-1;
% 		modelyears=years-startyear;
% 		smb=smb*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
% 		md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,nyears);
% 		md.smb.mass_balance(1:end-1,:)=smb;
% 		md.smb.mass_balance(end,:)=modelyears;
	end
	clear smbdata years startyear modelyears

        if ismonthlysmb,
            %Extract mean of modelled SMB 1961-1970 (startyear of SMB model is 1960-01)
            % 			startyear=2010;
            % 			years=2011:2020;
            startyear=1960;
            years=1960:1969;
%             startyear=1961;
%             years=1971:1980;
				disp(['Use SMB from ' num2str(years(1)) '-' num2str(years(end))])
            nyears=length(years);
            years=years-startyear;

            %Extract smb for all months for all years
            %Calculate yearly mean for each year, for each vertex
            smb_sum_yearly=NaN(md.mesh.numberofvertices,nyears);
            %Go through all years
            year=years(1);
            for i=1:nyears
                %Sum mass balance for each year
                for j=1:md.mesh.numberofvertices
						 smb_sum_yearly(j,i)= sum(md.smb.mass_balance(j,year*12+1:year*12+12));
                end
                year=years(1)+i;
            end

            %Take mean of all years, for each vertex
            smb_mean=NaN(md.mesh.numberofvertices,1);
            for i=1:md.mesh.numberofvertices
                smb_mean(i,1)=mean(smb_sum_yearly(i,:));
            end

            md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,1); %reset SMB
            md.smb.mass_balance=smb_mean; %only need one smb since it's the mean

        else
%             %Extract mean of modelled SMB 1961-1970 (startyear of SMB model output is 1957)
%             startyear=1957;
%             years=1961:1970;
%             years=years-startyear;
%             smb=squeeze(md.smb.mass_balance(1:end-1,years));
%             model_years=[1 nyears];
% 
%             md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,length(model_years)); %reset SMB
%             %Take mean of 1961-1970 for each vertex
%             for i=1:md.mesh.numberofvertices
%                 smb_mean=mean(smb(i,:));
%                 md.smb.mass_balance(i,:)=smb_mean;
%             end
% 
%             %Fill in forcing times
%             md.smb.mass_balance(end,:)=model_years;
% 
        end
	if issmbhistoric,
        %To be implemented
	 elseif istempanomaly,
		 disp('Running with SMB with an anomaly')
		ncdata	='./Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_100m_1960_1989_UTM32_mean_annual.nc';
		
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));
		md.smb.mass_balance=NaN(md.mesh.numberofvertices,1); %reset SMB
		md.smb.mass_balance = InterpFromGrid(xdata,ydata,smbdata',md.mesh.x,md.mesh.y); %works

		md.smb.mass_balance=md.smb.mass_balance/1000; %convert from mm to m w.e.
		md.smb.mass_balance=md.smb.mass_balance*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

    else
        %do nothing
    end

	if issmbanomaly,
		if issmblinearanomaly, %elevation-dependent SMB anomaly
			disp('linear anomaly active - careful!!!!!')
			smb_anomaly_upper=-0.5; %SMB anomaly at divide (m w.e.)
			smb_anomaly_sealevel=-1.0; %SMB anomaly at sea level (m w.e.)
			elevation_upper=max(md.geometry.surface); %max surface elevation
			elevation_sealevel=0;
			%Calculate elevation-dependency (slope k in linear function y = kx+m)
			smb_k=(smb_anomaly_sealevel-smb_anomaly_upper)/(elevation_sealevel-elevation_upper);
			%Calculate spatially variable smb anomaly 
			smb_anomaly=zeros(md.mesh.numberofvertices,1);
			smb_anomaly=smb_k.*md.geometry.surface + smb_anomaly_sealevel;

			%test plot
	%        scatter(smb_anomaly,md.results.TransientSolution(end).Surface);

			%Calculate altered SMB as requested
			if length(md.smb.mass_balance) > length(smb_anomaly) %if we're using a model with an existing time-series of anomalies, don't use last row in matrix computation (forcing times)
				smb=md.smb.mass_balance(1:end-1)+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
			else
				smb=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
			end

			smb = [smb];

			%Set SMB
			md.smb.mass_balance = [smb;1];
		elseif issmbnorthanomaly,
			
			disp('--- Adding SMB anomaly in the North ---')
			%Get smb without anomaly
			smb=md.smb.mass_balance;
			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostNorthSmall.exp',1);
			smb_anomaly=0.75;
			smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);

			smb = [smb];

			%Set SMB
			md.smb.mass_balance = [smb;1];

		else %spatially uniform anomaly
			smb_anomaly=1; % m w.e.
			%Calculate altered SMB as requested
			smb=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)

			smb = [smb];

			%Set SMB
			md.smb.mass_balance = [smb;1];
		end
	end

	if issmboldeanomaly,
		disp('--- Adding SMB anomaly in Oldedalen ---')
		%Get smb without anomaly
		smb=md.smb.mass_balance;
      in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OldedalenBasins.exp',1);
		smb_anomaly=-1;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);

		smb = [smb];

		%Set SMB
		md.smb.mass_balance = [smb];
	end

	if issmbsouthanomaly,
		disp('--- Adding SMB anomaly in the south ---')
		%Get smb without anomaly
		smb=md.smb.mass_balance;
      in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostSouth.exp',1);
		smb_anomaly=-0.5;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);

		smb = [smb];

		%Set SMB
		md.smb.mass_balance = [smb];
	end


	 if issmbkill,
		 smb_anomaly=-10;
		 smb_anomaly=smb_anomaly/md.materials.rho_ice*md.materials.rho_freshwater; %convert to ice eq.
		 pos=find(md.geometry.thickness<10);
		 md.smb.mass_balance(pos,1)=smb_anomaly;
	 end
    %%%% SMB END

    %Front and GL options
    md.transient.ismovingfront=0;
    md.transient.isgroundingline=0;
    if md.transient.isgroundingline,
        md.groundingline.migration='SubelementMigration';
        md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';
        md.groundingline.friction_interpolation='SubelementFriction1';
    end

    %Calving options
    if md.transient.ismovingfront==1
        md.calving=calvingvonmises(); %activate von mises calving law

        %Define calving rate and melt rate (only effective if ismovingfront==1)
        md.frontalforcings.meltingrate=30*ones(md.mesh.numberofvertices,1);

        %Levelset options
        md.levelset.kill_icebergs=1;
    end

    %Set initialization vel's
    md.initialization.vx=md.results.StressbalanceSolution.Vx;
    md.initialization.vy=md.results.StressbalanceSolution.Vy;
    md.initialization.vel=sqrt(md.initialization.vx.^2+md.initialization.vy.^2);

    %Dont use damage model
    md.damage.D=zeros(md.mesh.numberofvertices,1);
    md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

    %Dont use thermal model
    md.transient.isthermal=0;
    md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);

    %Additional options
    md.inversion.iscontrol=0;

    %Set transient options
    if istimesteppingadaptive,
        md.timestepping=timesteppingadaptive();
		  md.timstepping.time_step_min=0.005; %            -- minimum length of time step [yr]
        md.timestepping.time_step_max=0.2;  %  -- maximum length of time step [yr]
        md.timestepping.cfl_coefficient=0.5; %- coefficient applied to cfl condition, default = 0.5
        if isHO,
%             md.settings.output_frequency=1/time_step_min; %yearly output if tstep = 0.01, need to adjust acc. to model
            md.settings.output_frequency=20; %5 = yearly output if tstep = 0.04, need to adjust acc. to model
        end
    else
        md.timestepping=timestepping();
        md.timestepping.time_step=0.05;
        md.settings.output_frequency=1/md.timestepping.time_step; %yearly
        % 		md.settings.output_frequency=1; %every time step
		disp(['--- NB! Using tstep = ' num2str(md.timestepping.time_step) ' for spinup'])
    end

    md.timestepping.final_time=nyrs;

    md.settings.solver_residue_threshold=1e-3; %default 1e-6

    md.toolkits = toolkits();

    if md.transient.isgroundingline,
        %Specify specific non-default solver for grounding line advection
        md.toolkits  = addoptions(md.toolkits,'GLheightadvectionAnalysis',asmoptions());
        md.toolkits.DefaultAnalysis=rmfield(md.toolkits.DefaultAnalysis,field); %remove field
        md.toolkits.RecoveryAnalysis=rmfield(md.toolkits.RecoveryAnalysis,field); %remove field
    end

    %Specify solver option for PETSc version <3.12
    % 	md.toolkits=addoptions(md.toolkits,'DefaultAnalysis','pc_factor_mat_solver_type');
    md.toolkits.DefaultAnalysis.pc_factor_mat_solver_type='mumps';
    md.toolkits.RecoveryAnalysis.pc_factor_mat_solver_type='mumps';
    field='pc_factor_mat_solver_package'; %remove solver package option for PETSc >3.12 to avoid warnings

    %Faster solver?
    if isbcgsolver,
        md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
    end

    if ismallocdebug,
        %Specify debugging option
        md.toolkits.DefaultAnalysis.malloc_debug = [];
        md.toolkits.RecoveryAnalysis.malloc_debug = [];
    end

    if md.transient.ismovingfront==1
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceVolumeAboveFloatation',...
            'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
            'MaskOceanLevelset','MaskIceLevelset',...
            'FloatingAreaScaled','IceMass',...
            'GroundedArea','FloatingArea','TotalFloatingBmb',...
            'BasalforcingsFloatingiceMeltingRate',...
            'TotalCalvingFluxLevelset',... %Gt/r
            'GroundinglineMassFlux',... %Gt/yr
            'CalvingMeltingrate','TotalCalvingMeltingFluxLevelset','IcefrontMassFluxLevelset',...
            'TotalCalvingFluxLevelset','TotalGroundedBmb',...
            'Calvingratex','Calvingratey','CalvingCalvingrate'};
    elseif md.transient.isgroundingline,
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceVolumeAboveFloatation',...
            'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
            'FloatingAreaScaled','IceMass',...
            'GroundedArea','FloatingArea','TotalFloatingBmb',...
            'BasalforcingsFloatingiceMeltingRate',...
            'GroundinglineMassFlux'}; %Gt/yr
    else
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceMass'};
    end

    % 					'GroundinglineHeight',...
    % 							'StrainRatexx','StrainRateyy','StrainRatexy',...
    % 							'StrainRateparallel','StrainRateperpendicular'};

    md.verbose=verbose('solution',true,'module',false,'convergence',false);

    %Go solve
        disp(['	Solving on mimi using ' num2str(ncpus) ' cpus'])
        md.cluster=generic('name',oshostname,'np',ncpus);

        %for debugging
        % 		savemodel(org,md);

        %Solve
        md=solve(md,'Transient');

    savemodel(org,md);
    md=loadmodel([modelpath num2str(org.prefix) 'Transient']);

end%}}}
if perform(org,'TransientPerturb'),% {{{
    isrestart=0;
	isloadcustommodel=0;

    if isrestart, %start from perturbed results
    	disp('Loading TransientPerturb')
    	md = loadmodel(org,'TransientPerturb');
    else
		if isloadcustommodel,
			disp('---- NB! Loading with custom spinup model')
			modelpath='/uio/hypatia/geofag-felles/projects/jostice/issm/jost/Models/';
			%	name of model to be loaded and used as initial state (end of spinup)
			modelnamestring='Jost20240821b_';
			modelloadstring=[modelpath modelnamestring 'Transient'];

			%Load custom model
			md=loadmodel(modelloadstring);

			%Rename with new name to save run as
			md.miscellaneous.name=modelstring;

		else
			disp('Loading Transient')
			md = loadmodel(org,'Transient'); %load the Transient relax model of the current experiment
		end
	end

	if ishistoric,
		nyrs=nyrsHistoric;
	else
		nyrs=nyrsPerturb; %set simulation length number of years
	end

    %Initialize geometry and update mask from previous transient results
    md=transientrestart(md);

    %Make sure we're using correct flow equation
    if isHO,
    	md = setflowequation(md,'HO','all');
    else
    	md = setflowequation(md,'SSA','all');
    end

    %SMB settings
    ischangesmb=1;
    issmbhistoric=1;
	 ismonthlysmb=0;
    issmblinearanomaly=0; %elevation-dependent SMB anomaly
    issmbuniformanomaly=0; %spatial uniform SMB anomaly
	 issmbnorthanomaly=1; %anomaly only in .exp-defined area in the north
	 issmboldeanomaly=0; %anomaly for outlets feeding the Oldedalen catchment 
	 issmbsouthanomaly=0; %anomaly for southern ice cap
    ispresent=0; %SMB back to present-day

    %Set up forcing times
    smbForcingTimes=[1:nyrs];
    smb_nyrs=length(smbForcingTimes);

    if ischangesmb,
    	if issmbhistoric, %1960-2020 simulation
    		disp('   Loading historic SMB rates');

    		if ismonthlysmb,
    			%%% 100m model data %%%
				ncdata	= './Data/SMB/monthly_refmb_JOB_100m_1960_2020_UTM32_mean_annual_dec23.nc';
    			xdata		= double(ncread(ncdata,'X'));
    			ydata		= double(ncread(ncdata,'Y'));
    			smbdata	= double(ncread(ncdata,'mb_monthly'));

				%Interp SMB onto model mesh
					% 		nmonths=length(smbdata); %length = 732 (monthly data: 1960-2020)
				nmonths=nyrs*12;
					% 		nyrs=nmonths/12; %calculate number of years = 61
				smb=NaN(md.mesh.numberofvertices,nmonths);

				%Fill in SMB data for each month
				for i=1:nmonths
					smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
				end

            %setup arrays with model forcing times, with monthly forcing for each year
            forcingtimes=[];
            for i=1:nyrs
                for j=1:12
                    forcingtime=i-1+j/12;
                    forcingtimes=[forcingtimes forcingtime];
                end
            end

			else %yeary SMB
				ncdata	= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
    			xdata		= double(ncread(ncdata,'X'));
    			ydata		= double(ncread(ncdata,'Y'));
    			smbdata	= double(ncread(ncdata,'mb_annual'));

				nyears=length(smbdata(1,1,:))
				smb=NaN(md.mesh.numberofvertices,nyears);

				%interp SMB onto model mesh
				for i=1:nyears
					smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
				end

            forcingtimes=[];
            for i=1:nyears
					forcingtimes=[forcingtimes i];
            end

    		end


    		%Insert SMB as model forcing
            smb=smb/1000; %convert from mm to m w.e.
            smb=smb*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
            md.smb.mass_balance=NaN(md.mesh.numberofvertices+1,nyears);
            md.smb.mass_balance(1:end-1,:)=smb;
            md.smb.mass_balance(end,:)=forcingtimes;

    	elseif issmblinearanomaly, %elevation-dependent SMB anomaly
    		smb_anomaly_upper=0.2; %SMB anomaly at divide (m w.e.)
    		smb_anomaly_sealevel=1.0; %SMB anomaly at sea level (m w.e.)
    		elevation_upper=max(md.geometry.surface); %max surface elevation
    		elevation_sealevel=0;
    		%Calculate elevation-dependency (slope k in linear function y = kx+m)
    		smb_k=(smb_anomaly_sealevel-smb_anomaly_upper)/(elevation_sealevel-elevation_upper);
    		%Calculate spatially variable smb anomaly
    		smb_anomaly=zeros(md.mesh.numberofvertices,1);
    		smb_anomaly=smb_k.*md.geometry.surface + smb_anomaly_sealevel;

    		%test plot
            % 			scatter(smb_anomaly,md.results.TransientSolution(end).Surface);

    		%Calculate altered SMB as requested
    		if length(md.smb.mass_balance) > length(smb_anomaly) %if we're using a model with an existing time-series of anomalies, don't use last row in matrix computation (forcing times)
    			smb=md.smb.mass_balance(1:end-1)+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
    		else
    			smb=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
    		end

    	else
            % 		%spatially uniform anomaly
            % 		smb_anomaly=1; % m w.e.
            % 		%Calculate altered SMB as requested
            % 		smb=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
    	end

        % 	smb = [smb];
        %
        % 	%Set SMB
        % 	md.smb.mass_balance = [smb;1];
    end

    if issmbuniformanomaly,
    	disp('		Adding spatially uniform SMB anomaly')
        %spatially uniform anomaly
        % 		smb_anomaly=-1*ones(md.mesh.numberofvertices,1); %m w.e.
        smb_anomaly=-1; % m w.e.
        %Calculate altered SMB as requested
        if length(md.smb.mass_balance) > length(smb_anomaly) %if we're using a model with an existing time-series of anomalies, don't use last row in matrix computation (forcing times)
            md.smb.mass_balance(1:end-1,:)=md.smb.mass_balance(1:end-1,:)+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
        else
            md.smb.mass_balance=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)
        end
    end


	if issmboldeanomaly,
		disp('--- Adding SMB anomaly in Oldedalen ---')
		%Get smb without anomaly
		smb=md.smb.mass_balance;
      in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/OldedalenBasins.exp',1);
		smb_anomaly=-1;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);

		smb = [smb];

		%Set SMB
		md.smb.mass_balance = [smb];
	end

	if issmbsouthanomaly,
		disp('--- Adding SMB anomaly in the south ---')
		%Get smb without anomaly
		smb=md.smb.mass_balance;
      in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostSouth.exp',1);
		smb_anomaly=-0.5;
		smb(find(in))=smb(find(in))+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);

		smb = [smb];

		%Set SMB
		md.smb.mass_balance = [smb];
	end

	if issmbnorthanomaly,
		disp('--- Adding SMB anomaly in the North ---')
		%Get smb without anomaly
		smb=md.smb.mass_balance;
		in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostNorthSmall.exp',1);
		smb_anomaly=0.75;
		for i=1:nyears
			smb(find(in),i)=smb(find(in),i)+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);
		end

		smb = [smb];

		%Set SMB
		md.smb.mass_balance = [smb];
	end


    %Make sure bed is below base
    pos1=find(md.geometry.bed>md.geometry.base);
    md.geometry.base(pos1)=md.geometry.bed(pos1);

    %Recalculate surface
    md.geometry.surface=md.geometry.base+md.geometry.thickness;

    %Front and GL options
    md.transient.ismovingfront=0;
    md.transient.isgroundingline=0;
    md.groundingline.migration='SubelementMigration';
    md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';
    md.groundingline.friction_interpolation='SubelementFriction1';


    %Set transient options
    if istimesteppingadaptive,
    	md.timestepping=timesteppingadaptive();
    	md.timestepping.time_step_min=0.01; %            -- minimum length of time step [yr], 0.01 is a good first try
    	md.timestepping.time_step_max=0.2;  %  -- maximum length of time step [yr]
    	md.timestepping.cfl_coefficient=0.5; %- coefficient applied to cfl condition, default = 0.5
        % 	md.settings.output_frequency=1/time_step_min; %every ~1 yr if tsteps tend to time_step_min, need to adjust this for a specific run to get approx. yearly outputs
    	md.settings.output_frequency=1/0.5;
        % 	md.settings.output_frequency=1; %every tstep, mostly for debugging
    	disp(['Using adaptive timestepping, min = ' num2str(time_step_min) ' y, max = ' num2str(time_step_max) ' y'])
    else
    	md.timestepping=timestepping();
    	md.timestepping.time_step=0.02;
    	md.settings.output_frequency=1/md.timestepping.time_step; %yearly
        % 	md.settings.output_frequency=1; %1: every tstep; 5: every fifth tstep, etc (for debugging)
    	disp(['Setting fixed time step to ' num2str(md.timestepping.time_step) ' yrs'])

    end

    if md.settings.output_frequency==1;	disp('---- Warning: output will be saved every time step --- '); end
    md.timestepping.final_time=md.timestepping.start_time+nyrs;

    %Collapse 3d HO into 2d SSA model?
    if isHO,
    	if iscollapse,
    		disp('Collapsing 3d HO -> 2d model')
    		md=collapse(md);
    		md=setflowequation(md,'SSA','all');
    		isHO=0; %deactivate HO
    	end
    end


    if md.transient.ismovingfront==1
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceVolumeAboveFloatation',...
            'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
            'MaskOceanLevelset','MaskIceLevelset',...
            'FloatingAreaScaled','IceMass',...
            'GroundedArea','FloatingArea','TotalFloatingBmb',...
            'BasalforcingsFloatingiceMeltingRate',...
            'TotalCalvingFluxLevelset',... %Gt/r
            'GroundinglineMassFlux',... %Gt/yr
            'CalvingMeltingrate','TotalCalvingMeltingFluxLevelset','IcefrontMassFluxLevelset',...
            'TotalCalvingFluxLevelset','TotalGroundedBmb',...
            'Calvingratex','Calvingratey','CalvingCalvingrate'};
    elseif md.transient.isgroundingline,
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceVolumeAboveFloatation',...
            'IceVolumeAboveFloatationScaled','GroundedAreaScaled',...
            'FloatingAreaScaled','IceMass',...
            'GroundedArea','FloatingArea','TotalFloatingBmb',...
            'BasalforcingsFloatingiceMeltingRate',...
            'GroundinglineMassFlux'}; %Gt/yr
    else
        md.transient.requested_outputs={'TotalSmb','SmbMassBalance',...
            'IceVolume','IceMass'};
    end

    % 					'GroundinglineHeight',...
    % 							'StrainRatexx','StrainRateyy','StrainRatexy',...
    % 							'StrainRateparallel','StrainRateperpendicular'};

    md.verbose=verbose('solution',true,'module',false,'convergence',false);

    %Verbose for debugging
    % 	md.verbose = verbose('all');

    %Go solve
	  md.toolkits=toolkits();
	  md.cluster=generic('name',oshostname,'np',ncpus);

	  %Solve
	  md=solve(md,'Transient');

    savemodel(org,md);

    md=loadmodel([modelpath num2str(org.prefix) 'TransientPerturb']);

end%}}}
if perform(org,'TransientCoupledHist'),% {{{

    isrestart=0;
	 isloadcustommodel=1;

   if isrestart,
    	disp('Loading TransientPerturb')
    	md = loadmodel(org,'TransientPerturb');
	elseif isloadcustommodel,
			disp('---- NB! Loading with custom spinup model')
			modelpath='/uio/hypatia/geofag-felles/projects/jostice/issm/jost/Models/';
			%	name of model to be loaded and used as initial state (end of spinup)
			modelnamestring='Jost20240821b_';
			modelloadstring=[modelpath modelnamestring 'Transient'];

			%Load custom model
			md=loadmodel(modelloadstring);

			%Rename with new name to save run as
			md.miscellaneous.name=modelstring;
    else
    	disp('Loading Transient')
    	md = loadmodel(org,'Transient'); %load the Transient relax model of the current experiment
    end

	 nyrs=nyrsHistoric; %set simulation length number of years


	 %Time settings
	 md.timestepping.time_step=0.05;
	 disp(['--- NB! Using tstep = ' num2str(md.timestepping.time_step) ' for Coupled run'])

	%%% ISSM: interpolate simulated ice geometry from relaxation, from mesh onto 100 m regular grid, write to netcdf4
	disp(['-----------   Reading modelled surface']);
	srfMesh=md.results.TransientSolution(end).Surface; %save modelled surface

	%%Get xy coords from netcdf
	ncdata	='./Data/SMB/monthly_refmb_JOB_100m_1960_1989_UTM32_mean_annual.nc';
	xdata		= double(ncread(ncdata,'X'));
	ydata		= double(ncread(ncdata,'Y'));

	%Meshgrid
	[Xdata,Ydata]=meshgrid(xdata,ydata);
	ny=length(Xdata(:,1));
	nx=length(Ydata(1,:));

	disp(['-----------   Interpolating ice surface from to grid']);
	if isHO,
		srfMesh=project2d(md,srfMesh,1);
		srfGrid = InterpFromMeshToGrid(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,srfMesh,xdata,ydata,NaN);
	else
		srfGrid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,srfMesh,xdata,ydata,NaN);
	end
	clear srfMesh

	%% Write to netcdf
	filename=['./Coupling/Output/' md.miscellaneous.name 'issmoutput.nc'];
	if exist(filename)
		delete(filename)
	end
	
	%Create file and set format
	nccreate(filename,'X','dimensions',{'X',1,nx})
	nccreate(filename,'Y','dimensions',{'Y',1,ny})
	nccreate(filename,'surface_modelled','dimensions',{'Y','X'})

	ncwrite(filename,'X',xdata) % add data to variables
	ncwrite(filename,'Y',ydata) % add data to variables
	ncwrite(filename,'surface_modelled',srfGrid) % add data to variables

	S = ncinfo(filename); %Read variable, dimension, group definitions. This information defines the file's schema.
	S.Format = 'netcdf4';

    %Initialize geometry and update mask from previous transient results
    md=transientrestart(md);
% 
    %Solution options
    md.toolkits=toolkits();
    md.cluster=generic('name',oshostname,'np',ncpus);

    modelstring=[modelpath org.prefix]; %name of model

	%define filepath to smb here
	pathsmb= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_1km_1960_2020_UTM32.nc';
	
	isconstantsmb=0;

    %Call coupler
    md=issmcoupler(md,nyrs,iscoupledfuture,coupling_interval,ncpus,isHO,modelstring,org,pathsmb,isconstantsmb);

    %Save results
    savemodel(org,md);

	 md=loadmodel([modelpath num2str(org.prefix) 'TransientCoupledHist']);


end%}}}
if perform(org,'TransientCoupledFuture'),% {{{

	isloadcustommodel=0;

	disp('Loading historic model')
	if isloadcustommodel,
		disp(['--- NB! Loading custom model as initial conditions!!'])
		modelpath='/uio/hypatia/geofag-felles/projects/jostice/issm/jost/Models/';
		%	name of model to be loaded and used as initial state (end of historic run)
		modelnamestring='Jost20240821b_'; %historic, high-res

		if iscommitnow,
			modelnamestring='Jost20240821b_'; %future coupled, larger domain in SE
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
			isconstantsmb=1;
		elseif iscommitrcp45,
			modelnamestring='Jost20240813a_'; %future coupled - RCP4.5
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		elseif iscommitrcp85,
	 		modelnamestring='Jost20240813e_'; %future coupled - RCP8.5
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		elseif isregrownow,
			modelnamestring='Jost20240821b_'; %future coupled, larger domain in SE
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
			isconstantsmb=1;
		elseif isregrowrcp45,
			modelnamestring='Jost20240813a_'; %future coupled - RCP4.5 - high-resolution
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		elseif isregrowrcp85,
	 		modelnamestring='Jost20240813e_'; %future coupled - RCP8.5 - high-resolution
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		elseif isreversercp45,
			modelnamestring='Jost20240813a_'; %future coupled - RCP4.5
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		elseif isreversercp85,
	 		modelnamestring='Jost20240813e_'; %future coupled - RCP8.5
			modelloadstring=[modelpath modelnamestring 'TransientCoupledFuture'];
			isconstantsmb=1;
		else
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
			isconstantsmb=0;
		end
		%Load custom model
		md=loadmodel(modelloadstring);

		%Rename with new name to save future run as
		md.miscellaneous.name=modelstring;

	else
		md = loadmodel(org,'TransientPerturb');
	end

    %Set transient options
    if istimesteppingadaptive,
        md.timestepping=timesteppingadaptive();
			md.timestepping=timesteppingadaptive();
			md.timestepping.time_step_min=0.005; %            -- minimum length of time step [yr], 0.01 is a good first try
			md.timestepping.time_step_max=0.2;  %  -- maximum length of time step [yr]
			md.timestepping.cfl_coefficient=0.5; %- coefficient applied to cfl condition, default = 0.5
        if isHO,
%             md.settings.output_frequency=1/time_step_min; %yearly output if tstep = 0.01, need to adjust acc. to model
            md.settings.output_frequency=100; %100 = yearly output if tstep = 0.01, need to adjust acc. to model
        end
		disp(['--- NB! Using adaptive timestepping for future coupled'])
    else
        md.timestepping=timestepping();
        md.timestepping.time_step=0.05;
        md.settings.output_frequency=1/md.timestepping.time_step; %yearly
        % 		md.settings.output_frequency=1; %every time step
		disp(['--- NB! Using tstep = ' num2str(md.timestepping.time_step) ' for future'])
    end

	 %Time settings
	 nyrs=nyrsFuture;


    modelstring=[modelpath org.prefix]; %name of model

	%Define filepath to coarse resolution SMB here, to be used in downscaling
	% The following SMB model output can be obtained upon request to the corresponding author (Henning Akesson)
	if isparis1,
		pathsmb	='./Data/SMB/future/future_ECEARTH_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif isparis2
		pathsmb	='./Data/SMB/future/future_ECEARTH_HIRHAM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif isparis3,
		pathsmb  ='./Data/SMB/future/future_CNRM_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif isparis4,
		pathsmb	='./Data/SMB/future/future_MPI_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif iswarm1,
		pathsmb  ='./Data/SMB/future/future_ECEARTH_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif iswarm2,
		pathsmb	='./Data/SMB/future/future_ECEARTH_HIRHAM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif iswarm3,
		pathsmb	='./Data/SMB/future/future_CNRM_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif iswarm4,
		pathsmb	='./Data/SMB/future/future_MPI_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
	elseif iscommitnow,
		
		%Get last 20yrs of SMB from historic run
		ncdata	= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_1km_1960_2020_UTM32.nc';
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nyears=20; %2000-2020: 20 yrs x 12 months per yr
		nmonths=nyears*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

		%We now have an array of SMB data of montly resolution for 2000-2020. Would like to repeat this forcing for 2101-2300,
		%meaning nyrsFuture/20 times (= 200/20 = 10 times)
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);

		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureCommitNow_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;
		
	elseif iscommitrcp45,

		%% - extract last 20yrs of 2021_2100 smb .nc 1 km
      pathsmb  ='./Data/SMB/future/future_ECEARTH_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
		ncdata	= pathsmb; 
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
% 		smbdata = permute(smbdata, [2, 1, 3]); %switch the x and y dimensions

		%Created repeated smb forcing 2101-2300
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);


		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureCommit45ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;

	elseif iscommitrcp85,
		%% - extract last 20yrs of 2021_2100 smb .nc 1 km
      pathsmb  ='./Data/SMB/future/future_ECEARTH_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
		ncdata	= pathsmb; 
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
% 		smbdata = permute(smbdata, [2, 1, 3]); %switch the x and y dimensions

		%Created repeated smb forcing 2101-2300
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);


		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureCommit85ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;

	elseif isregrownow,
		%Get last 20yrs of SMB from historic run
		ncdata	= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_1km_1960_2020_UTM32.nc';
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nyears=20; %2000-2020: 20 yrs x 12 months per yr
		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

		%We now have an array of SMB data of montly resolution for 2000-2020. Would like to repeat this forcing for 2101-2300,
		%meaning nyrsFuture/20 times (= 200/20 = 10 times)
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);

		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureRegrowNow_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;

		%Collapse current model 3d -> 2d
		md=collapse(md);

		%Set ice thickness to zero
		md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
		md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

		 %Dont use damage model
		 md.damage.D=zeros(md.mesh.numberofvertices,1);
		 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

		%Reextrude mesh vertically
		md=extrude(md,nlayers,1);


	elseif isregrowrcp45,
		%% - extract last 20yrs of 2021_2100 smb .nc 1 km
      pathsmb  ='./Data/SMB/future/future_ECEARTH_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
		ncdata	= pathsmb; 
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
% 		smbdata = permute(smbdata, [2, 1, 3]); %switch the x and y dimensions

		%Created repeated smb forcing 2101-2300
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);


		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureRegrow45ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;


		%Collapse current model 3d -> 2d
		md=collapse(md);

		%Set ice thickness to zero
		md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
		md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

		 %Dont use damage model
		 md.damage.D=zeros(md.mesh.numberofvertices,1);
		 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

		%Reextrude mesh vertically
		md=extrude(md,nlayers,1);

	elseif isregrowrcp85,
		%% - extract last 20yrs of 2021_2100 smb .nc 1 km
      pathsmb  ='./Data/SMB/future/future_ECEARTH_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_1km_2021_2100_UTM32.nc';
		ncdata	= pathsmb; 
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
% 		smbdata = permute(smbdata, [2, 1, 3]); %switch the x and y dimensions

		%Created repeated smb forcing 2101-2300
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);


		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureRegrow85ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;

		%Collapse current model 3d -> 2d
		md=collapse(md);

		%Set ice thickness to zero
		md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
		md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

		 %Dont use damage model
		 md.damage.D=zeros(md.mesh.numberofvertices,1);
		 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

		%Reextrude mesh vertically
		md=extrude(md,nlayers,1);


	elseif isreversercp45,

		%Get last 20yrs of SMB from historic run
		ncdata	= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_1km_1960_2020_UTM32.nc';
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nyears=20; %2000-2020: 20 yrs x 12 months per yr
		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

		%We now have an array of SMB data of montly resolution for 2000-2020. Would like to repeat this forcing for 2101-2300,
		%meaning nyrsFuture/20 times (= 200/20 = 10 times)
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);

		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureReverse45ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;


	elseif isreversercp85,

		%Get last 20yrs of SMB from historic run
		ncdata	= './Data/SMB/final/3.1_pcorr_tcorr_corr_2000_new/monthly_refmb_JOB_1km_1960_2020_UTM32.nc';
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		smbdata	= double(ncread(ncdata,'mb_monthly'));

		nyears=20; %2000-2020: 20 yrs x 12 months per yr
		nmonths=20*12; %2081-2100: 20 yrs x 12 months per yr
		nx=length(xdata);
		ny=length(ydata);
		smbdata=squeeze(smbdata(:,:,end-nmonths+1:end));
		smbdata=smbdata*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

		%We now have an array of SMB data of montly resolution for 2000-2020. Would like to repeat this forcing for 2101-2300,
		%meaning nyrsFuture/20 times (= 200/20 = 10 times)
		smb = repmat(smbdata,[1, 1, nyrsFuture/20]);

		%% - create new .nc repeating smb for these 20yrs after each other 2101-2300 (1 km)
		%Write SMB forcing to netcdf (no variation over time, but with explicit time vector to work w/ downscale_smb_spline_future.py
		ismeancalc=0; %turn off calculation of mean SMB, now we just want to write the SMB array with length nmonths to netcdf

		%Create time vector of strings in the format YEAR-MONTH, e.g. 2021-01, 2021-02, etc
		nmonths=nyrsFuture*12;
		tstrings=strings(1,nmonths); %Create string vector
		tdata=linspace(1,nmonths,nmonths); %Create time vector

		%fill in string vector for each year and month
		startyear=2100;
		k=0;
		for i=1:nyrsFuture
			for j=1:12
				k=k+1;
				if j<10
					tstrings(1,k)=[num2str(startyear+i) '-0' num2str(j)];
				else
					tstrings(1,k)=[num2str(startyear+i) '-' num2str(j)];
				end
			end
		end

		pathsmblong=['./Coupling/Output/smbFutureReverse85ECC_1km' num2str(startyear) '-'  num2str(startyear+nyrsFuture) md.miscellaneous.name '.nc'];
		createNetcdfSMB(md,pathsmblong,smb,xdata,ydata,tdata,tstrings,isHO,ismeancalc)

		%Set path of newly created coarse-res SMB, to be read by downscaling script in issmcoupler.m further down
		pathsmb=pathsmblong;

	else
	end


	 disp(['--- NB! Using tstep = ' num2str(md.timestepping.time_step) ' for Coupled run'])

	%%% ISSM: interpolate simulated ice geometry, from mesh onto 100 m regular grid, write to netcdf4
	disp(['-----------   Reading modelled surface']);
	if isregrownow==0 && isregrowrcp45==0 && isregrowrcp85==0
		srfMesh=md.results.TransientSolution(end).Surface; %save modelled surface

		 %Initialize geometry and update mask from previous transient results
		 md=transientrestart(md);

	else %if we're starting from no ice
		srfMesh=md.geometry.surface;

		 md.timestepping.start_time=0;
		 md.timestepping.final_time=md.timestepping.start_time+nyrsFuture;
	end

	%%Get xy coords from netcdf
	ncdata	='./Data/SMB/monthly_refmb_JOB_100m_1960_1989_UTM32_mean_annual.nc';
	xdata		= double(ncread(ncdata,'X'));
	ydata		= double(ncread(ncdata,'Y'));

	%Meshgrid
	[Xdata,Ydata]=meshgrid(xdata,ydata);
	ny=length(Xdata(:,1));
	nx=length(Ydata(1,:));

	disp(['-----------   Interpolating ice surface from to grid']);
	if isHO,
		srfMesh=project2d(md,srfMesh,1);
		srfGrid = InterpFromMeshToGrid(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,srfMesh,xdata,ydata,NaN);
	else
		srfGrid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,srfMesh,xdata,ydata,NaN);
	end
	clear srfMesh Xdata Ydata

	%% Write to netcdf
	filename=['./Coupling/Output/' md.miscellaneous.name 'issmoutput.nc'];
	if exist(filename)
		delete(filename)
	end

	%Create file and set format
	nccreate(filename,'X','dimensions',{'X',1,nx})
	nccreate(filename,'Y','dimensions',{'Y',1,ny})
	nccreate(filename,'surface_modelled','dimensions',{'Y','X'})

	ncwrite(filename,'X',xdata) % add data to variables
	ncwrite(filename,'Y',ydata) % add data to variables
	ncwrite(filename,'surface_modelled',srfGrid) % add data to variables

	S = ncinfo(filename); %Read variable, dimension, group definitions. This information defines the file's schema.
	S.Format = 'netcdf4';


	 %test plot
% 	 plotmodel(md,'data',md.geometry.surface,'layer',1)

% 
    %Solution options
    md.toolkits=toolkits();
    md.cluster=generic('name',oshostname,'np',ncpus);

	isconstantsmb=0;

    %Call coupler
    md=issmcoupler(md,nyrs,iscoupledfuture,coupling_interval,ncpus,isHO,modelstring,org,pathsmb,isconstantsmb);

    %Save results
    savemodel(org,md);

	 md=loadmodel([modelpath num2str(org.prefix) 'TransientCoupledFuture']);

    %  savemodel(org,md);
end%}}}
if perform(org,'TransientFuture'),% {{{

	isloadcustommodel=0;

	disp('Loading historic model')
	if isloadcustommodel,
		modelpath='/uio/hypatia/geofag-felles/projects/jostice/issm/jost/Models/';
		%	name of model to be loaded and used as initial state
% 		modelnamestring='Jost20240523a_'; %RCP4.5, large domain
		modelnamestring='Jost20240523b_'; %RCP8.5, large domain

		if iscommitnow,
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
		elseif iscommitrcp45,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		elseif iscommitrcp85,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		elseif isregrownow,
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
		elseif isregrowrcp45,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		elseif isregrowrcp85,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		elseif isreversercp45,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		elseif isreversercp85,
			modelloadstring=[modelpath modelnamestring 'TransientFuture'];
		else
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
		end
		%Load custom model
		md=loadmodel(modelloadstring);

		%Rename with new name to save future run as
		md.miscellaneous.name=modelstring;

	else
		md = loadmodel(org,'TransientPerturb');
	end

    %Initialize geometry and update mask from previous transient results
    md=transientrestart(md);
	
	 %Time settings
    %Set transient options
    if istimesteppingadaptive,
        md.timestepping=timesteppingadaptive();
        time_step_min=0.005; %            -- minimum length of time step [yr]
        time_step_max=0.2;  %  -- maximum length of time step [yr]
        cfl_coefficient=0.5; %- coefficient applied to cfl condition, default = 0.5
        if isHO,
%             md.settings.output_frequency=1/time_step_min; %yearly output if tstep = 0.01, need to adjust acc. to model
            md.settings.output_frequency=100; %100 = yearly output if tstep = 0.01, need to adjust acc. to model
        end
    else
        md.timestepping=timestepping();
        md.timestepping.time_step=0.05;
        md.settings.output_frequency=1/md.timestepping.time_step; %yearly
        % 		md.settings.output_frequency=1; %every time step
		disp(['--- NB! Using tstep = ' num2str(md.timestepping.time_step) ' for future'])
    end

	 nyrs=nyrsFuture;
	 md.timestepping.start_time=0;
	 md.timestepping.final_time=md.timestepping.start_time+nyrs;

	%SMB settings
	issmbmonthly=0;
	issmbanomaly=0;
	isloadhistsmb=0;
	issmblinearanomaly=0;
	issmbnorthanomaly=0;
	
	%%Load SMB from climate model output
	%Load data
	disp('Loading future SMB data')
	if issmbmonthly,
		if isparis1,
			ncdata	='./Data/SMB/future/future_ECEARTH_CCLM_rcp45/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_100m_2021_2100_UTM32.nc';
		elseif iswarm1,
			ncdata	='./Data/SMB/future/future_ECEARTH_CCLM_rcp85/june_2024_downscaled_to_DTM2020/monthly_refmb_JOB_100m_2021_2100_UTM32.nc';
		else
		end
	else
		if isparis1,
			ncdata	='./Data/SMB/future/future_ECEARTH_CCLM_rcp45/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc'; %misleading file name, but this is correct
		elseif isparis2
			ncdata	='./Data/SMB/future/future_ECEARTH_HIRHAM_rcp45/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif isparis3,
			ncdata	='./Data/SMB/future/future_CNRM_CCLM_rcp45/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif isparis4,
			ncdata	='./Data/SMB/future/future_MPI_CCLM_rcp45/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif iswarm1,
			ncdata	='./Data/SMB/future/future_ECEARTH_CCLM_rcp85/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc'; %misleading file name, this is correct
		elseif iswarm2,
			ncdata	='./Data/SMB/future/future_ECEARTH_HIRHAM_rcp85/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif iswarm3,
			ncdata	='./Data/SMB/future/future_CNRM_CCLM_rcp85/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif iswarm4,
			ncdata	='./Data/SMB/future/future_MPI_CCLM_rcp85/june_2024_downscaled_to_DTM2020/annual_refmb_JOB_100m_1960_2020_UTM32.nc';
		elseif iscommitnow,
			%Get last 20yrs of SMB from historic run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2000-2020
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;
			
		elseif iscommitrcp45,
			%Get last 20yrs of SMB from 2020-2100 run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2080-2100
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		elseif iscommitrcp85,
			%Get last 20yrs of SMB from 2020-2100 run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2080-2100
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		elseif isregrownow,
			%Collapse current model 3d -> 2d
			md=collapse(md);

			%Set ice thickness to zero
			md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
			md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

			 %Dont use damage model
			 md.damage.D=zeros(md.mesh.numberofvertices,1);
			 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

			%Reextrude mesh vertically
			md=extrude(md,nlayers,1);

			%Get last 20yrs of SMB from historic run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2000-2020
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;


		elseif isregrowrcp45,
			%Collapse current model 3d -> 2d
			md=collapse(md);

			%Set ice thickness to zero
			md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
			md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

			 %Dont use damage model
			 md.damage.D=zeros(md.mesh.numberofvertices,1);
			 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

			%Reextrude mesh vertically
			md=extrude(md,nlayers,1);

			%Get last 20yrs of SMB from 2020-2100 run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2080-2100
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		elseif isregrowrcp85,
			%Collapse current model 3d -> 2d
			md=collapse(md);

			%Set ice thickness to zero
			md.geometry.thickness=ones(md.mesh.numberofvertices,1); %set thickness to 1 m everywhere
			md.geometry.surface=md.geometry.bed+md.geometry.thickness; %recalculate ice surface

			 %Dont use damage model
			 md.damage.D=zeros(md.mesh.numberofvertices,1);
			 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

			%Reextrude mesh vertically
			md=extrude(md,nlayers,1);

			%Get last 20yrs of SMB from 2020-2100 run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2080-2100
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		elseif isreversercp45,
			%Temporarily save future model
			mdF=md;

			%load historic model for extraction of SMB
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
			md=loadmodel(modelloadstring);

			%Collapse current model 3d -> 2d
			md=collapse(md);

			 %Dont use damage model
			 md.damage.D=zeros(md.mesh.numberofvertices,1);
			 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

			%Reextrude mesh vertically
			md=extrude(md,nlayers,1);

			%Get last 20yrs of SMB from historic run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2000-2020
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

			%Reinsert future model as current model
			md=mdF;
			clear mdF

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		elseif isreversercp85,
			%Temporarily save future model
			mdF=md;

			%load historic model for extraction of SMB
			modelloadstring=[modelpath modelnamestring 'TransientPerturb'];
			md=loadmodel(modelloadstring);

			%Collapse current model 3d -> 2d
			md=collapse(md);

			 %Dont use damage model
			 md.damage.D=zeros(md.mesh.numberofvertices,1);
			 md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

			%Reextrude mesh vertically
			md=extrude(md,nlayers,1);

			%Get last 20yrs of SMB from historic run
			smb=md.smb.mass_balance(1:end-1,end-20:end);

			%Calculate mean SMB 2000-2020
			smbMean=NaN(md.mesh.numberofvertices,1);

			%loop over all vertices, calculate mean over all years for that vertex
			for i=1:md.mesh.numberofvertices
				smbMean(i,1)=mean(smb(i,:));
			end

			%Reinsert future model as current model
			md=mdF;
			clear mdF

    		%Insert SMB as model forcing (no variation over time)
			md.smb.mass_balance=smbMean;

		else

		end
	end


	if isparis1==1 || iswarm1==1 || isparis2==1 || iswarm2==1 || isparis3==1 || iswarm3==1 || isparis4==1 || iswarm4==1
		xdata		= double(ncread(ncdata,'X'));
		ydata		= double(ncread(ncdata,'Y'));
		if issmbmonthly,
			smbdata	= double(ncread(ncdata,'mb_monthly'));
		else
			smbdata	= double(ncread(ncdata,'mb_annual'));
		end

		if issmbmonthly,
			nmonths=length(smbdata(1,1,:));
			smb=NaN(md.mesh.numberofvertices,nmonths);

			%interp SMB onto model mesh
			for i=1:nmonths
				smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
			end
		else
			nyears=length(smbdata(1,1,:))
			smb=NaN(md.mesh.numberofvertices,nyears);

			%interp SMB onto model mesh
			for i=1:nyears
				smb(:,i) = InterpFromGrid(xdata,ydata,squeeze(smbdata(:,:,i))',md.mesh.x,md.mesh.y); %works
			end
		end

		smb=smb/1000; %convert from mm to w.e.
		smb=smb*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.

		%Set time
		forcingtimes=[];
		if issmbmonthly,
			for i=1:nyrsFuture
				for j=1:12
					forcingtime=i-1+j/12;
					forcingtimes=[forcingtimes forcingtime];
				end
			end
		else
			for i=1:nyrsFuture
				forcingtimes=[forcingtimes i];
			end
		end
	end


	
	%Manual SMB anomalies
	if issmbanomaly,
		if isloadhistsmb,
			%load mean historic smb
			ncdata	='./data/smb/202402/noprop/monthly_refmb_job_100m_1960_2020_utm32_mean_annual.nc';
			xdata		= double(ncread(ncdata,'x'));
			ydata		= double(ncread(ncdata,'y'));
			smbdata	= double(ncread(ncdata,'mb_monthly'));
			md.smb.mass_balance=nan(md.mesh.numberofvertices,1); %reset smb
			md.smb.mass_balance = interpfromgrid(xdata,ydata,smbdata',md.mesh.x,md.mesh.y); %works
			md.smb.mass_balance=md.smb.mass_balance/1000; %convert from mm to m w.e.
			md.smb.mass_balance=md.smb.mass_balance*md.materials.rho_water/md.materials.rho_ice; % convert to ice eq.
		end

		if issmblinearanomaly, %elevation-dependent SMB anomaly
			disp('linear smb anom active!!!!!')
			smb_anomaly_upper=0.2; %SMB anomaly at divide (m w.e.)
			smb_anomaly_sealevel=-4.0; %SMB anomaly at sea level (m w.e.)
			elevation_upper=max(md.geometry.surface); %max surface elevation
			elevation_sealevel=0;
			%Calculate elevation-dependency (slope k in linear function y = kx+m)
			smb_k=(smb_anomaly_sealevel-smb_anomaly_upper)/(elevation_sealevel-elevation_upper);
			%Calculate spatially variable smb anomaly 
			smb_anomaly=zeros(md.mesh.numberofvertices,1);
			smb_anomaly=smb_k.*md.geometry.surface + smb_anomaly_sealevel;

			%test plot
	%        scatter(smb_anomaly,md.results.TransientSolution(end).Surface);

			%Calculate altered SMB as requested
			smb=md.smb.mass_balance+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice); %add anomaly to current SMB field, convert anomaly to ice eq. (NB that md.smb.mass_balance is already in ice eq.)

			%Set SMB
			forcingtimes=[1];
			md.smb.mass_balance = [smb;forcingtimes];
		end

		if issmbnorthanomaly,
			disp('--- Adding SMB anomaly in the North ---')
			%Get smb without anomaly
% 			smb=md.smb.mass_balance(1:end-1,1);
			in=ContourToNodes(md.mesh.x,md.mesh.y,'./Exp/Others/JostNorthSmall.exp',1);
			smb_anomaly=0.5;
			if issmbmonthly,
				for i=1:nmonths
					smb(find(in),i)=smb(find(in),i)+smb_anomaly/12*(md.materials.rho_freshwater/md.materials.rho_ice);
				end
			else
				for i=1:nyears
					smb(find(in),i)=smb(find(in),i)+smb_anomaly*(md.materials.rho_freshwater/md.materials.rho_ice);
				end
			end

			smb = [smb];

		end
	end

	if isparis1==1 || iswarm1==1 || isparis2==1 || iswarm2==1 || isparis3==1 || iswarm3==1 || isparis4==1 || iswarm4==1
		%Set SMB
		md.smb.mass_balance=[];
		md.smb.mass_balance = [smb;forcingtimes];
	end

	 %Go solve
	 disp(['	Solving on mimi using ' num2str(ncpus) ' cpus'])
	 md.cluster=generic('name',oshostname,'np',ncpus);
	 md=solve(md,'Transient');

	 %Save and load model
	 savemodel(org,md);
	 md=loadmodel([modelpath num2str(org.prefix) 'TransientFuture']);

end%}}}

