function fe = feSet(fe,param,val,varargin)
% Set fascicle evaluation parameters.
%
%   fe = feSet(fe,param,val,varargin)
%
%----------
% feSet(fe,'bvecsindices');
%----------   
% Set the a second diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'diffusion signal repeat',image_vals);
%----------
% Set the the s0 (no diffusion direction measurement) for a second
% diffusion signla measuremnt (repeat).
% Image_vals is the image array with the dwi signal taken from all the
% voxels in the new dataset the fe.roi
% fe = feSet(fe,'b0signalrepeat',image_vals);
%----------
%
%
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com

% Check for input parameters
if notDefined('fe'),    error('fe structure required'); end
if notDefined('param'), error('param required'); end
if ~exist('val','var'), error('Value required'); end

% Squeeze out spaces and force lower case
param = mrvParamFormat(param);

%%
switch param
  % Book-keeping
  case 'name'
    fe.name  = val;        % This one's name
  case 'type'
    fe.type  = 'faseval'; % Always fascicle evaluation
  case 'savedir'
    fe.path.savedir  = val; % Always fascicle evaluation
  
    % Set top level structure, not just single slot
  case 'life'
    fe.life  = val;  % Structure of parameters and results from LIFE analysis
    
  case 'fgfromacpc'
    % Fiber group candidate fascicles, Connectome.
    % Everything isin img coordinates in LiFE
    fe.fg  = dtiXformFiberCoords(val, fe.life.xform.acpc2img,'img');
    % Must clear when we change the fg
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);
    
  case {'fgimg', 'fg'}
    % Fiber group candidate fascicles, Connectome.
    % Everything is in img coordinates in LiFE
    %
    fe.fg  = val;
    % Must clear when we change the fg
    fe = feSet(fe,'voxel 2 fiber node pairs',[]);

  case {'fgtensors','tensors','fgq','q'}
    % fe = feSet(fe, 'tensors', tensors);    - set the passed tensors
    fe.life.fibers.tensors = val;

  case 'roi'
    % Cell array of regions of interest where we evaluate
    % Must clear v2fnp
    fe.roi   = val;
    fe       = feSet(fe,'v2fnp',[]);
    
  case {'roifromfg','fgroi','roifg'}
    name   = sprintf('roi_%s',fe.fg.name); 
    randColor = rand(1,3);
    fe.roi = dtiNewRoi(name,randColor,fefgGet(feGet(fe,'fg img'),'unique image coords'));
    
  case 'xform'
    fe.life.xform = val;  % Transforms between coords for fg, roi, and dwi data
     
  %% Diffusion data related parameters
  case {'bvecs','diffusionbvecs'}
    % feSet(fe,'bvecs');
    fe.life.bvecs = val;
  case {'bvecsindices','diffusionimagesindicesindwivolume'}
    % feSet(fe,'bvecsindices');
    fe.life.bvecsindices = val;
  case {'bvals','diffusionbvals'}
    fe.life.bvals = val;    
  case {'diffusionsignalimage','dsi', 'diffusion_signal_img'}
    fe.life.diffusion_signal_img = val;
  case {'b0signalimage','b0img', 'diffusion_s0_im','s0image'}
    fe.life.diffusion_S0_img = val;
  case {'usedvoxels'}
    fe.life.usedVoxels = val;
  case {'modeltensor'}
    fe.life.modelTensor = val;
  case {'roivoxels','roicoords'}
    % What space?  What form for the coords?
    % Always in IMG coords in LiFE.
    fe.roi.coords = val;
  
    
    %% The LiFE model
  case 'mfiber'
    fe.life.Mfiber = val;             % Fiber portion of M matrix
  case {'measuredsignalfull', 'dsigmeasured'}      % Measured signal in ROI
    fe.life.dSig  = val;
  case 'fit'
    fe.life.fit = val;
  case 'voxfit'
    fe.life.voxfit = val;
  case 'xvalfit'
    fe.life.xvalfit = val;
  
    %% Connectome fibers information.
  case {'numberofuniquefibersineachvoxel','uniquefibersnum','numberofuniquefibers','numuniquef'}
    fe.life.fibers.unique.num = val;
  case {'indextouniquefibersineachvoxel','uniquefibersindex','uniqueindex','indexesofuniquefibers','indexuniquef','uniquefibers'}
    fe.life.fibers.unique.index = val;
  case {'numberoftotalfibersineachvoxel','totalfibernmber','fibersnum','numberoffibers','numf','numfibers'}
    fe.life.fibers.total.num = val;
  case {'indexoftotalfibersineachvoxel','totalfiberindex','fibersbyvox','fibersinvox'}
    fe.life.fibers.total.index = val;
  case {'voxel2fibernodepairs','v2fnp'}
    % This has to be cleared whenever we change fg or roi
    fe.life.voxel2FNpair = val;
    % Spatial coordinate transforms for voxels and fg to coordinate frames
  case {'xformimg2acpc','img2acpc','img2acpcxform'}
    fe.life.xform.img2acpc = val;
  case {'xformacpc2img','acpc2img','acpc2imgxform'}
    fe.life.xform.acpc2img = val;
  case {'size','imgsize','volumesize','dims','dim'}
    fe.life.imagedim = val;
    
    %% Diffusion data reapeted measure parameters
  case 'dwirepeatfile'
    fe.path.dwifilerep = val;  % Diffusion weighted file used for testing results
  case {'diffusionsignalimagerepeat'}
    % Set the a second diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'diffusion signal repeat',image_vals);
    fe.rep.diffusion_signal_img = val;
  case {'s0imagerepeat'}
    % Set the the s0 (no diffusion direction measurement) for a second
    % diffusion signla measuremnt (repeat).
    %
    % Image_vals is the image array with the dwi signal taken from all the
    % voxels in the new dataset the fe.roi
    %
    % fe = feSet(fe,'b0signalrepeat',image_vals);
    fe.rep.diffusion_S0_img = val;
  case {'bvecsrepeat','diffusionbvecsrepeat'}
    % feSet(fe,'bvecsrepeat');
    fe.rep.bvecs = val;
  case {'bvecsindicesrepeat','diffusionimagesindicesindwivolumerepeat'}
    % feSet(fe,'bvecsindicesrepeat');
    fe.rep.bvecsindices = val;
  case {'bvalsrepeat','diffusionbvalsrepeat'}
    % fe = feSet(fe,'bvalsrepeat')
    fe.rep.bvals = val;
  case {'imgsizerepeat'}
    fe.rep.imagedim = val;
    
  case {'anatomyfile'}
    fe.path.anatomy = val;
  case 'dwifile'
    fe.path.dwifile = val;  % Diffusion weighted file used for testing results
  case 'dtfile'
    fe.path.dtfile = val;  % Diffusion weighted file used for testing results
  case 'dictionaryparameters'
    fe.life.M.Nphi = val{1};
    fe.life.M.Ntheta = val{2};
    fe.life.M.orient = val{3};
    fe.life.M.DictSig = val{4}; 
    fe.life.M.DictFull = val{5};
    fe.life.M.DictMean = val{6};
    fe.life.M.DictIso = val{7};
  case 'gradients'
    fe.life.fibers.grad = val{1};
    fe.life.fibers.Nelem = val{2};
  case 'indicationtensor'
    fe.life.M.Phi = val;
  case 'atomslist'   %% To revise if can be avoided
    fe.life.fibers.atoms_list = val;
  case 'voxelslist'  %% To revise if can be avoided
    fe.life.fibers.voxels_list = val;
  case 'nelem'   %% To revise if can be avoided
    fe.life.fibers.Nelem = val;
  case 'ms0'
    fe.life.M.S0 = val;
  case 'curvature'
    fe.fg.Cur = val{1};
    fe.fg.Indication = val{2};
  case 'torsion'
    fe.fg.Tor = val{1};  
    fe.fg.Indication = val{2}; 
  case {'indxbvalues','indxbvaluesrepeat'}
      fe.life.bvalues_centers = val;
      for n = 1:length(fe.life.bvalues_centers)
          bval = fe.life.bvalues_centers(n);
          fe.life.bvals_ind{n} = find(abs(fe.life.bvals - bval/1000) < 100/1000);
      end
  case 'subsampledirs'
    fe.life.diffusion_signal_img = fe.life.diffusion_signal_img(:,val);
    fe.life.bvecs = fe.life.bvecs(val,:);
    fe.life.bvals = fe.life.bvals(val);
    fe.life.bvecsindices = fe.life.bvecsindices(val);
    fe.life.M.DictSig = fe.life.M.DictSig(val,:);
    fe.life.M.DictFull = fe.life.M.DictFull(val,:);
    fe.life.M.DictMean = fe.life.M.DictMean(val,:);
    fe.life.M.DictIso = fe.life.M.DictIso(val,:);   
    
    case 'tracts_info'
        Ntracts = size(val.names,2);
        for n=1:Ntracts
            fe.life.M.tracts{n}.ind = find(val.index==n);
            fe.life.M.tracts{n}.name = val.names{n};
        end
        fe.life.M.tracts{Ntracts+1}.ind = find(val.index==0);
        fe.life.M.tracts{Ntracts+1}.name = 'not a tract';
        
    case 'dict2sph'
        aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; %convert dire from Matlab azimuth-elevetion to azimuth-inclination
        
        bvecs_sym = [fe.life.bvecs; -fe.life.bvecs];
        bvals_sym = [fe.life.bvals; fe.life.bvals];
        nDict = size(fe.life.M.Dictionaries,2);
        [nTheta,nAtoms] = size(fe.life.M.DictSig);
        Lmax = val;
        %nHarm = (Lmax+1)*(Lmax+2)/2;
        nHarm = (Lmax+1)^2;
        for n=1:nDict
            disp(['Dict=',num2str(n)]); 
            SPHcoeffs = zeros(nHarm,nAtoms);
            error_fit = zeros(1,nAtoms);
            mus = zeros(1,nAtoms);
            Dict_with_mean = zeros(size(fe.life.M.Dictionaries{n}));
            orient = fe.life.M.orient;
            Dict = fe.life.M.Dictionaries{n};
            parfor a=1:nAtoms
                %disp([num2str(a),'/',num2str(nAtoms)]);
                [Rot,~, ~] = svd(orient(:,a)); % Compute rotation to allign atoms
                d = Dict(:,a);
                d = [d; d]; % Make atom symmetric
                [ Q{a}, mu, error_fit(a)] = Diff2Tensors( d, 50, bvecs_sym, bvals_sym); % Fit Tensor to Atom and compute its mean mu
                d = d + mu;
                d = 5*d/norm(d);
                
                % save mean
                mus(a) = mu;
                
                % save atom with mean
                Dict_with_mean(:,a) = Dict(:,a) + mu;
                
                % Compute orientation of tensor
                [U,l] = eig(Q{a},'vector');
                [l,ind1] = sort(l,'descend');
                v = U(:,ind1(1));
                
                rho = abs(dot(orient(:,a),v));
%                 title_text = ['  Dict ', num2str(n), '  Atom ', num2str(a), '  Corr=', num2str(rho), '  Error fit=',num2str(error_fit(a))];
%                 disp(title_text)

                ind_nnz = find(fe.life.fit.weights~=0);
                atom_usage = nnz(fe.life.M.Phi(a,fe.life.M.ind_vox{n},ind_nnz));
                
                if (rho > 0.75) && (atom_usage > 0)
                    bvecs_rot=bvecs_sym*Rot;
                    [azi, elev, rho] = cart2sph(bvecs_rot(:,1),bvecs_rot(:,2),bvecs_rot(:,3));
                    [points] = aziElev2aziIncl([azi, elev]);
                    %SPHcoeffs(:,a) = leastSquaresSHTeven(Lmax, d, bvecs_sym*Rot, 'real', []);
                    [SPHcoeffs(:,a), ~ , e]= leastSquaresSHT(Lmax, d, points, 'real');
                else
                    SPHcoeffs(:,a) = NaN(nHarm,1);
                end

            end
            fe.life.M.SPHcoeffs{n} = SPHcoeffs;
            fe.life.M.Qtensors{n} = Q;
            fe.life.M.Qten_fit_error{n} = error_fit;
            fe.life.M.mus{n} = mus;
            fe.life.M.Dict_with_mean{n} = Dict_with_mean;
        end
        
    case 'sph2fa_md'
        bvecs_sym = [fe.life.bvecs; -fe.life.bvecs];
        bvals_sym = [fe.life.bvals; fe.life.bvals];
        nDict = size(fe.life.M.Dictionaries,2);
        [nTheta,nAtoms] = size(fe.life.M.DictSig);
        for n=1:nDict
            disp(['Dict=',num2str(n)]); 
            FAs = zeros(1,nAtoms);
            MDs = zeros(1,nAtoms);
            GFAs = zeros(1,nAtoms);
            parfor a=1:nAtoms
                s = sort(eigs(fe.life.M.Qtensors{n}{a}),'descend');
                sm = mean(s);
                % Compute FA and MD
                FAs(a) = sqrt(1.5*((s(1)-sm)^2 + (s(2)-sm)^2 + (s(3)-sm)^2)/(s(1)^2 + s(2)^2 + s(3)^2));
                MDs(a) = (s(1) + s(2) + s(3))/3; 
            end
            fe.life.M.FAs{n} = FAs;
            fe.life.M.MDs{n} = MDs;
        end
        
            
  otherwise
    error('Unknown parameter %s\n',param);
end

end

function [ Qmin, mu, error_min] = Diff2Tensors( d, nRuns, bvecs, bvals)
alpha = norm(d);
d = d/alpha;

error_min = Inf;
for r=1:nRuns
    Q0 = randQa(1);
    [ Qest, mu, error ] = Estimate_atom_tensor_demeaned( d, alpha, bvecs, bvals, 'unconstrained' ,Q0);
    if error < error_min
        error_min = error;
        Qmin = Qest;
    end
end

end

function [Qa] = randQa(sigma)
x = sigma*randn(6,1);
X = [x(1), x(4), x(6); 0, x(2), x(5); 0, 0 ,x(3)];
Qa =  X'*X;
end
