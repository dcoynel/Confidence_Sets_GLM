function Algorithm_2_GLM(contrast4d, brainmask, design, contrastvector, thr, outputfolder)
% Inputs:
% - contrast4d (string):        path to a 4d nii file of individual 
%                               contrast estimates. gz is not supported.
% - brainmask (string):         path to a binary brain mask. Dimensions 
%                               must be identical to contrast4d.
% - design (string):            path to the group-level design matrix in a
%                               csv format, with named columns. Should 
%                               include an intercept: a column named "CONSTANT" 
%                               consisting of 1s.
% - contrastvector (vector):    contrast representing which covariate to test 
%                               (e.g. [1,0,0]). The order should be the
%                               same as in the design matrix.
% - thr (float, optional):      cohen's d threshold, default 0.5.
% - outputfolder (string, optional):    path to a folder where the ouputs
%                                       will be written. Default "output"
%                                       in the current directory.
%
% Author: David Coynel, University of Basel
% Adpated from https://github.com/AlexBowring/Confidence_Sets_Manuscript/blob/master/Biobank_simulation/Algorithm_2_scripts/Algorithm_2.m

%% initialization
if (~exist('spm', 'file'))
    error('please add SPM12 to your path')
end
if (nargin<5)
  thr = 0.5;
end
if (nargin<6)
  outputfolder = fullfile(pwd, 'output');
end 

%% Define parameters
nBoot = 5000;

%% read design matrix
designmatrix = readtable(design);
X = table2array(designmatrix);
if ~ismember('CONSTANT',designmatrix.Properties.VariableNames)
    error('the design matrix shoud include a column of 1s name "CONSTANT"')
end
% check that we have as many elements in the design matrix as 
% in the contrast vector
if (length(contrastvector)~=size(X,2))
    error('contrast vector and design matrix don''t have the same number of elements')
end
% check that the contrast vector has the proper dimension, otherwise flip
if (size(contrastvector,2)==size(X,2))
    contrastvector = contrastvector';
end

%% read data & mask
VY = spm_vol(contrast4d);
nSubj=length(VY);
p = length(contrastvector);

V = spm_vol(brainmask);
Mask = spm_read_vols(V);
dim = size(Mask);
Mask(Mask==0)=NaN;
Masklin = reshape(Mask, [dim(1)*dim(2)*dim(3) 1]);

if (~all(V.dim==VY(1).dim))
    error('brain mask and 4d data should have identical x/y/z dimensions')
end

Y = zeros([dim nSubj]);
% Load up each of the 3D images into a massive 4D image
for i=1:nSubj
  % Take care to mask each image as its read
  Y(:,:,:,i) = spm_read_vols(VY(i)).*Mask;
end
% transform and exclude voxels with no/missing data
Ylin = reshape(Y, [dim(1)*dim(2)*dim(3) nSubj]);
excludevoxels = ones([dim,1]);
excludevoxels(isnan(sum(Ylin,2))) = NaN;
excludevoxels(isnan(Masklin)) = NaN;
Ylin(isnan(excludevoxels),:) = NaN;

%% estimate Cohen's d
% normalized variance of the contrast vector
vw2 = contrastvector'*inv(X'*X)*contrastvector;
% observed betas
beta_observed = zeros(length(contrastvector), size(Ylin,1));
for voxel = 1:size(Ylin,1)
    beta_observed(:,voxel) = inv(X'*X)*X'*Ylin(voxel,:)';
end
% observed residuals
epsilon_observed = Ylin - (X*beta_observed)';
% observed variance
std_observed = zeros(size(Ylin,1),1);
for voxel = 1:size(Ylin,1)
    std_observed(voxel,1) = sqrt( (1/(nSubj-p)) * epsilon_observed(voxel,:)*epsilon_observed(voxel,:)' );
end
% observed Cohen's d
cohen_d_observed = contrastvector'*beta_observed./std_observed';


%% Cohen's d residuals
observed_resid = epsilon_observed./std_observed - .5*cohen_d_observed'.*( (epsilon_observed./std_observed).^2 -1 );

% normalize the residuals
observed_resid_std = std(observed_resid,0,2);
observed_resid_norm = observed_resid./observed_resid_std;

%% Wild t-bootstrap approximating field
observed_boundary_edges         = getBdryparams(reshape(cohen_d_observed, dim), (1 - (3/(4*nSubj - 5)))^(-1)*thr);
n_observed_boundary_edges       = size(getBdryvalues(reshape(cohen_d_observed, dim), observed_boundary_edges),1);
observed_resid_boundary_values  = zeros([n_observed_boundary_edges nSubj]);

for i=1:nSubj
    subject_resid_field                 = reshape(observed_resid_norm(:,i), [dim 1]);
    observed_resid_boundary_values(:,i) = getBdryvalues(subject_resid_field, observed_boundary_edges);
end

% Implementing the Multiplier Boostrap to obtain confidence intervals
supG_observed = zeros([nBoot 1]);
for k=1:nBoot
    % Applying the bootstrap using Rademacher variables (signflips)
    signflips                         = randi(2,[nSubj,1])*2-3;
    % Estimated boundary
    observed_boundary_bootstrap       = observed_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
    observed_boundary_resid_field     = sum(observed_boundary_bootstrap, 2)/sqrt(nSubj); 
    % Re-standardizing by bootstrap standard deviation
    observed_boot_std                 = std(observed_boundary_bootstrap, 0, 2);
    observed_boundary_resid_field     = observed_boundary_resid_field./observed_boot_std;
    supG_observed(k)                  = max(abs(observed_boundary_resid_field));
end

%% outputs: coverage probability level of 0.95
observed_cohen_d = reshape(cohen_d_observed, dim);
observed_resid_std = reshape(observed_resid_std, dim);
supGa_observed_95 = prctile(supG_observed, 95);

% Observed boundary variables
lower_contour_observed_95        = observed_cohen_d >= (1 - (3/(4*nSubj - 5)))^(-1)*thr - supGa_observed_95*sqrt(vw2)*observed_resid_std;
upper_contour_observed_95        = observed_cohen_d >= (1 - (3/(4*nSubj - 5)))^(-1)*thr + supGa_observed_95*sqrt(vw2)*observed_resid_std;

%% write parameters and nifti images
mkdir(outputfolder)
fileID = fopen(fullfile(outputfolder, 'contrastvector.txt'),'w');
for i=1:length(designmatrix.Properties.VariableNames)
    fprintf(fileID, '%s: %d\n', designmatrix.Properties.VariableNames{i}, contrastvector(i));
end
fclose(fileID);

fileID = fopen(fullfile(outputfolder, 'parameters.txt'),'w');
fprintf(fileID, 'c threshold = %s\n', num2str(thr));
fprintf(fileID, 'n bootstrap = %d\n', nBoot);
fclose(fileID);

V.fname = fullfile(outputfolder, 'observed_cohen_d.nii');
spm_write_vol(V, observed_cohen_d);

V.fname = fullfile(outputfolder, 'lower_contour_observed_95.nii');
spm_write_vol(V, lower_contour_observed_95);

V.fname = fullfile(outputfolder, 'upper_contour_observed_95.nii');
spm_write_vol(V, upper_contour_observed_95);


