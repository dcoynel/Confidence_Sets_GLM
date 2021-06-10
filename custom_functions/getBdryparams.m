function [bdry_params] = getBdryparams(field, c)

% Finds edges and weights needed to linearly interpolate a discretized field 
% to find the locations where the true continuous field = c
%
% Input:
% field:     random field over a domain in R^2, it is an (2+1)-dimensional array,
%            where the last dimension enumerates realisations
% c:         threshold value for excursion
%
%Output:
% bdry_params is a struct containing the edge locations either side of the
% true continuous boundary where the field = c, as well as the weights that
% the adjA_cent voxels need to be multiplied by to give the precise location where the 
% field = c assuming that the gradient of the signal is linear between adjA_cent
% voxels. 
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Alex Bowring (alex.bowring@bdi.ox.A_c.uk)
% Last changes: 10/25/2018
%__________________________________________________________________________


A_c = (field >= c);
dim = size(A_c);
D   = length(dim);

switch D
    case 2
        %%%%%%%%%%%%%%%% Case 2D random field with 4-connectivity %%%%%%%%%%%%%%%%%
        %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the second component might be called vertical)
        horz = A_c(:,2:end) | A_c(:,1:end-1);
        %%% Compute the left shifted horizontal edges
        lshift            = A_c; % initialize
        lshift(:,1:end-1) = horz;
        lshift            = lshift & ~A_c;
        % Signal locations to be used with lshift edges
        lshift_sig        = lshift(:,[dim(2) 1:dim(2)-1]);
        %%% Compute the right shifted horizontal edges
        rshift          = A_c; % initialize
        rshift(:,2:end) = horz;
        rshift          = rshift & ~A_c;
        % Signal locations to be used with rshift edges
        rshift_sig      = rshift(:,[2:dim(2) 1]);
        %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the first component might be called horizontal)
        vert = A_c(1:end-1,:) | A_c(2:end,:);
        %%% Compute the right shifted horizontal edges
        ushift = A_c;
        ushift(1:end-1,:) = vert;
        ushift = ushift & ~A_c;
        % Signal locations to be used with ushift edges
        ushift_sig      = ushift([dim(1) 1:dim(1)-1],:);
        %%% Compute the down shifted vertical edges
        dshift = A_c;
        dshift(2:end,:)   = vert;
        dshift = dshift & ~A_c;
        % Signal locations to be used with dshift edges
        dshift_sig      = dshift([2:dim(1) 1],:);
        
        % Computing the weights for the weighted linear boundary of the
        % field
        lshift_w1 = abs(field(lshift_sig) - c)./abs(field(lshift) - field(lshift_sig));
        lshift_w2 = abs(field(lshift) - c)./abs(field(lshift) - field(lshift_sig));

        rshift_w1 = abs(field(rshift_sig) - c)./abs(field(rshift) - field(rshift_sig));
        rshift_w2 = abs(field(rshift) - c)./abs(field(rshift) - field(rshift_sig));

        ushift_w1 = abs(field(ushift_sig) - c)./abs(field(ushift) - field(ushift_sig));
        ushift_w2 = abs(field(ushift) - c)./abs(field(ushift) - field(ushift_sig));

        dshift_w1 = abs(field(dshift_sig) - c)./abs(field(dshift) - field(dshift_sig));
        dshift_w2 = abs(field(dshift) - c)./abs(field(dshift) - field(dshift_sig));
        
        % Compute the length of the boundary
        len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1);
                
        % Storing parameters in a structure
        bdry_params = struct('length', len, ...
                             'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                             'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                             'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                             'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2));
          
    case 3

        %%%%%%%%%%%%%%%% Case 3D random field with 6-connectivity %%%%%%%%%%%%%%%%%
        %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the second component might be called vertical)
        horz = A_c(:,2:end,:) | A_c(:,1:end-1,:);
        %%% Compute the left shifted horizontal edges
        lshift              = A_c; % initialize
        lshift(:,1:end-1,:) = horz;
        lshift              = lshift & ~A_c;
        % Signal locations to be used with lshift edges
        lshift_sig          = lshift(:,[dim(2) 1:dim(2)-1],:);
        %%% Compute the right shifted horizontal edges
        rshift            = A_c; % initialize
        rshift(:,2:end,:) = horz;
        rshift            = rshift & ~A_c;
        % Signal locations to be used with rshift edges
        rshift_sig        = rshift(:,[2:dim(2) 1],:);
        %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the first component might be called horizontal)
        vert = A_c(1:end-1,:,:) | A_c(2:end,:,:);
        %%% Compute the up shifted vertical edges
        ushift = A_c;
        ushift(1:end-1,:,:) = vert;
        ushift = ushift & ~A_c;
        % Signal locations to be used with ushift edges
        ushift_sig = ushift([dim(1) 1:dim(1)-1],:,:);
        %%% Compute the down shifted vertical edges
        dshift = A_c;
        dshift(2:end,:,:)   = vert;
        dshift = dshift & ~A_c;
        % Signal locations to be used with dshift edges
        dshift_sig = dshift([2:dim(1) 1],:,:);
        %%%%%% depth edges
        depth = A_c(:,:,1:end-1) | A_c(:,:,2:end);
        %%% Compute the bA_ck shifted depth edges
        bshift = A_c;
        bshift(:,:,1:end-1) = depth;
        bshift = bshift & ~A_c;
        % Signal locations to be used with bshift edges
        bshift_sig = bshift(:,:,[dim(3) 1:dim(3)-1]);
        %%% Compute the front shifted depth edges
        fshift = A_c;
        fshift(:,:,2:end)   = depth;
        fshift = fshift & ~A_c;
        % Signal locations to be used with fshift edges
        fshift_sig = fshift(:,:,[2:dim(3) 1]);

        % Computing the weights for the weighted linear boundary of the
        % field
        lshift_w1 = abs(field(lshift_sig) - c)./abs(field(lshift) - field(lshift_sig));
        lshift_w2 = abs(field(lshift) - c)./abs(field(lshift) - field(lshift_sig));

        rshift_w1 = abs(field(rshift_sig) - c)./abs(field(rshift) - field(rshift_sig));
        rshift_w2 = abs(field(rshift) - c)./abs(field(rshift) - field(rshift_sig));

        ushift_w1 = abs(field(ushift_sig) - c)./abs(field(ushift) - field(ushift_sig));
        ushift_w2 = abs(field(ushift) - c)./abs(field(ushift) - field(ushift_sig));

        dshift_w1 = abs(field(dshift_sig) - c)./abs(field(dshift) - field(dshift_sig));
        dshift_w2 = abs(field(dshift) - c)./abs(field(dshift) - field(dshift_sig));

        bshift_w1 = abs(field(bshift_sig) - c)./abs(field(bshift) - field(bshift_sig));
        bshift_w2 = abs(field(bshift) - c)./abs(field(bshift) - field(bshift_sig));

        fshift_w1 = abs(field(fshift_sig) - c)./abs(field(fshift) - field(fshift_sig));
        fshift_w2 = abs(field(fshift) - c)./abs(field(fshift) - field(fshift_sig));


        % Compute the length of the boundary
        len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1) ...
                   + length(bshift_w1) + length(fshift_w1);


        % Create structure for storing parameters
        bdry_params = struct('length', len, ...
                             'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                             'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                             'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                             'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2), ...
                             'bshift', struct('edges', bshift, 'sig_edges', bshift_sig, 'w1', bshift_w1,'w2', bshift_w2), ...
                             'fshift', struct('edges', fshift, 'sig_edges', fshift_sig, 'w1', fshift_w1,'w2', fshift_w2)); 
end

