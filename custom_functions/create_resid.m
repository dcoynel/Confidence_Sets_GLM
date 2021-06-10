function Y = create_resid(observed_data, observed_mean, observed_std, resid_type)

% Generates Signal
% 	Input:
%		observed_data:		Array of all subjects observed data.
%		observed_mean:		Array of the mean of all subjects data.
%		observed_std:		Array of the standard devations of all subjects data.
%		resid_type:         A value to specify the type of residual -
%                           0 (default) = standardized residuals (Y - mu^)/std^
%                           1 = unstandardized residuals (Y - mu^)
%                           2 = cohens d transformed residuals
%
%	Output
%		Y:					Vector of standardized residuals, of size (prod(dim), nSubj). 
%

dim = size(observed_mean);
temp = size(observed_data);
nSubj = temp(length(temp));

resid = bsxfun(@minus, observed_data, observed_mean);

% Reshaping to save memory
resid = reshape(resid, [prod(dim) nSubj]);
observed_mean = reshape(observed_mean, [prod(dim) 1]);
observed_std = reshape(observed_std, [prod(dim) 1]);


% Cohens d case
if resid_type == 2
	cohen_d = observed_mean./observed_std;
	Y = spdiags(1./observed_std, 0, prod(dim), prod(dim))*resid - cohen_d/2.*((spdiags(1./observed_std, 0, prod(dim), prod(dim))*resid).^2-1);
elseif resid_type == 1
    Y = resid;
else    
	Y = spdiags(1./observed_std, 0, prod(dim), prod(dim))*resid;
end 