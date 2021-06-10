function Y = create_noise(wdim, noise_type, var, smo, trnind)

% Generates either homogeneous or heterogeneous smoothing isotropic noise
% 	Input:
%		wdim:				Working dimension.
%		noise_type:			2D options are 'homo', creating homogeneous i.i.d. noise, or 'ramp', creating heterogenous noise increasing as a ramp.
%		sd:					Standard deviation of the noise. 
%							If 'noise_type' is 'homo', then param is a number X.
%							If 'noise_type' is 'ramp', then param is vector [X,Y], and the standard devation of the noise increases as a ramp
%							from X to Y in the y direction. 
%		smo:				A value X giving the smoothing FWHM for isotropic smoothing
%		trnind: 			A cell array for truncating the noise back to original dimensions
%
%	Output
%		Y:					Truncated Array of noise. 
%

noise = randn(wdim);
D     = length(size(noise));

switch D
    case 2
        % Smoothing the noise, note: this uses functions from SPM8
        [Noises, tt] = spm_conv(noise, smo, smo);
        Noises       = Noises/sqrt(tt);
        tNoises      = Noises(trnind{1}, trnind{2});
        dim          = size(tNoises);

        switch(noise_type)
            case 'homo'
                Y = var.*tNoises;
            case 'ramp'
                non_stationary_sd = repmat(linspace(var(1), var(2)), dim(2), 1)';
                Y = non_stationary_sd.*tNoises;
        end
    case 3
        % Smoothing the noise, note: this uses functions from SPM8
        Noises = zeros(wdim);
        tt     = spm_smooth(noise,Noises,smo*ones(1,3));
        Noises = Noises/sqrt(tt);      
        tNoises = Noises(trnind{1},trnind{2},trnind{3});  
        
        switch(noise_type)
            case 'homo'
                Y = var.*tNoises;
        end
end
