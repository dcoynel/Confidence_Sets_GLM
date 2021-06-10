function Y = create_signal(dim, signal_type, param, smo, trnind)

% Generates Signal
% 	Input:
%		dim:				Dimension.
%		signal_type:		Shape of the signal. 2D options are 'ramp' or 'circle'. 3D options are 'sphere'.
%		param:				Parameters of the signal. 
%							If 'signal_type' is 'ramp', then param is either a number X, in which case the signal
%							is a ramp from 0 to X, or a vector [X,Y], in which case the signal a ramp from X to Y in the x direction.
%							If 'signal_type' is 'circle', then param is a vector [M, R], where M is the magnitude of the signal and R is the radius.
%							If 'signal_type' is 'sphere', then param is a vector [M, R], where M is the magnitude of the signal and R is the radius.
%                           If 'signal_type' is 'multi', then param is a value M, which is the magnitude of the three spherical signals created. 
%       smo:                A value X giving the smoothing FWHM for isotropic smoothing
%       trnind: 			A cell array for truncating the signal back to original dimensions
%
%	Output
%		Y:					Array of size dim containing the signal
%

if (nargin<4)
  smo  = 0; 
end


switch(signal_type)
	case 'ramp'
		if size(param,2) == 1
			Y = repmat(linspace(0, param), dim(1), 1);
		elseif size(param,2) == 2
			Y = repmat(linspace(param(1), param(2)), dim(1), 1);
		else
			error('param must be a scalar X or vector [X,Y]')
		end
	case 'circle'
        % Create smooth circular signal and truncates, therefore wdim
        % must be entered as dim
		Sig = CircularSignal(dim, param(2), param(1), 0);
        [Sigs, ss]  = spm_conv(Sig, smo, smo);
        % Truncate to avoid edge effects
        tSigs       = Sigs(trnind{1}, trnind{2});
        maxtSigs    = max(tSigs(:));
        Y           = (param(1)/maxtSigs)*tSigs;
    case 'sphere'
        % Create smooth spherical signal and truncates, therefore wdim
        % must be entered as dim
		Sig = SpheroidSignal(dim, param(2), param(1), 0);
        % Smoothing the signal
        Sigs = zeros(dim);
        ss   = spm_smooth(Sig,Sigs,smo*ones(1,3));
        % Truncate to avoid edge effects
        tSigs          = Sigs(trnind{1}, trnind{2}, trnind{3});
        maxtSigs       = max(tSigs(:));
        Y              = (param(1)/maxtSigs)*tSigs;    
    case 'multi'
        % Create smooth spherical signal and truncates, therefore wdim
        % must be entered as dim
        Sig = imtranslate(SpheroidSignal(dim, 20, param, 0),[-10,-25]) + imtranslate(SpheroidSignal(dim, 15, param, 0),[20, -10]) + imtranslate(SpheroidSignal(dim, 10, param, 0),[25, 13]) + ...
              imtranslate(SpheroidSignal(dim, 18, param, 0),[-15,25]);
        % Smoothing the signal
        Sigs = zeros(dim);
        ss   = spm_smooth(Sig,Sigs,smo*ones(1,3));
        % Truncate to avoid edge effects
        tSigs              = Sigs(trnind{1}, trnind{2}, trnind{3});
        tSigs(tSigs>=param)= param;
        maxtSigs           = max(tSigs(:));
        Y                  = (param/maxtSigs)*tSigs;   
end

