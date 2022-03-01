function [diagonal_components,non_diagonal_components] = create_tridiagonal_systems(analysis_window,time_skip,frequency_type)
%CREATE_TRIDIAGONAL_SYSTEMS creates tridiagonal_systems for inverse FUSTFT
%
%% -- Inputs ---------------------------------------------------------
% analysis_window:   nonzero real-valued analysis window (window_length x 1)
% (window_length must be a multiple of 4)
% time_skip:    shift size of time frame (positive integer)
% (time_skip must be equal to or smaller than (window_length/2))
% frequency_type:   1 = extract even frequency indices (omega_k =  2*pi*2*k/window_length)
%                   2 = extract odd frequency indices (omega_k = 2*pi*(2*k + 1)/window_length)
%                   3 = extract even and odd frequency indices alternately
%
%% -- Temporary Variable ---------------------------------------------
% window_length:   (support) length of the analysis and dual windows (1 x 1)
% 
%% -- Outputs --------------------------------------------------------
% diagonal_components:   diagonal components of tridiagonal Toeplitz systems (time_skip x 1)
% non_diagonal_components:   non-diagonal components of tridiagonal Toeplitz systems 
%                            frequency_type = 1 or 2 --> (time_skip x 1)
%                            frequency_type = 3 --> (2*time_skip x 1)
%
%% start program

window_length = length(analysis_window);
diagonal_components = [analysis_window; zeros(time_skip * ceil(window_length/time_skip) - window_length,1)];
diagonal_components = reshape(diagonal_components,time_skip,[]);
diagonal_components = sum(abs(diagonal_components).^2,2)/2;
non_diagonal_components = [analysis_window(1:(window_length/2)).*analysis_window((window_length/2 + 1):end); zeros(time_skip * ceil(window_length/(2*time_skip)) - window_length/2,1)];
non_diagonal_components = reshape(non_diagonal_components,time_skip,[]);
if frequency_type <= 2
    non_diagonal_components = sum(non_diagonal_components,2)/2;
    if frequency_type == 2
        non_diagonal_components = - non_diagonal_components;
    end
elseif frequency_type == 3
    non_diagonal_components(:,2:2:end) = -non_diagonal_components(:,2:2:end);
    non_diagonal_components = sum(non_diagonal_components,2)/2;
    s = (-1).^(floor(((0:(2*time_skip-1))' + window_length - time_skip) / time_skip));
    s = circshift(s,rem(window_length,time_skip));
    non_diagonal_components = s .* repmat(non_diagonal_components,2,1);
end

end