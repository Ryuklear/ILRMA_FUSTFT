function S = FUSTFT(signal,analysis_window,time_skip,frequency_type,real_valued_flag,phase_flag)
%FUSTFT computes a complex spectrogram (number of frequencies = (1/2) * window_length)
%  
%% -- Inputs ---------------------------------------------------------
% signal:   input signal (signal_length x 1)
% analysis_window:   nonzero real-valued analysis window (window_length x 1)
% (window_length must be a multiple of 4)
% time_skip:    shift size of time frame (positive integer)
% (time_skip must be equal to or smaller than (window_length/2))
% frequency_type:   1 = extract even frequency indices (omega_k =  2*pi*2*k/window_length)
%                   2 = extract odd frequency indices (omega_k = 2*pi*(2*k + 1)/window_length)
%                   3 = extract even and odd frequency indices alternately
% real_valued_flag:   0 = signal is complex-valued, 1 = signal is real-valued
% phase_flag:   0 = simple STFT, 1 = phase-aware STFT
% (difference between simple and phase-aware STFTs is phase spectrogram)
%
%% -- Temporary Variables --------------------------------------------
% signal_length:   length of the orignal input signal (1 x 1)
% window_length:   (support) length of the analysis window (1 x 1)
% 
% To make the (pseudo)inverse computation easy, (window_length - time_skip)
% or more zeros are needed to be padded before and after the input signal
% As a result, the number of time frames is
%        ceil((signal_length + window_length - time_skip)/time_skip)
%
% total_length:   length of the signal after zero padding (1 x 1)
% n_f:   number of frequencies befor undersampling when real_valued_flag = 1
%
%% -- Output ---------------------------------------------------------
% S:   complex spectrogram (frequencies x time frames)
%
%% start program

signal_length = length(signal);
window_length = length(analysis_window);
total_length = ceil((signal_length + window_length - time_skip)/time_skip  - 1) * time_skip + window_length;

%% zero padding

signal = [zeros(window_length - time_skip,1); signal; zeros(total_length - signal_length - (window_length - time_skip),1)];

%% compute complex spectrogram

index = (1:window_length)' + (0:time_skip:(total_length-window_length));
S = fft(signal(index).* analysis_window); % discrete STFT
if real_valued_flag == 1
    n_f = window_length/2 + 1;
    S = S(1:n_f,:); % truncate complex conjugate frequency components
end

%% phase rotation for phase-aware STFT

if phase_flag == 1
    S = S.*exp(-2i*pi*(mod((0:(size(S,1)-1))'*(0:(size(S,2)-1))*time_skip, window_length)/window_length)); % phase-aware STFT
end

%% undersampling of frequency components

if frequency_type == 1
    S = S(1:2:end,:); % Type-I FUSTFT
elseif frequency_type == 2
    S = S(2:2:end,:); % Type-II FUSTFT
elseif frequency_type == 3
    S(2:2:end,1:2:end) = NaN;
    S(1:2:end,2:2:end) = NaN; % Type-III FUSTFT
end

end