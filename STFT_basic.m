function S = STFT_basic(signal,analysis_window,time_skip,flag)
%STFT_BASIC computes a complex spectrogram (the number of frequencies = the length of window)
%  
%% -- Inputs ---------------------------------------------------------
% signal:   input (real-valued) signal (signal_length x 1)
% analysis_window:   (positive) analysis window (window_length x 1)
% time_skip:    shift size of time frame (1 x 1) (time_skip must be equal to or smaller than window_length)
% flag:   1 = standard STFT, 2 = theoretical STFT
% (difference between standard and theoretical STFTs is phase spectrogram)
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
% n_f:   number of frequencies
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
S = fft(signal(index).* analysis_window);% / sqrt(window_length); % / sqrt(window_length) is added for the unitarity of the Fourier transform
n_f = floor(window_length/2) + 1;

%% truncate negative frequency components

if flag == 1
    S = S(1:n_f,:); % standard STFT
elseif flag == 2
    S = S(1:n_f,:).*exp(-2i*pi*(mod((0:(n_f-1))'*(0:(size(S,2)-1))*time_skip, window_length)/window_length)); % theoretical STFT
end

end