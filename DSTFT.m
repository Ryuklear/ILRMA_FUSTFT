function S = DSTFT(signal,analysis_window,time_skip,number_of_padded_zeros,frequency_type,real_valued_flag,phase_flag)
%DSTFT computes a complex spectrogram (number of frequencies >= window_length)
%  
%% -- Inputs ---------------------------------------------------------
% signal:   input signal (signal_length x 1)
% analysis_window:   nonzero real-valued analysis window (window_length x 1)
% time_skip:    shift size of time frame (positive integer)
% (time_skip must be equal to or smaller than window_length)
% number_of_padded_zeros:   number of padded zeros right before FFT (nonnegative integer)
% (0 --> discrete STFT, positive integer --> FOSTFT)
% frequency_type:   1 = standard sampling (omega_k =  2*pi*k/(window_length + number_of_padded_zeros))
%                   2 = half-slide sampling (omega_k = 2*pi*(k + 0.5)/(window_length + number_of_padded_zeros))
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
% n_f:   number of frequencies when real_valued_flag = 1
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
if frequency_type == 1
    S = fft([signal(index).* analysis_window; zeros(number_of_padded_zeros, size(index,2))]); % simple STFT
    if real_valued_flag == 1
        n_f = floor((window_length + number_of_padded_zeros)/2) + 1;
        S = S(1:n_f,:); % truncate complex conjugate frequency components
    end
elseif frequency_type == 2
    S = fft([signal(index).* analysis_window; zeros(2 * number_of_padded_zeros + window_length, size(index,2))]);
    S = S(2:2:end,:); % simple STFT
    if real_valued_flag == 1
        n_f = ceil((window_length + number_of_padded_zeros)/2);
        S = S(1:n_f,:); % truncate complex conjugate frequency components
    end
end

%% phase rotation for phase-aware STFT
if phase_flag == 1
    if frequency_type == 1
        S = S.*exp(-2i*pi*(mod((0:(size(S,1)-1))'*(0:(size(S,2)-1))*time_skip, window_length + number_of_padded_zeros)/(window_length + number_of_padded_zeros))); % phase-aware STFT
    elseif frequency_type == 2
        S = S.*exp(-2i*pi*(mod(((0:(size(S,1)-1))'+ 0.5)*(0:(size(S,2)-1))*time_skip, window_length + number_of_padded_zeros)/(window_length + number_of_padded_zeros))); % phase-aware STFT
    end
end

end