function signal = ISTFT_basic(S,dual_window,time_skip,signal_length,flag)
%ISTFT_BASIC recovers a time-domain signal from a complex spectrogram (the number of frequencies = the length of window)
% 
%% -- Inputs ---------------------------------------------------------
% S:   complex spectrogram (frequencies x time frames)
% dual_window:   (nonzero real-valued) canonical dual window (window_length x 1)
% time_skip:    shift size of time frame (1 x 1) (time_skip must be equal to or smaller than window_length)
% signal_length:   length of the orignal signal before STFT (1 x 1)
% flag:   1 = standard STFT, 2 = theoretical STFT
% (difference between standard and theoretical STFTs is phase spectrogram)
%
%% -- Temporary Variables --------------------------------------------
% window_length:   (support) length of the analysis and dual windows (1 x 1)
% total_length:   length of the signal after zero padding (1 x 1)
%
%% -- Output ---------------------------------------------------------
% signal:   time-domain (real-valued) signal (signal_length x 1)
%
%% start program

window_length = length(dual_window);
total_length = (size(S,2) - 1) * time_skip + window_length;

%% phase rotation for theoretical STFT spectrogram

if flag == 2
    S = S.*exp(2i*pi*(mod((0:(size(S,1)-1))'*(0:(size(S,2)-1))*time_skip, window_length)/window_length));
end

%% applying the dual window after the inverse Fourier transform

if rem(window_length,2) == 0
    %signal = real(ifft([S;conj(flipud(S(2:(end-1),:)))])).* dual_window * sqrt(window_length); % * sqrt(window_length) is added for the unitarity of the inverse Fourier transform
    signal = ifft([S;zeros(size(S)-[2,0])],'symmetric').* dual_window;% * sqrt(window_length); % this code is faster than the upper code
else
    %signal = real(ifft([S;conj(flipud(S(2:end,:)))])).* dual_window * sqrt(window_length); % * sqrt(window_length) is added for the unitarity of the inverse Fourier transform
    signal = ifft([S;zeros(size(S)-[1,0])],'symmetric').* dual_window;% * sqrt(window_length); % this code is faster than the upper code
end

%% stack each frame and extract time-domain signal

index = (1:window_length)' + (0:time_skip:(total_length-window_length));
frame_index = repmat(1:size(S,2),window_length,1);
signal = full(sum(sparse(index(:),frame_index(:),signal(:)),2));
signal = signal((window_length - time_skip + 1):(window_length - time_skip + signal_length)); % remove the padded zeros

end