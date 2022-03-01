function signal = Inv_DSTFT(S,dual_window,time_skip,signal_length,number_of_padded_zeros,frequency_type,real_valued_flag,phase_flag)
%INV_DSTFT recovers a time-domain signal from a complex spectrogram (number of frequencies >= window_length)
% 
%% -- Inputs ---------------------------------------------------------
% S:   complex spectrogram (frequencies x time frames)
% dual_window:   canonical dual window (window_length x 1)
% time_skip:    shift size of time frame (positive integer)
% (time_skip must be equal to or smaller than window_length)
% signal_length:   length of the orignal signal before STFT (positive integer)
% number_of_padded_zeros:   number of padded zeros right before FFT (nonnegative integer)
% (0 --> discrete STFT, positive integer --> FOSTFT)
% frequency_type:   1 = standard sampling (omega_k =  2*pi*k/(window_length + number_of_padded_zeros))
%                   2 = half-slide sampling (omega_k = 2*pi*(k + 0.5)/(window_length + number_of_padded_zeros))
% real_valued_flag:   0 = signal is complex-valued, 1 = signal is real-valued
% phase_flag:   0 = simple STFT, 1 = phase-aware STFT
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

%% phase rotation for phase-aware STFT spectrogram

if phase_flag == 1
    if frequency_type == 1
        S = S.*exp(2i*pi*(mod((0:(size(S,1)-1))'*(0:(size(S,2)-1))*time_skip, window_length + number_of_padded_zeros)/(window_length + number_of_padded_zeros)));
    elseif frequency_type == 2
        S = S.*exp(2i*pi*(mod(((0:(size(S,1)-1))'+ 0.5)*(0:(size(S,2)-1))*time_skip, window_length + number_of_padded_zeros)/(window_length + number_of_padded_zeros)));
    end
end

%% applying the dual window after the inverse Fourier transform
if real_valued_flag == 1
    if frequency_type == 1
        if rem(window_length + number_of_padded_zeros,2) == 0
%             signal = real(ifft([S;conj(flipud(S(2:(end-1),:)))]));
            signal = ifft([S;zeros(size(S)-[2,0])],'symmetric'); % this code is faster than the upper code
        else
%             signal = real(ifft([S;conj(flipud(S(2:end,:)))]));
            signal = ifft([S;zeros(size(S)-[1,0])],'symmetric'); % this code is faster than the upper code
        end
    elseif frequency_type == 2
        S_including_zeros = zeros(2*size(S,1),size(S,2));
        S_including_zeros(2:2:end,:) = 2*S;
        if rem(window_length + number_of_padded_zeros,2) == 0
            signal = ifft([S_including_zeros;zeros(size(S_including_zeros))],'symmetric');
        else
            signal = ifft([S_including_zeros;zeros(size(S_including_zeros)-[2,0])],'symmetric');
        end
    end
else
    if frequency_type == 1
        signal = ifft(S);
    elseif frequency_type == 2
        S_including_zeros = zeros(2*size(S,1),size(S,2));
        S_including_zeros(2:2:end,:) = 2*S;
        signal = ifft(S_including_zeros);
    end
end
signal = signal(1:window_length,:).* dual_window; % remove the padded zeros and apply the dual_window

%% stack each frame and extract time-domain signal

index = (1:window_length)' + (0:time_skip:(total_length-window_length));
frame_index = repmat(1:size(S,2),window_length,1);
signal = full(sum(sparse(index(:),frame_index(:),signal(:)),2));
signal = signal((window_length - time_skip + 1):(window_length - time_skip + signal_length)); % remove the padded zeros

end