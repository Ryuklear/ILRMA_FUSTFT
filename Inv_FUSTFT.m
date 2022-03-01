function signal = Inv_FUSTFT(S,analysis_window,diagonal_components,non_diagonal_components,time_skip,signal_length,frequency_type,real_valued_flag,phase_flag,periodicity_flag)
%INV_IFUSTFT recovers a time-domain signal from a complex spectrogram (number of frequencies = (1/2) * window length)
% 
%% -- Inputs ---------------------------------------------------------
% S:   complex spectrogram (frequencies x time frames)
% analysis_window:   nonzero real-valued analysis window (window_length x 1)
% (window_length must be a multiple of 4)
% diagonal_components:   diagonal components of tridiagonal systems (time_skip x 1)
% non_diagonal_components:   non-diagonal components of tridiagonal systems
%                            frequency_type = 1 or 2 --> (time_skip x 1)
%                            frequency_type = 3 --> (2*time_skip x 1)
% time_skip:    shift size of time frame (positive integer)
% (time_skip must be equal to or smaller than (window_length/2))
% signal_length:   length of the orignal signal before STFT (positive integer)
% frequency_type:   1 = extract even frequency indices (omega_k =  2*pi*2*k/window_length)
%                   2 = extract odd frequency indices (omega_k = 2*pi*(2*k + 1)/window_length)
%                   3 = extract even and odd frequency indices alternately
% real_valued_flag:   0 = signal is complex-valued, 1 = signal is real-valued
% phase_flag:   0 = simple STFT, 1 = phase-aware STFT
% (difference between simple and phase-aware STFTs is phase spectrogram)
% periodicity_flag:   0 = non_periodic (pseudoinverse)
%                     1 = using the periodic assumption (NOT pseudoinverse)
%
%% -- Temporary Variables --------------------------------------------
% window_length:   (support) length of the analysis window (1 x 1)
% total_length:   length of the signal after zero padding (1 x 1)
% S_including_zeros:   complex spectrogram after zero padding ((window_length/2 + 1) x time frames)
%
%% -- Output ---------------------------------------------------------
% signal:   time-domain (real-valued) signal (signal_length x 1)
%
%% start program

window_length = length(analysis_window);
total_length = (size(S,2) - 1) * time_skip + window_length;

%% add zeros at under-sampled frequency components

if real_valued_flag == 1
    S_including_zeros = zeros(window_length/2 +1,size(S,2));
else
    S_including_zeros = zeros(window_length, size(S,2));
end

if frequency_type == 1
    S_including_zeros(1:2:end,:) = S;
elseif frequency_type == 2
    S_including_zeros(2:2:end,:) = S;
elseif frequency_type == 3
    S_including_zeros(~isnan(S)) = S(~isnan(S));
end

%% phase rotation for phase-aware STFT spectrogram

if phase_flag == 1
    S_including_zeros = S_including_zeros.*exp(2i*pi*(mod((0:(size(S_including_zeros,1)-1))'*(0:(size(S_including_zeros,2)-1))*time_skip, window_length)/window_length));
end

%% applying the analysis window after the under-sampled inverse Fourier transform

if real_valued_flag == 1
    %signal = real(ifft([S_including_zeros;conj(flipud(S_including_zeros(2:(end-1),:)))])).* analysis_window;
    signal = ifft([S_including_zeros;zeros(size(S_including_zeros)-[2,0])],'symmetric').* analysis_window; % this code is faster than the upper code
else
    signal = ifft(S_including_zeros).* analysis_window;
end

%% stack each frame and extract time-domain signal

index = (1:window_length)' + (0:time_skip:(total_length-window_length));
frame_index = repmat(1:size(S,2),window_length,1);
signal = full(sum(sparse(index(:),frame_index(:),signal(:)),2));

%% solve tridiagnoal systems
if periodicity_flag == 0 % ISTFT is the pseudoinverse
    signal = signal((window_length - time_skip + 1):(window_length - time_skip + signal_length)); % remove the padded zeros
    if (real_valued_flag == 1) && (rem(window_length/2,time_skip) == 0) % solve tridiagnoal Toeplitz systems by using DST
        equation_index = rem(0:(window_length/2-1), time_skip) + 1;
        if frequency_type == 1 || frequency_type == 2
            equation_index_2 = equation_index;
        elseif frequency_type == 3
            equation_index_2 = rem(0:(window_length/2-1), 2*time_skip) + 1;
        end
        if rem(signal_length,window_length/2) == 0
            signal = (reshape(signal,window_length/2,[]))';
            if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            end
            p = zeros(2*(size(signal,1)+1), window_length/2);
            p(2:(size(signal,1)+1),:) = signal;
            q = fft(p);
            q = -imag(q(2:(size(signal,1)+1),:)); % DST
            q = q./(diagonal_components(equation_index)' + 2 * non_diagonal_components(equation_index_2)' .* cos(pi*(1:size(signal,1))'/(size(signal,1)+1)));
            p(2:(size(signal,1)+1),:) = q;
            q = fft(p);
            signal = (-2*imag(q(2:(size(signal,1)+1),:)) / (size(signal,1)+1)); % IDST
            if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            end
            signal = signal';
            signal = signal(:);

        else
            r = rem(signal_length,window_length/2);
            signal = [signal; zeros((window_length/2) * ceil(signal_length/(window_length/2)) - signal_length,1)];
            signal = (reshape(signal,window_length/2,[]))';
            if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            end
            p = zeros(2*(size(signal,1)+1), r);
            p(2:(size(signal,1)+1),:) = signal(:,1:r);
            q = fft(p);
            q = -imag(q(2:(size(signal,1)+1),:)); % DST
            q = q./(diagonal_components(equation_index(1:r))' + 2 * non_diagonal_components(equation_index_2(1:r))' .* cos(pi*(1:size(signal,1))'/(size(signal,1)+1)));
            p(2:(size(signal,1)+1),:) = q;
            q = fft(p);
            signal(:,1:r) = -2*imag(q(2:(size(signal,1)+1),:)) / (size(signal,1)+1); % IDST

            p = zeros(2*size(signal,1), window_length/2 - r);
            p(2:size(signal,1),:) = signal(1:(end-1),(r+1):end);
            q = fft(p);
            q = -imag(q(2:size(signal,1),:)); % DST
            q = q./(diagonal_components(equation_index((r+1):end))' + 2 * non_diagonal_components(equation_index_2((r+1):end))' .* cos(pi*(1:(size(signal,1)-1))'/size(signal,1)));
            p(2:size(signal,1),:) = q;
            q = fft(p);
            signal(1:(end-1),(r+1):end) = -2*imag(q(2:size(signal,1),:)) / size(signal,1); % IDST
            if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            end
            signal = signal';
            signal = reshape(signal(1:signal_length),[],1);

        end
        
    else % solve tridiagnoal systems by using the Thomas algorithm (LU decomposition)
        equation_index = rem(rem(window_length,time_skip) + (0:(window_length/2-1)), time_skip) + 1;
        if frequency_type == 1 || frequency_type == 2
            equation_index_2 = equation_index;
        elseif frequency_type == 3
            equation_index_2 = rem(rem(window_length,time_skip) + (0:(window_length/2-1)), 2*time_skip) + 1;
        end
        if rem(signal_length,window_length/2) == 0
            signal = (reshape(signal,window_length/2,[]))';       
            u = zeros(size(signal));
            l = zeros(size(signal)-[1,0]);
            y = zeros(size(signal));
            u(1,:) = diagonal_components(equation_index)';
            y(1,:) = signal(1,:); 
            for k = 2:size(signal,1)
                equation_index = rem(equation_index + window_length/2 - 1, time_skip) + 1;
                k_1 = k - 1;
                l(k_1,:) = non_diagonal_components(equation_index_2)' ./ u(k_1,:);
                u(k,:) = diagonal_components(equation_index)' - l(k_1,:) .* non_diagonal_components(equation_index_2)';
                y(k,:) = signal(k,:) - l(k_1,:) .* y(k_1,:);
                if frequency_type == 1 || frequency_type == 2
                    equation_index_2 = equation_index;
                elseif frequency_type == 3
                    equation_index_2 = rem(equation_index_2 + window_length/2 - 1, 2*time_skip) + 1;
                end
            end

            signal(end,:) = y(end,:) ./ u(end,:);
            for k = (size(signal,1)-1):-1:1
                if frequency_type == 1 || frequency_type == 2
                    equation_index_2 = mod(equation_index_2 - window_length/2 - 1, time_skip) + 1;
                elseif frequency_type == 3
                    equation_index_2 = mod(equation_index_2 - window_length/2 - 1, 2*time_skip) + 1;
                end
                signal(k,:) = (y(k,:) - non_diagonal_components(equation_index_2)' .* signal(k+1,:)) ./ u(k,:);
            end
            signal = signal';
            signal = signal(:);
            
        else
            r = rem(signal_length,window_length/2);
            signal = [signal; zeros((window_length/2) * ceil(signal_length/(window_length/2)) - signal_length,1)];
            signal = (reshape(signal,window_length/2,[]))';
            u = zeros(size(signal));
            l = zeros(size(signal)-[1,0]);
            y = zeros(size(signal));
            u(1,:) = diagonal_components(equation_index)';
            y(1,:) = signal(1,:); 
            for k = 2:(size(signal,1)-1)
                equation_index = rem(equation_index + window_length/2 - 1, time_skip) + 1;
                k_1 = k - 1;
                l(k_1,:) = non_diagonal_components(equation_index_2)' ./ u(k_1,:);
                u(k,:) = diagonal_components(equation_index)' - l(k_1,:) .* non_diagonal_components(equation_index_2)';
                y(k,:) = signal(k,:) - l(k_1,:) .* y(k_1,:);
                if frequency_type == 1 || frequency_type == 2
                    equation_index_2 = equation_index;
                elseif frequency_type == 3
                    equation_index_2 = rem(equation_index_2 + window_length/2 - 1, 2*time_skip) + 1;
                end
            end
            k = size(signal,1);
            equation_index = rem(equation_index + window_length/2 - 1, time_skip) + 1;
            k_1 = k - 1;
            l(k_1,1:r) = non_diagonal_components(equation_index_2(1:r))' ./ u(k_1,1:r);
            u(k,1:r) = diagonal_components(equation_index(1:r))' - l(k_1,1:r) .* non_diagonal_components(equation_index_2(1:r))';
            y(k,1:r) = signal(k,1:r) - l(k_1,1:r) .* y(k_1,1:r);
            if frequency_type == 1 || frequency_type == 2
                equation_index_2 = equation_index;
            elseif frequency_type == 3
                equation_index_2 = rem(equation_index_2 + window_length/2 - 1, 2*time_skip) + 1;
            end
            
            signal(end,1:r) = y(end,1:r) ./ u(end,1:r);
            for k = (size(signal,1)-1):-1:1
                if frequency_type == 1 || frequency_type == 2
                    equation_index_2 = mod(equation_index_2 - window_length/2 - 1, time_skip) + 1;
                elseif frequency_type == 3
                    equation_index_2 = mod(equation_index_2 - window_length/2 - 1, 2*time_skip) + 1;
                end
                signal(k,:) = (y(k,:) - non_diagonal_components(equation_index_2)' .* signal(k+1,:)) ./ u(k,:);
            end
            signal = signal';
            signal = reshape(signal(1:signal_length),[],1);
            
        end
    end
    
elseif periodicity_flag == 1 % ISTFT is NOT the pseudoinverse
    total_length = total_length - window_length + time_skip;
    if frequency_type == 1 || frequency_type == 2
        while rem(total_length,window_length/2) ~= 0
            total_length = total_length + time_skip;
        end
    elseif frequency_type == 3
        while (rem(total_length,window_length/2) ~= 0) || (rem(total_length,2*time_skip) ~= 0)
            total_length = total_length + time_skip;
        end
    end
    if (real_valued_flag == 1) && (rem(window_length/2,time_skip) == 0) % solve tridiagnoal circulant systems by using FFT
        equation_index = rem(0:(window_length/2-1), time_skip) + 1;
        if frequency_type == 1 || frequency_type == 2
            equation_index_2 = equation_index;
        elseif  frequency_type == 3
            equation_index_2 = rem(0:(window_length/2-1), 2*time_skip) + 1;
        end
        signal = [signal;zeros(total_length - length(signal) + window_length - time_skip,1)];
        signal((end - window_length + time_skip + 1):end) = signal((end - window_length + time_skip + 1):end) + signal(1:(window_length - time_skip));
        signal = signal((window_length - time_skip + 1):end); 
        signal = (reshape(signal,window_length/2,[]))';
        if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
            if rem(2*total_length/window_length,4) == 0
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            else
                signal(2:2:end,:) = -1i*signal(2:2:end,:); % multiplication process by imaginary unit
            end
        end
        signal = fft(signal);
        if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0) && (rem(2*total_length/window_length,4) == 2))
            q = fft([diagonal_components(equation_index)';
            -1i*(non_diagonal_components(equation_index_2)');
            zeros(size(signal) - [3,0]);
            1i*(non_diagonal_components(equation_index_2)')]);
            signal = signal./q;
            signal = ifft(signal);
        else
            q = fft([diagonal_components(equation_index)';
            non_diagonal_components(equation_index_2)';
            zeros(size(signal) - [3,0]);
            non_diagonal_components(equation_index_2)']);
            signal = signal./q;
            signal = ifft(signal,'symmetric');
        end
        
        if ((frequency_type == 3) && (rem(window_length/2,2*time_skip) ~= 0))
            if rem(2*total_length/window_length,4) == 0
                signal(3:4:end,:) = -signal(3:4:end,:); % sign reversal process
                signal(4:4:end,:) = -signal(4:4:end,:); % sign reversal process
            else
                signal(2:2:end,:) = 1i*signal(2:2:end,:); % multiplication process by imaginary unit
                signal = real(signal); 
            end
        end
        signal = signal';
        signal = reshape(signal(1:signal_length),[],1);
        
    else % solve periodic tridiagonal systems by using the LU decomposition
        equation_index = rem(rem(window_length,time_skip) + (0:(window_length/2-1)), time_skip) + 1;
        if frequency_type == 1 || frequency_type == 2
            equation_index_2 = equation_index;
        elseif frequency_type == 3
            equation_index_2 = rem(rem(window_length,time_skip) + (0:(window_length/2-1)), 2*time_skip) + 1;
        end
        signal = [signal;zeros(total_length - length(signal) + window_length - time_skip,1)];
        signal((end - window_length + time_skip + 1):end) = signal((end - window_length + time_skip + 1):end) + signal(1:(window_length - time_skip));
        signal = signal((window_length - time_skip + 1):end); 
        signal = (reshape(signal,window_length/2,[]))';
        last_index = rem(equation_index + (size(signal,1)-1) * window_length/2 - 1, time_skip) + 1;
        if frequency_type == 1 || frequency_type == 2
            last_index_2 = last_index;
        elseif frequency_type == 3
            last_index_2 = rem(equation_index_2 + (size(signal,1)-1) * window_length/2 - 1, 2*time_skip) + 1;
        end
        u = zeros(size(signal)-[1,0]);
        l = zeros(size(signal)-[2,0]);
        ur = zeros(size(signal));
        lb = zeros(size(signal)-[1,0]);
        y = zeros(size(signal));
        u(1,:) = diagonal_components(equation_index)';
        y(1,:) = signal(1,:);
        ur(1,:) = non_diagonal_components(last_index_2)';
        lb(1,:) = non_diagonal_components(last_index_2)' ./ u(1,:);

        for k = 2:(size(signal,1)-2)
            equation_index = rem(equation_index + window_length/2 - 1, time_skip) + 1;
            k_1 = k - 1;
            l(k_1,:) = non_diagonal_components(equation_index_2)' ./ u(k_1,:);
            u(k,:) = diagonal_components(equation_index)' - l(k_1,:) .* non_diagonal_components(equation_index_2)';
            y(k,:) = signal(k,:) - l(k_1,:) .* y(k_1,:);
            ur(k,:) = - l(k_1,:) .* ur(k_1,:);
            lb(k,:) = - (lb(k_1,:) .* non_diagonal_components(equation_index_2)') ./ u(k,:);
            if frequency_type == 1 || frequency_type == 2
                equation_index_2 = equation_index;
            elseif frequency_type == 3
                equation_index_2 = rem(equation_index_2 + window_length/2 - 1, 2*time_skip) + 1;
            end
        end

        k = size(signal,1)-1;
        equation_index = rem(equation_index + window_length/2 - 1, time_skip) + 1;
        k_1 = k - 1;
        l(k_1,:) = non_diagonal_components(equation_index_2)' ./ u(k_1,:);
        u(k,:) = diagonal_components(equation_index)' - l(k_1,:) .* non_diagonal_components(equation_index_2)';
        y(k,:) = signal(k,:) - l(k_1,:) .* y(k_1,:);
        pre_equation_index_2 = equation_index_2;
        if frequency_type == 1 || frequency_type == 2
            equation_index_2 = equation_index;
        elseif frequency_type == 3
            equation_index_2 = rem(equation_index_2 + window_length/2 - 1, 2*time_skip) + 1;
        end
        ur(k,:) = non_diagonal_components(equation_index_2)' - l(k_1,:) .* ur(k_1,:);
        lb(k,:) = (non_diagonal_components(equation_index_2)'- lb(k_1,:) .* non_diagonal_components(pre_equation_index_2)') ./ u(k,:);
        y(end,:) = signal(end,:) - sum(lb(1:k,:).*y(1:k,:),1);
        ur(end,:) = diagonal_components(last_index)' - sum(lb(1:k,:).*ur(1:k,:),1);
        
        
        signal(end,:) = y(end,:) ./ ur(end,:);
        signal(k,:) = (y(k,:) -  ur(k,:) .* signal(end,:)) ./ u(k,:);
        for k = (size(signal,1)-2):-1:1
            if frequency_type == 1 || frequency_type == 2
                equation_index_2 = mod(equation_index_2 - window_length/2 - 1, time_skip) + 1;
            elseif frequency_type == 3
                equation_index_2 = mod(equation_index_2 - window_length/2 - 1, 2*time_skip) + 1;
            end
            signal(k,:) = (y(k,:) - non_diagonal_components(equation_index_2)' .* signal(k+1,:) - ur(k,:) .* signal(end,:)) ./ u(k,:);
        end
        signal = signal';
        signal = reshape(signal(1:signal_length),[],1);
        
    end
end

end