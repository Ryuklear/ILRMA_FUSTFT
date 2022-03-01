%実験１の音源
[x,~] = audioread('drum_loop.wav');
window_length = 4096*3;

%analysis_window = (0.5 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length))/sqrt(window_length); % normalized symmetric Hann window
analysis_window = exp(-4*pi*((0:(window_length-1))' - (window_length-1)/2).^2/(window_length^2))/sqrt(window_length); % 正規化打ち切り対称ガウス窓
%analysis_window = (0.42 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length) +  0.08*cos(4*pi*((0:(window_length-1))'+0.5)/window_length))/sqrt(window_length); % 正規化対称ブラックマン窓

signal_length = length(x);


%実験２の音源
[x,~] = audioread('piano_c5_mono.wav');
window_length = 512;

%analysis_window = sin(pi*((0:(window_length-1))'+0.5)/window_length)/sqrt(window_length);  % normalized symmetric sine window
%analysis_window = exp(-4*pi*((0:(window_length-1))' - (window_length-1)/2).^2/(window_length^2))/sqrt(window_length); % 正規化打ち切り対称ガウス窓
%analysis_window = (0.42 - 0.5*cos(2*pi*((0:(window_length-1))'+ 0.5)/window_length) +  0.08*cos(4*pi*((0:(window_length-1))'+0.5)/window_length))/sqrt(window_length); % 正規化対称ブラックマン窓

x = x(1:end - rem(size(x,1),window_length));
signal_length = length(x);


number_of_padded_zeros = window_length ;
real_valued_flag = 1; %0だと複素数値、１だと実数値 グラフ表示は実数値で表示した方が良い
phase_flag = 1;
time_skip = window_length/8;

STFT_X = DSTFT(x,analysis_window,time_skip,0,1,real_valued_flag,phase_flag);
STFT_X_Type2 = DSTFT(x,analysis_window,time_skip,0,2,real_valued_flag,phase_flag);
FOSTFT_X  = DSTFT(x,analysis_window,time_skip,number_of_padded_zeros,1,real_valued_flag,phase_flag);
FOSTFT_X_Type2 = DSTFT(x,analysis_window,time_skip,number_of_padded_zeros,2,real_valued_flag,phase_flag);
FUSTFT_X = FUSTFT(x,analysis_window,time_skip,1,real_valued_flag,phase_flag);
FUSTFT_X_Type2 = FUSTFT(x,analysis_window,time_skip,2,real_valued_flag,phase_flag);
FUSTFT_X_Type3 = FUSTFT(x,analysis_window,time_skip,3,real_valued_flag,phase_flag);

fprintf('size of STFT_X is %d times %d \n',size(STFT_X,1),size(STFT_X,2));
fprintf('size of STFT_X_Type2 is %d times %d \n',size(STFT_X_Type2,1),size(STFT_X_Type2,2));
fprintf('size of FOSTFT_X is %d times %d \n',size(FOSTFT_X,1),size(FOSTFT_X,2));
fprintf('size of FOSTFT_X_Type2 is %d times %d \n',size(FOSTFT_X_Type2,1),size(FOSTFT_X_Type2,2));
fprintf('size of FUSTFT_X is %d times %d \n',size(FUSTFT_X,1),size(FUSTFT_X,2));
fprintf('size of FUSTFT_X_Type2 is %d times %d \n',size(FUSTFT_X_Type2,1),size(FUSTFT_X_Type2,2));
fprintf('size of FUSTFT_X_Type3 is %d times %d \n',size(FUSTFT_X_Type3,1),size(FUSTFT_X_Type3,2));

dual_window = create_canonical_dual_window(analysis_window,time_skip);
x_hat = Inv_DSTFT(STFT_X,dual_window,time_skip,signal_length,0,1,real_valued_flag,phase_flag);
x_hat_Type2 = Inv_DSTFT(STFT_X_Type2,dual_window,time_skip,signal_length,0,2,real_valued_flag,phase_flag);
FO_x_hat = Inv_DSTFT(FOSTFT_X ,dual_window,time_skip,signal_length,number_of_padded_zeros,1,real_valued_flag,phase_flag);
FO_x_hat_Type2 = Inv_DSTFT(FOSTFT_X_Type2,dual_window,time_skip,signal_length,number_of_padded_zeros,2,real_valued_flag,phase_flag);

fprintf('reconstruction error of STFT_X is %d \n',norm(x-x_hat)/norm(x));
fprintf('reconstruction error of STFT_X_Type2 is %d \n',norm(x-x_hat_Type2)/norm(x));
fprintf('reconstruction error of FOSTFT_X is %d \n',norm(x-FO_x_hat)/norm(x));
fprintf('reconstruction error of FOSTFT_X_Type2 is %d \n',norm(x-FO_x_hat_Type2)/norm(x));

[diagonal_components,non_diagonal_components_1] = create_tridiagonal_systems(analysis_window,time_skip,1);
[~,non_diagonal_components_2] = create_tridiagonal_systems(analysis_window,time_skip,2);
[~,non_diagonal_components_3] = create_tridiagonal_systems(analysis_window,time_skip,3);

FU_x_hat = Inv_FUSTFT(FUSTFT_X,analysis_window,diagonal_components,non_diagonal_components_1,time_skip,signal_length,1,real_valued_flag,phase_flag,0);
FU_x_hat_Type2 = Inv_FUSTFT(FUSTFT_X_Type2,analysis_window,diagonal_components,non_diagonal_components_2,time_skip,signal_length,2,real_valued_flag,phase_flag,0);
FU_x_hat_Type3 = Inv_FUSTFT(FUSTFT_X_Type3,analysis_window,diagonal_components,non_diagonal_components_3,time_skip,signal_length,3,real_valued_flag,phase_flag,0);
FU_x_hat_p = Inv_FUSTFT(FUSTFT_X,analysis_window,diagonal_components,non_diagonal_components_1,time_skip,signal_length,1,real_valued_flag,phase_flag,1);
FU_x_hat_Type2_p = Inv_FUSTFT(FUSTFT_X_Type2,analysis_window,diagonal_components,non_diagonal_components_2,time_skip,signal_length,2,real_valued_flag,phase_flag,1);
FU_x_hat_Type3_p = Inv_FUSTFT(FUSTFT_X_Type3,analysis_window,diagonal_components,non_diagonal_components_3,time_skip,signal_length,3,real_valued_flag,phase_flag,1);

fprintf('reconstruction error of FUSTFT_X is %d \n',norm(x-FU_x_hat)/norm(x));
fprintf('reconstruction error of FUSTFT_X_Type2 is %d \n',norm(x-FU_x_hat_Type2)/norm(x));
fprintf('reconstruction error of FUSTFT_X_Type3 is %d \n',norm(x-FU_x_hat_Type3)/norm(x));
fprintf('reconstruction error of FUSTFT_X under periodicity is %d \n',norm(x-FU_x_hat_p)/norm(x));
fprintf('reconstruction error of FUSTFT_X_Type2 under periodicity is %d \n',norm(x-FU_x_hat_Type2_p)/norm(x));
fprintf('reconstruction error of FUSTFT_X_Type3 under periodicity is %d \n',norm(x-FU_x_hat_Type3_p)/norm(x));







%実験２
 %figure(1)
 %surf(20*log10(abs(STFT_X)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(STFT_X,2) 1 size(STFT_X,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('STFT')
 
 %figure(2)
 %surf(20*log10(abs(STFT_X_Type2)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(STFT_X_Type2,2) 1 size(STFT_X_Type2,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('Type-II STFT')
 
 %figure(11)
 %surf(20*log10(abs(FOSTFT_X)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(FOSTFT_X,2) 1 size(FOSTFT_X,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('FOSTFT')
 
 %figure(12)
 %surf(20*log10(abs(FOSTFT_X_Type2)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(FOSTFT_X_Type2,2) 1 size(FOSTFT_X_Type2,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('Type-II FOSTFT')
 
 %figure(21)
 %surf(20*log10(abs(FUSTFT_X)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(FUSTFT_X,2) 1 size(FUSTFT_X,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('FUSTFT')
 
 %figure(22)
 %surf(20*log10(abs(FUSTFT_X_Type2)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(FUSTFT_X_Type2,2) 1 size(FUSTFT_X_Type2,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('Type-II FUSTFT')
 
 %F = FUSTFT_X_Type3;
 %F(isnan(F)) = min(abs(FUSTFT_X_Type3(:)));
 %figure(23)
 %surf(20*log10(abs(F)),'EdgeColor','none');
 %view([0 90])
 %xlabel('Time Frame')
 %ylabel('Frequency')
 %axis([1 size(F,2) 1 size(F,1)])
 %colormap jet
 %colorbar
 %caxis([ -100    0])
 %title('Type-III FUSTFT')