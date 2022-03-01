function [estimate] = Supervised_separation(mix,Param,Sparse,STFT,Spv,ext,L_type,musicSS,drawConv,normalize)
%Learning basis matrix T and adapt ILRMA and ILRMA+SP
%%%%%%%%%%%%%%% Start learning %%%%%%%%%%%%%%%%%
for n = 1:Param.ns
Learning_count = Spv(n).tone;
R = 0;   
 fprintf('Learning count:   ');
 for i = 1:Spv(n).tone
    % make learning sound file name
    filename = [Spv(n).file,Spv(n).rootname,num2str(i),ext];
    %学習
    switch L_type
        case 1
        [B,r] = rank1_learning(filename,Param,STFT); %rank-1学習
        case 2
        [B,r] = rankr_learning(filename,Raram,STFT); %rank-r学習
    end 
    
    if i == 1
       B_2_all = B;
    else
       B_2_all = horzcat(B_2_all,B); %1音ずつ学習した基底を繋げる
    end
    R = R + r; 
   
    fprintf('\b\b\b\b%4d', Learning_count);
    Learning_count = Learning_count -1;
 end
Spv(n).R = R;
Spv(n).B = B_2_all;
end

Param.R_x = zeros(1,Param.ns); %R_xは各音源の基底数を集めたベクトル
for n = 1:Param.ns
   Param.R_x(:,n) = Spv(n).R;
end
Param.R_max = max(Param.R_x); %各音源で最大の基底数に合わせる

if Param.R_max ~= 0
    B_all = ones(size(Spv(1).B,1),Param.R_max,Param.ns)/sqrt(size(Spv(1).B,1));    % 使用する基底行列の初期化
%     B_all = zeros(size(B_1,1),R_max,N);  %周波数＊基底＊音源数 0で初期化
end
%B_allにすべての基底行列を格納
for n = 1:Param.ns 
    B_all(:,1:Spv(n).R,n) = Spv(n).B;
end
%%%%%%%%%%%%%%%%%%%%%%% end learning %%%%%%%%%%%%%%%%%%%%%%%%%%%

 switch STFT.type
    case 1
        [S_1,~] = STFT_1(musicSS(1,:)',fftSize, shiftSize,STFT.w_name);
    case 2
        S_1 = STFT_basic(musicSS(1,:)',STFT.analysis_window, STFT.shiftSize, STFT.frag);
end
I=size(S_1,1);
J=size(S_1,2);

% Short-time Fourier transform
X=zeros(I,J,Param.mic);
for m=1:Param.mic
    switch STFT.type
        case 1
           [X(:,:,m), window] = STFT_1(mix(:,m), fftSize, shiftSize,STFT.w_name);
        case 2
            X(:,:,m) = STFT_basic(mix(:,m),STFT.analysis_window,STFT.shiftSize,STFT.frag);
    end
end 
% Whitening
Xwhite = whitening( X, Param.ns ); % decorrelate input multichannel signal by applying principal component analysis

V = rand(Param.R_max,J,Param.ns);

% 分離アルゴリズム
% [Y,cost,~] = ILRMA_readable(Xwhite, type, it, Param.R_max, drawConv, normalize,B_all,V);   
[Y1,~,A1,~,~] = ILRMA1_initial(X,B_all,V,Param);
[Y_SP,A_SP] = ILRMA1_Sparse(X,B_all,V,Param,Sparse,STFT.fftSize);

% For superised ILRMA
Z_PB = projectionBack(Y1,A1);
for m = 1:Param.mic
    estimate(m).Y_ILRMA = zeros(I,J,Param.ns);
    for n = 1:Param.ns
        estimate(m).Y_ILRMA(:,:,n) = Z_PB(:,:,n,m);
    end
end

% For superised ILRMA+SP
Z_SP_PB = projectionBack(Y_SP,A_SP);
for m = 1:Param.mic
    estimate(m).Y_ILRMA_SP = zeros(I,J,Param.ns);
    for n = 1:Param.ns
        estimate(m).Y_ILRMA_SP(:,:,n) = Z_SP_PB(:,:,n,m);
    end
end

%%%%% ISTFT %%%%%
for m=1:Param.mic
    estimate(m).T_signal_ILRMA = zeros(Param.sig_length,Param.ns);
    estimate(m).T_signal_ILRMA_SP = zeros(Param.sig_length,Param.ns);
    switch STFT.type
        case 1
          estimate(m).T_signal_ILRMA = ISTFT(estimate(m).Z_PB,shiftSize,window,size(mix,1));
          estimate(m).T_signal_ILRMA_SP = ISTFT(estimate(m).Z_SP_PB,shiftSize,window,size(mix,1));
        case 2
          for n = 1:Param.ns
          estimate(m).T_signal_ILRMA(:,n) = ISTFT_basic(estimate(m).Y_ILRMA(:,:,n),STFT.dual_window,STFT.shiftSize,Param.sig_length,STFT.frag);
          estimate(m).T_signal_ILRMA_SP(:,n) = ISTFT_basic(estimate(m).Y_ILRMA_SP(:,:,n),STFT.dual_window,STFT.shiftSize,Param.sig_length,STFT.frag);
          end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%