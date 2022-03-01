function[B,r] = rank1_learning_FUSTFT(filename,Param,FUSTFTpara)
%楽器の単音を基底行列として抽出
r = 1;         %ランク１ 
fs = Param.fsResample;
L = FUSTFTpara.fftSize;
eta = FUSTFTpara.shiftSize;

[music1,Fs] = audioread(filename);
% T = size(music1,1);  Fs時のサンプル数（時間はT/Fs[s]）

music1 = resample(music1,fs,Fs);
T1 = size(music1,1);   %music1の総サンプル数
%music1のSTFTちょうどの長さ
T2 = fix((T1 - L + eta)/eta) * eta + L;

music1 = music1(1:T2);

switch FUSTFTpara.type
    case 1
    [X_comp,~] = FUSTFT(music1',L,eta,FUSTFTpara.w_name); %複素数値が求まる
    case 2
     X_comp = STFT_basic(music1',FUSTFTpara.analysis_window, eta, FUSTFTpara.frag);
end
X_mono = abs(X_comp).^2;   %我々の研究ではパワースペクトル!

I = size(X_mono,1); %行列の列の長さ（縦）
J = size(X_mono,2); %行の長さ（横）

%SVDを施す
[U,~,~] = svd(X_mono); 

%基底行列Bの設計
B = zeros(size(X_mono,1),r);
%Bの1列目はUの1列目
B(:,1) = abs(U(:,1));

%Bの負の値を0に
B(B<0)=0;
%norm(B,'fro')^2

    B(:,r) = B(:,r)/norm(B(:,r));
end