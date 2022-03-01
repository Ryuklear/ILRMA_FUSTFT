function[B,r] = rankr_learning_FUSTFT(filename,Param,FUSTFTpara)
%楽器の単音を基底行列として抽出(ほぼ論文通り) 

fs = Param.fsResample;
L = FUSTFTpara.fftSize;
eta = FUSTFTpara.shiftSize;

[music1,Fs] = audioread(filename);
music1 = resample(music1,fs,Fs);
T1 = size(music1,1);  
T2 = fix((T1 - L + eta)/eta) * eta + L;
music1 = music1(1:T2);

switch FUSTFTpara.type
    case 1
    [X_comp,~] = FUSTFT(music1',L,eta,FUSTFTpara.w_name); 
    case 2
     X_comp = STFT_basic(music1',FUSTFTpara.analysis_window, eta, FUSTFTpara.frag);   
end
X_mono = abs(X_comp).^2;   %我々の研究ではパワースペクトル!
J = size(X_mono,2); %行の長さ（横）

%SVDを施す
[U,S,V] = svd(X_mono); 
%diag()では対角成分を取得（ここでは特異値取得）
s = diag(S);
%ランク数を決めるための誤差範囲(10%以下)
epsilon = norm(X_mono,'fro')/10;

%ランク数ｒ(Bベクトル数)を決定
for r = 1:(length(s)-1)
    %何本目で誤差範囲を満たすか
    %[.]は各要素ごとの意味
    if sqrt(sum(s(r+1:end).^2)) <= epsilon   %フロベルニウスノルム?
        break
    end
end

%基底行列Bの設計
B = zeros(size(X_mono,1),r);
%Bの1列目はUの1列目
B(:,1) = abs(U(:,1));

%X_monoのランクr近似
X_mono_r = U(:,1:r)*S(1:r,1:r)*(V(:,1:r)');

%Bの2列目以降はX_mono_rから取得する
for k = 2:r
    PB_x = zeros(size(X_mono,1),J);
    Y = zeros(1,J);
    
   for i = 1:J
       PB_x(:,i) = B(:,1:k-1)*((B(:,1:k-1)'*B(:,1:k-1))\(B(:,1:k-1)'*X_mono_r(:,i)));
       Y(1,i) = abs(pi/2 - acos(PB_x(:,i)'*X_mono_r(:,i)/(norm(PB_x(:,i))*norm(X_mono_r(:,i)))));
   end
    [~, I] = min(Y);
    B(:,k)= X_mono_r(:,I)/norm(X_mono_r(:,I));
end

%Bの負の値を0に
B(B<0)=0;
%norm(B,'fro')^2

for k=2:r
    B(:,k) = B(:,k)/norm(B(:,k));
end