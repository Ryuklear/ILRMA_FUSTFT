function [Y_hat] = projectionBack(Y,A)

I = size(Y,1);
J = size(Y,2);
N = size(A,2);
M = size(A,1);

Y_hat =zeros(I,J,N,M);
for i = 1:I
    for j = 1:J
        for n = 1:N
            Y_hat(i,j,n,:) = A(:,n,i)* Y(i,j,n);
        end
    end
end



% % function [Z, D] = projectionBack(Y,X)
% % %
% % % Corded by D. Kitamura (d-kitamura@nii.ac.jp) on 9 Aug. 2016.
% % %
% % % This function restores the scale of the signals estimated by ICA-based
% % % blind source separation techniques.
% % %
% % % see also
% % % http://d-kitamura.sakura.ne.jp/index.html
% % %
% % % [syntax]
% % %   [Z, D] = projectionBack(Y, X)
% % %
% % % [inputs]
% % %   Y: estimated (separated) signals (freq. x frames x sources)
% % %   X: observed (mixture) signal with desired channel (freq. x frames x 1)
% % %      or observed multichannel signals (freq. x frames x channels)
% % %
% % % [outputs]
% % %   Z: scale-fitted estimated signals (freq. x frames x sources)
% % %      or scale-fitted estimated source images (freq. x frames x sources x channels)
% % %   D: transformation matrix
% % %
% % 
% % % check errors
% % if (nargin<2)
% %     error('Too few input arguments.\n');
% % end
% % [I,J,M] = size(Y); % freq. x frames x channels
% % 
% % % projection back
% % if (size(X,3)==1) % calculate scale-fitted estimated signals to X(:,:,1)
% %     A = zeros(1,M,I);
% %     D = zeros(M,M,I);
% %     Z = zeros(I,J,M);
% %     for i=1:I
% %         Yi = squeeze(Y(i,:,:)).'; % channels x frames (M x J)
% %         A(1,:,i) = X(i,:,1)*Yi'/(Yi*Yi');
% %     end
% %     for m=1:M
% %         for i=1:I
% %             Z(i,:,m)=A(1,m,i)*Y(i,:,m);
% %             D(m,m,i)=A(1,m,i);
% %         end
% %     end
% % elseif (size(X,3)==M) % calculate source image of the estimated signals
% %     A = zeros(M,M,I);
% %     Z = zeros(I,J,M,M); % freq. x frames x sources x channels
% %     for i=1:I
% %         for m=1:M
% %             Yi = squeeze(Y(i,:,:)).'; % channels x frames (M x J)
% %             A(m,:,i) = X(i,:,m)*Yi'/(Yi*Yi');
% %         end
% %     end
% %     for n=1:M
% %         for m=1:M
% %             for i=1:I
% %                 Z(i,:,n,m)=A(m,n,i)*Y(i,:,n);
% %             end
% %         end
% %     end
% %     D = A;
% % else
% %     error('The number of channel in X must be 1 or equal to that in Y.\n');
% % end
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%