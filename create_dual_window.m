function dual_window = create_dual_window(analysis_window,time_skip)
%CREATE_DUAL_WINDOW creates the canonical dual window for an analytic window
%
%% -- Inputs ---------------------------------------------------------
% analysis_window:   (positive) analysis window (window_length x 1)
% time_skip:   shift size of time frame (1 x 1) (time_skip must be equal to or smaller than window_length)
%
%% -- Temporary Variable ---------------------------------------------
% window_length:   (support) length of the analysis and dual windows (1 x 1)
% 
%% -- Output ---------------------------------------------------------
% dual_window:   canonical dual window (window_length x 1)
%
%% start program

window_length = length(analysis_window);
dual_window = [analysis_window; zeros(time_skip * ceil(window_length/time_skip) - window_length,1)];
dual_window = reshape(dual_window,time_skip,[]);
dual_window = dual_window./sum(abs(dual_window).^2,2);
dual_window = reshape(dual_window(1:window_length),[],1);

end