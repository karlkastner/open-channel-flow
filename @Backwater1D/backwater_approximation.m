% 2017-05-07 12:01:26.644611139 +0200
% Karl Kastner, Berlin
%% approximation of the backwater curve by an exponential function
%% note: this is not necessarily a good approximation
%% in the case of tide, Qt can be given
function [h dh_dx L Ha] = backwater_approximation(Qr,C,W,S,H0,x,Qt)
	[L dh_dx0 Ha] = backwater_length(Qr,C,W,S,H0,Qt);
	% depth along channel (the backwater)
	h = H0 + (Ha-H0)*(1-exp(-x/L));
	% slope along channel
	dh_dx = (Ha-H0)/L*exp(-x/L);
end

