% 2016-04-08 11:15:55.270978358 +0200
% Karl Kastner, Berlin
%
%% change of depth along channel for the backwater equation
%
%  S0   : bed slope
%% beta : momentum coefficient
%% this is effectively an equation in h^3
%
% TODO improve for cases of high froude number and strong tidal flow
%      froude number correction only exact if tidal flow zero
%      and tidal friction only exact when froude number is low
% TODO interaction of tide-width is neglected
function dh_dx = dh_dx(obj,x,h)
	Q0  = obj.Q0;
	Qt  = obj.Qt(x);
	% TODO, proper tidal range
	Qhr  = sum(abs(Qt),2);
	Qmid = Q0;
	C    = obj.chezy(x,h);
	W    = obj.width(x);
	% dw_dx = obj.dw_dx(x);
	dzb_dx = obj.dzb_dx(x);
	%dh_dx = obj.dh_dx_(obj,h,Q0,Qt,C,W,dzb_dx);
	dh_dx = obj.dh_dx_(h,Q0,Qt,Qmid,Qhr,C,W,dzb_dx);
end


