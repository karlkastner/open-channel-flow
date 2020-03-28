% 2017-11-05 21:50:32.999822433 +0100
%
%% source term of the SWE caused by a change of the bed level
%%
%% Note: this term causes splitting and averaging methods to fail to
%%       give accurate predictions of the smooth surface at steps of the bed
%
% TODO store x in SWE?
% TODO pass hql and hqr, xl and xr?
% TODO pass left operator instead of just performing left differences
% TODO allow for dynamic depth (exner)
%
% y_bar : roe average at cell interfaces
%
% function [f obj] = source_bed_level(obj, t, x, y_bar, ~)
%
function [f obj] = source_bed_level(obj, t, x, y_bar, ~)
	g     = obj.g;
	zb    = obj.zb;
	n     = length(y_bar)/2;

	% note:
	% if unknowns are QA, then y1 is A, else it is h
	% the implementation remains identical
	y1     = y_bar(1:n);

	% ldiff, not cdiff!
	dzbdx = leftdiff(zb)./leftdiff(x);

	% source term	
	f      = [zeros(n,1);
	          -g.*y1.*dzbdx];
end % source bed level

