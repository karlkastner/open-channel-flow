% 2016-04-08 11:15:55.270978358 +0200
%
%% change of surface elevation along channel
function dzs_dx = dzs_dx(obj,x,zs,Q,C,S0,W)
	% bed level
	zb = S0*x;

	% flow depth
	h  = zs-zb;
		
	% change of flow depth
	dh_dx = obj.dh_dx(x,h,Q,C,S0,W);

	% change of surface elevation
	dzs_dx = dh_dx + S0;
end

