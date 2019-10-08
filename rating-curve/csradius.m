% Thu  9 Feb 13:37:15 CET 2017
% Karl Kastner, Berlin
%
%% compute hydraulic radius of the cross section
%
% function Rh = csradius(z_s, z_b, dn)
function Rh = csradius(z_s, z_b, dn)
	A = csarea(z_s,z_b,dn);
	P = csperimeter(z_s,z_b,dn);
	Rh = A./P;
end

