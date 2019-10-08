% 2017-05-07 12:01:26.644611139 +0200
%% backwater length
function [L, dh_dx0, Ha] = backwater_length(Qr,C,W,S,H0,Qt)
	g = Constant.gravity;
	% asymptotic depth (limit of depth for x->infinity)
	Ha = (Qr^2/(C^2*S*W^2))^(1/3);
	if (nargin() < 6)
		% friction for river flow only
		Sf = Qr^2/(C^2*H0^3*W^2);
	else
		% friction for river and tidal flow
		% TODO this is not a good approximation
		% Qt decreases along channel
		Sf = 8/(3*pi)*Qr*(Qr+Qt)/(C^2*H0^3*W^2);
	end
	% TODO what is with the Froude number?
	dh_dx0 = (Sf - S)/(1-Qr^2/(g*W^2*H0^3));
	%dh_dx0 = (Q^2/(C^2*H0^3*W^2) - S)/(1-Q^2/(g*W^2*H0^3))

	% backwater length (equal to length where slope at origin intersects asymptotic depth)
	L = (Ha - H0)/dh_dx0;
end

