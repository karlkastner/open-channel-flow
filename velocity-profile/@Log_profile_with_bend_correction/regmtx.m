% Thu 20 Sep 18:29:00 CEST 2018
%% regression matrix
function A = regmtx(obj,z,h)
	s = z/h;
%	p = obj.p;
	A = [log(z), ones(size(z)), ...
	     2*s-4*s.^3, 3*s.^2-4*s.^3];
%(2/h^p*Z-3./h^(p+1)*Z.^2)];
%			A = [ln_Z(binmask(:,idx),idx), ...
%			     o(binmask(:,idx)), ...
%			     (2/h(idx)^p*Z(binmask(:,idx),idx)-3/h(idx)^(p+1)*Z(binmask(:,idx),idx).^2)];
%			      2*s - 3*s^2
end
