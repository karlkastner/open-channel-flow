% Sat  7 Jan 14:24:42 CET 2017
%% log law with wake
%% u = us/k ln(z) - us/k ln(z0) + us/k (2/H^2 z  - 3/H^3 z^2)
% TODO use orthogonal (trans)-form
function A = regmtx(obj,Z,h)
%	p = obj.p;
	p = 1;
	A = [log(Z), ones(size(Z)), (2/h^p*Z-3./h^(p+1)*Z.^2)];
%			A = [ln_Z(binmask(:,idx),idx), ...
%			     o(binmask(:,idx)), ...
%			     (2/h(idx)^p*Z(binmask(:,idx),idx)-3/h(idx)^(p+1)*Z(binmask(:,idx),idx).^2)];
end
