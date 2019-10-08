% Thu Jul 17 13:16:00 WIB 2014
%% regression matrix
function A = regmtx(obj,Z,~)
	% log law: u = us/k ln(z) - us/k ln(z0)
	% TODO use S instead of Z and h
	A     = [log(Z), ones(size(Z))];
	%A     = [ln_Z(binmask(:,idx),idx), ...
       	%         o(binmask(:,idx)) ];
end
