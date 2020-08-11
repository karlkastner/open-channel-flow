% Mon 19 Feb 16:07:09 CET 2018
% Mon 21 Oct 17:42:14 +08 2019
%% convert chezy to manning
function n = chezy2manning(C,R)
	n = 1./C.*R.^(1/6);
	if (~issym(R))
		n(R<0) = NaN;
	end
end


