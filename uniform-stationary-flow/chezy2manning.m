% Mon 21 Oct 17:42:14 +08 2019
function n = chezy2manning(C,h)
	n = 1./C.*h.^(1/6);
	if (~issym(C))
		n(h<0) = NaN;
	end
end


