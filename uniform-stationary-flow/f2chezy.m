% Wed 25 Mar 10:05:49 +08 2020
function C = f2chezy(f)
	if (issym(f))
		syms g positive
	else
		g = Constant.gravity;
	end
	C = sqrt(2*g./f);
end
