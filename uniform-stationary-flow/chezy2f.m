% Wed 25 Mar 10:03:31 +08 2020
function f = chezy2f(C)
	if (issym(C))
		syms g
	else
		g = Constant.gravity;
	end
	f = 2*g./C.^2;
end


