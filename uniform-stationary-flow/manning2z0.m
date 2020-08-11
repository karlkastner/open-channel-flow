% Wed  8 Jan 18:24:59 +08 2020
function z0 = manning2z0(n,h)
	C = manning2chezy(n,h);
	z0 = chezy2z0(C,h); 
end


