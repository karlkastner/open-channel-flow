% Wed 29 Apr 19:57:13 +08 2020
function kc = manning2kc(n,h)
	C  = manning2chezy(n,h);
	z0 = chezy2z0(C,h); 
	kc = z02ks(z0);
end
