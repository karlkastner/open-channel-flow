% Tue 29 Oct 10:02:22 +08 2019
function cD = manning2drag(n,h)
	C = manning2chezy(n,h);
	cD = chezy2drag(C);
end

