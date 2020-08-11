% Mon 21 Oct 17:42:06 +08 2019
% 2018-02-19 14:58:48.380392520 +0100
%function C = manning2chezy(n,h)
%% manning to chezy conversion
function C = manning2chezy(n,R)
	C = 1./n.*R.^(1/6);
end

