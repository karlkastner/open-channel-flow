% 2019-07-28 03:40:54.927742302 +0800
function c_d = chezy2drag(C)
	g   = Constant.gravity;
	c_d = g./C.^2;
end

