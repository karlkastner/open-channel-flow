% 2016-02-20 22:41:54.532065611 +0800
% wrapper
% TODO quick fix
function Q = stage2discharge(l,p)
	Q = p(:,1)*(max(0,l-p(:,2))).^p(:,3);
end

