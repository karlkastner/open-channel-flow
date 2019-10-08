% 2016-02-20 22:42:05.379762544 +0800
%% wrapper function
function l = discharge2stage(Q,p)
	l = (Q./p(:,1)).^(1./p(:,3)) + p(:,2);
end

