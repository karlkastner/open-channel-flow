% 2017-11-06 00:22:36.341370088 +0100
%
%% friction source term of the SWE
% y_bar : roe average at cell interfaces
function f = source_friction(obj, t, x, y_bar, w_bar)
	n = length(y_bar)/2;

	% TODO this should be c_bar (!)
	cd = obj.cd;
	
	% extract
	% if QAflag, then y1 = A, else h
	y1 = y_bar(1:n);
	q  = y_bar(n+1:end);
	
	if (obj.QAflag)
		w = w_bar;	
	else
		w = 1;
	end

	% friction forcing term
	% TODO this is a wide channel approximation
	f = [zeros(n,1);
             -cd.*w.*q.*abs(q)./y1.^2];
end % source_friction

