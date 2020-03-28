% Thu Apr 28 20:29:22 MSD 2011
% Karl Kästner
%
%% roe average for the SWE
%
% TODO, should be part of the SWE
% TODO take left as index or matrix to avoid double computation
function [y_bar w_bar obj] = roe_average(obj,y,yl)
	g = obj.g;
	n = length(y)/2;
	
	if (obj.QAflag)
		% area
		A  = y(1:n);
		Al = yl(1:n);

		% depth
		h  = A./obj.w;
		% TODO this should not be automatically left, left should be an index or extrapolation matrix
		hl = Al./up(obj.w);

		% velocity
	        u  = y(n+1:end)./A;
        	ul = yl(n+1:end)./Al;
	else
		% depth
	        h  = y(1:n);
	        hl = yl(1:n);
		% velocity
	        u  = y(n+1:end)./h;
        	ul = yl(n+1:end)./hl;
	end

	% celerity
	c  = sqrt(g*h);
	cl = sqrt(g*hl);

	% roe average of velocity
        u_bar = ( cl.*ul + c.*u ) ./ ( cl + c );

	% roe averge of depth (15.35, 21.31 in randall-leveque)
	% follows from 15.28, 15.32
	if (obj.QAflag)
		A_bar  = 0.5*(Al+A);
		% roe average of discharge
		Q_bar  = A_bar.*u_bar;
		% assemble
		y_bar  = [A_bar;Q_bar];
	else
		h_bar  = 0.5*(hl+h);
        	q_bar  = h_bar.*u_bar;
		% assemble
		y_bar  = [h_bar; q_bar];
	end

	% average of width
	if (nargout()>1)
		% TODO this should not be automatically left, left should be an index or extrapolation matrix
		w_bar = leftmean(obj.w);
	end
end % roeAverage

