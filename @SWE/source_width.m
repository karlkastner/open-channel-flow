% Sat 11 Nov 23:15:30 CET 2017
%% source term (reaction term) for channels with variable width
% y_bar : roe average at cell interfaces
% TODO take left as index, do not automatically assume up(left)
% TODO, this does not seem right, the width source term should be zero,
%       if in left hand side derivatives are with respect to A and Q
function [f obj] = source_width(obj,t,x,y_bar,w_bar)
	n  = length(x);

	% ldiff, not cdiff!
	w = obj.w;
	dw_dx = leftdiff(w)./leftdiff(x);
	if (~obj.QAflag)
%		h  = y(1:n);
%		q  = y(n+1:2*n);

		% extract
		h_bar = y_bar(1:n);
		q_bar = y_bar(n+1:2*n);


		f  = [-q_bar.*dw_dx./w_bar;              % source term for continuity equation
        	      -q_bar.*q_bar./h_bar.*dw_dx./w_bar % source term for for momentum equation
		     ];
	else
		g = obj.g;
		A_bar = y_bar(1:n);
		%Q = y_bar(n+1:2*n);
		% c.f zoppu 2003
%		f = [zeros(n,1);
%                     A_bar./w_bar.*dwdx];
%                     %A_bar./w_bar.*dwdx];

%	f = [-Q./wbar.*dw_dx;
%	     -(c.^2.*A - Q^2./A).*(dw_dx./w_bar)];
		f = [zeros(n,1);
                     g*(A_bar./w_bar).^2.*dw_dx];
		     
	end
end

