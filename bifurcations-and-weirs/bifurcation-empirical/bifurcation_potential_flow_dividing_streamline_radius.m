% Sat 28 Jul 12:25:47 CEST 2018
% van der mark 2013
%
% function R = bifurcation_separation_line_radius(w0,Qs_div_Q0)
function R = bifurcation_separation_line_radius(w0,Qs_div_Q0)
	% eq 6
	lambda   = w0/pi;
	% following 12
	eta      = 2*lambda;
	% eq 8
	b2       = w0*Qs_div_Q0;
	% following eq 7
	b_eta    = b2;
	% following eq 7
	b_lambda = 0.63*b_eta;
	% eq 7
	R = sqrt( (b_lambda/2 + lambda/(2*b_lambda)*(2*eta + lambda)).^2 + eta^2 );
end

