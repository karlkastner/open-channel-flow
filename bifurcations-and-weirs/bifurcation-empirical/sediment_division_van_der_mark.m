% Wed  8 Aug 20:33:05 CEST 2018
% ps = S_s/S_0
function [ps,r,theta,R,d_] = sediment_division_van_der_mark(Q0,Qs,w0,ws,h)
	b = 5;
	c = 1.35;
	A = 10;

	lambda = w0/pi;
	eta    = 2*lambda;
%	R      = bifurcation_separation_line_radius(w0,Qs/Q0);
	R      = bifurcation_potential_flow_dividing_streamline_radius(w0,Qs/Q0);

	theta  = asin( (lambda+eta)/R ) - asin(eta/R);
	delta  = atan(A*h/R);
	% ps = S_s/S_0
	ps     = Qs/Q0 + ((c*Qs*w0)/(Q0*ws))^(b/3)*(R*theta)/w0*(1/cos(delta)-1);
	d_     = (1/cos(delta)-1);
	% r
	% dummy concentration
	%c   = 1;
	%S0  = c*Q0;
	%Ss  = S0*ps;
	%r   = (Q0/S0)*(Ss/Qs);
	% r   = (Ss/Qs)*(S0/Q0) = ps/pq
	pq = Qs/Q0;
	% r = c
	r  = ps/pq;
	%c*Q0/Qs)*ps;
end

