% Wed 15 Mar 09:07:25 CET 2017
%% reflection coefficients of shallow water waves at a sudden change of the
%% cross section (sudden change of admittance)
%% c.f. lighthill, ippen-harleman
function [ar, at] = sw_reflection(h1,w1,h2,w2)
	g = Constant.gravity;
	% wave numbers
	c1 = sqrt(g*h1);
	c2 = sqrt(g*h2);
	% reflection coefficient
	%ar = (c2*w2-c1*w1)/(c1*w1+c2*w2);
	ar = (c2*w2-c1*w1)/(c1*w1+c2*w2);
	%ar = -abs(c2*w2-c1*w1)/(c1*w1+c2*w2);
	%ar = abs(c2*w2-c1*w1)/(c1*w1+c2*w2);
	% transmission cofficient
	at = (2*c1*w1)/(c1*w1+c2*w2);
end

