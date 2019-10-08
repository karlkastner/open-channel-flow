% Sa 6. Feb 23:10:24 CET 2016
% Karl Kastner, Berlin
%
%% linearised SWE
%% width variation not included, goes into rhs force term
%%
%% [      0,  1] [A]    = [    Q]
%% [ -u^2+gH, 2u] [Q]_dx   [Q^2/A+1/2gA^2/w]_dx - 1/2gA^2/w^2 dw/dx
%%                                               force term
%%
function L = lindot(obj,t,A,Q)
	U = Q./A;
	H = A./obj.width;
	g = Constant.g;
	L = zeros(2,2,obj.nx);
	L(1,2,:) =          1;
	L(2,1,:) =  -U.*U+g*H;
	L(2,2,:) =        2*U;

