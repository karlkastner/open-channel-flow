% Sun  7 May 12:09:35 CEST 2017
% Karl Kastner, Berlin
%%
%% analytic solution of the gradually varied flow equation
%% c.f. Bresse, Chow
function [x h dh_dx] = backwater_analytic(Q,C,W,S,H0,X,nn)
	g = 9.81;

	x = linspace(X(1),X(2),nn)';

	dx = x(2)-x(1);

	% integration constant:
	c = 1/4*H0^4*W^2 - H0*Q^2/g; % + S*H0^3*W^2;

%	figure(10)
%	clf

	maxiter = 100;
	h = normal_flow_depth(Q,W,C,S)*ones(nn,1);
	% relaxation constant
	p = 0.5;
	for jdx=1:maxiter
		%Hc = (1-p)*H__ + p*H_;
		% central differences
		%Hc = 0.5*([0; Hc(1:end-1)] + [Hc(1:end)]);
		% upper differences
		Hc = [0;h(2:end)];
		for idx=1:length(x)
			rhs = x.*Q^2/C^2 - S*dx*sum(Hc(1:idx).^3)*W^2;
			% this is only exact for S==0 and and has to be iterated otherwise (because rhs (constant term) is not independent on H for S not 0)
			h4(idx,:) = roots([1/4*W^2-1/4*S*W^2,0,0,-Q^2/g,-rhs(idx)-c]);
		end
	%	figure(11);
	%	subplot(4,4,jdx)
	%	plot(real(h))
	%	title(jdx)
	
	%	figure(10)
	%	subplot(4,4,jdx)
	%	cla
	%	h = real(h);
	%	h = sort(h')';
	%	h = h(:,end);
		%h_ = h(:,1);
		h_ = max(max(real(h4),0),[],2);
	%	plot(x,real(h(:,1)));
	%	hold on
	%	H(:,jdx) = h_(:,1);
	%	H__ = H_;
	%	H_ = h_(:,1);

		h = (1-p)*h + p*h_;
	end % for jdx
	h = h_;
end

