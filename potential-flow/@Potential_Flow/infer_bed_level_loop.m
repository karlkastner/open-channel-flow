% Thu 15 Mar 17:36:43 CET 2018
%
%% the bed level does not completely converge but starts to oscillate,
%% this is presumably due to the non-compact kernel implementation of the laplacian oberator
% 
function obj = infer_bed_level_loop(obj,h0,bcfun,cutfun)
%function obj = infer_bed_level_loop(obj,h0,u0,len,R0,bcfun,cutfun)
% ,phi_side)
	p     = 0.05;
	k_max = 30;

	% TODO use pade
	k = 0;
	%h = obj.h;
	%bc(1) = struct();
	%bc(1).x0 = [obj.mesh.X(1,1),obj.mesh.X(end,1)];
	%bc(1).y0 = [obj.mesh.Y(1,1),obj.mesh.Y(end,1)];
	% bc(1).f  = @phi_in;
	n = obj.mesh.n;

	h = h0*ones(prod(n),1);
	obj.h = h;

	figure(1e4)
	clf
	delta_old = inf;
	while (1)
		k = k+1;	
		% set inflow bc according to bed level
		
		obj.assemble_potential_matrix(bcfun, cutfun);
		%obj.assemble_discretization_matrix(
		% bc);
		%obj = side_inflow(obj,phi_side,obj.h);

		% recmompute potential
		obj.solve_potential();

		% recompute bed level
		h_old = obj.h;

		h = obj.infer_bed_level(h0);
		delta = norm(h-h_old)

		% update bed-level with over-relaxation
		obj.h = (1-p)*obj.h + p*h;

		% TODO hotfix
		obj.h = max(0.1,obj.h);
		obj.h = flat(cmean(cmean(reshape(obj.h,obj.mesh.n)')'));
		obj.h = flat(cmean(cmean(reshape(obj.h,obj.mesh.n)')'));

if (0)

		u = obj.u;
		v = obj.v;
		[uS, uN] = obj.mesh.xy2sn(u,v);

		figure(400+k);
		clf
		subplot(2,2,1)
		obj.plot(h);
		title('h')
		colorbar

		subplot(2,2,2)
		obj.plot(obj.h);
		title('relaxed h')
		colorbar

		subplot(2,2,3)
		obj.plot(h_old);
		title('h_old')
		colorbar()
		
		dh = h-h_old;
		subplot(2,2,4)
		%obj.plot(flat(dh)); %h-h_old);
		%title('h-h_old')
		obj.plot(uN);
		title('uN');
		colorbar()


		figure(1e4);
%		subplot(2,2,3)
		subplot(2,2,1)
		hh = reshape(h,n);
%		hhold = reshape(hold,n);
%		plot(-[mean(hh,2),mean(hhold,2)])
		plot(-mean(hh,2))
		hold on
		subplot(2,2,2)
		uS = reshape(uS,n);
		uS_ = mean(uS,2);
		plot(uS_)

		subplot(2,2,3)
		%uS_ = mean(uS,1);
		plot(uS(end-1,:))
		hold on
end



		if (delta > delta_old)
			break
			%pause
		end
		delta_old = delta
		
		if (k>=k_max)


		%umag = obj.mag('u');
%		clf();
%		subplot(2,2,2)
%		obj.plot(obj.phi);
%		title('phi')

%		subplot(2,2,3)

%		subplot(2,2,4)
		%obj.plot(uN);

		

	
			break;
		end
%		pause();
	end

%	function phi = phi_in(x,y)
%		 r     = hypot(x,y);
%		 h_in  = h(1:n(1),1);
%		 phi   = (h_in./h0).^0.*u0.*len.*(r./R0).^0;
%	end
end

