% Thu 15 Mar 16:47:25 CET 2018
%
%% note: this is pretty much a broken function for the inference of stationary
%%       morphology
%%
%% Missing:
%% - rolling down of transverse slope to balance secondary flow in bends
%% - quasi time steippong
%%
%% at stationary state:
%% - changes of discharge along the streamlines of discharge are balanced
%%   by a change in depth, to keep the velocity and sediment transport constant along the streamline
%%
%% dz_b/dt = dqs/dx + dqs/dn = 0                    (i)
%% TODO this only true for infinite bends, as sediment can also move to the side
%% dqs/ds  = d/s(q/h) = 1/h dq/ds - q/h^2 dh/ds = 0
%% TODO this is only true in an ifinite bend (ikeda)
%% dqs/dn = 0
%% streamlines along discharge or velocity -> does not matter eq (i) is direction independent
function h = infer_bed_level(obj,h0)
%	qs = h0*obj.umag;

	% for continuous transport of sediment along streamline
	% streamline is hear in the direction of discharge
	% du/ds = 0
	% 1/h dq/ds - q/h^2 dh/ds = 0
	% dh/ds = h/q dq/ds

	% gradient of discharge
%	[dq_dx, dq_dy] = obj.grad_q();

%	[dUs_dx, dUs_dy] = obj.grad_umag();
%	dq_dx = h0*dUs_dx;
%	dq_dy = h0*dUs_dy;

	[dumag_dx, dumag_dy] = obj.grad_mag('u');

	% gradient along streamline of velocity
	[du_ds, du_dn] = obj.xy2sn(dumag_dx,dumag_dy,'u');

	% required change of depth to keep transport in balance along streamline of velocity
	dh_ds = 1./obj.mag('u').*du_ds;

	% for no transport across streamline
	% F_slope = F_secondary_flow
	R      = obj.streamline_radius_of_curvature();
	dh_dn  = dh_dn_bend_ikeda(obj.mag('u'));
	
	% transform coordinates
	[dh_dx, dh_dy] = obj.sn2xy(dh_ds, dh_dn, 'u');
	[dh_dS, dh_dN] = obj.mesh.xy2sn(dh_dx,dh_dy);

	n     = obj.mesh.n;
	dh_dS = reshape(dh_dS,n);
	dh_dN = reshape(dh_dN,n);

	hS = 1;
	hN = 1;

	dh_dS(:,1)   = dh_dS(:,2);
	dh_dS(:,end) = dh_dS(:,end-1);

	if (0)
		% pseudo time-step
		% TODO this reqires sediment transport, not dh_dS
		dt = 0.1;
		h = reshape(obj.h,obj.mesh.n);
		h = h-dt*dh_dt;
		% integration constant at inflow
		h = h + (h0 - h(end/2,1));
		h = flat(h);

	else 
	
		% TODO do not just select left right but also top-bottom by idir
		if (0) %nargin()> 1 && idir)
			% determine integration constants at inflow boundary
			% set level at grid point 1,1 to zero
			h1 = cumintL(dh_dN(:,1), hN(:,1));
			%h1 = h1 + (h0-h1(end/2));
	
			% integrate along mesh-SN coordinates
			h  = cumintL(dh_dS.', hS.').';
		else
			% determine integration constants at inflow boundary
			% set level at grid point 1,1 to zero
			h1 = cumintL(dh_dS(1,:)',1)';
	%		h1 = h1 + (h0-h1(end/2));
			h  = cumintL(dh_dN,1);		
		end
	
	
		% apply integration constants
		h  = bsxfun(@plus,h1,h);
		% determine integration constant
		h = h + (h0-h(ceil(end/2),1));
	
		h = h(:);
	end
end

