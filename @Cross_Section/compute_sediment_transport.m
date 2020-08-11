% Mon  9 Dec 19:05:11 +08 2019

function obj = compute_sediment_transport(obj,fdx_zs)
	g  = Constant.gravity;

	eps_h = 0; %sqrt(eps);
	
	nt = length(obj.zs);
	nn = obj.nn;

	if (size(obj.n.q,1) ~= nt)
		fdx_zs = isfinite(obj.zs);
		fdx_   = repmat(fdx_zs,nn,1);
	else
		fdx_ = true(nt,nn);
	end
	if (size(obj.slope,2) ~= obj.nn)
		obj.slope = repmat(cvec(obj.slope),1,obj.nn);
	end

	% prepare computation
	obj.n.h            = max(eps_h,   (repmat(cvec(obj.zs),1,nn) ...
				   - (obj.n.zb)));
	obj.flat(fdx_);

%	tcs = tcs.(name);
	% velocity
	fdx              = obj.n.h<=0;
	obj.n.u          = obj.n.q./(obj.n.h); %.*width);
	obj.n.u(fdx)     = 0;
	obj.n.C(fdx)   = 0; % sqrt(eps);

	% froude number
	obj.n.Fr      = obj.n.u./sqrt(g.*obj.n.h);
	obj.n.Fr(fdx) = 0;

		obj.n.T  = transport_stage_rijn(obj.n.d50,obj.n.d90,obj.n.h,obj.n.u,obj.T_C);

		% bed form dimensions
		[obj.n.dune.vr.h,obj.n.dune.vr.l]  = bedform_dimension_rijn(obj.n.h,obj.n.d50,obj.n.T);
		rgh.dune       = bedform_roughness_rijn(obj.n.dune.vr.h,obj.n.dune.vr.l,obj.n.h);                                           
		rgh.grain      = grain_roughness_rijn(obj.n.d90,obj.n.h); 
		rgh.total      = total_roughness_rijn(obj.n.dune.vr.h,obj.n.dune.vr.l,obj.n.d90,obj.n.h);
		obj.n.vr.C.total   = double(rgh.total.C);
		obj.n.vr.C.dune    = double(rgh.dune.C);
		obj.n.vr.C.grain   = double(rgh.grain.C);

	% shear velocity
	obj.n.us      = shear_velocity(obj.n.u,obj.n.C);
	obj.n.us(fdx) = 0;

	width = 1;

	% sediment rating curves
	
		obj.n.qs.karim = total_transport_karim(  obj.n.C, ...
							obj.n.u, ...
							obj.n.h, ...
							width, ...
							obj.n.d50, ...
							obj.T_C);
		
		obj.n.dune.karim.h = dune_height_karim(obj.n.C,obj.n.u,obj.n.h, ...
							obj.n.d50,obj.T_C);
		obj.n.qs.wp = suspended_transport_wright_parker( ...
							obj.n.C, ...
							obj.n.u, ...
							obj.n.h, ...
							width, ...
							obj.n.d50, ...
							obj.n.d90, ...
							[],[],obj.T_C);

		% van rijn transport
		obj.n.a.vr    = reference_height_rijn(obj.n.dune.vr.h,obj.n.h);
		obj.n.ca.vr   = reference_concentration_rijn(obj.n.d50,obj.n.a.vr,obj.n.T,obj.T_C);

		obj.n.dsus.vr = suspended_grain_size_rijn(obj.n.d50,obj.dsd,obj.n.T);
		obj.n.ws.vr   = settling_velocity(obj.n.dsus.vr,obj.T_C,'rijn')

		% rouse number
		obj.n.psus.rijn   = suspension_parameter(obj.n.ws.vr,obj.n.us,obj.n.ca.vr,'rijn');
		obj.n.iF.vr        = reference_to_flux_averaged_concentration_rijn(obj.n.a.vr./obj.n.h,obj.n.psus.rijn);
		obj.n.qs.vr.total = double(total_transport_rijn(obj.n.C,obj.n.d50,obj.n.d90,obj.dsd,obj.n.u,obj.n.h,width,obj.T_C));
		obj.n.qs.vr.sus   = double(suspended_transport_rijn(obj.n.C,obj.n.d50,obj.n.d90,obj.dsd,obj.n.u,obj.n.h,width,obj.T_C));
		obj.n.qs.vr.bed   = double(bed_load_transport_rijn(obj.n.C,obj.n.d50,obj.n.d90,obj.n.u,obj.n.h,width,obj.T_C));

		obj.n.ws.default  = settling_velocity(obj.n.d50,obj.T_C);
		obj.n.psus.rouse  = suspension_parameter(obj.n.ws.default,obj.n.us,[],'rouse');

		% this used to be t.dm, but EH use d50
		obj.n.qs.eh   = total_transport_engelund_hansen(obj.n.C,obj.n.d50,obj.n.u,obj.n.h,width);
		obj.n.qs.mpm  = bed_load_transport_mpm(obj.n.u,obj.n.C,obj.n.d50,width);
		obj.n.qs.wu   = total_transport_wu(obj.n.C,obj.n.d50,obj.n.u,obj.n.h,width,obj.T_C);
		obj.n.qs.aw   = total_transport_ackers_white(obj.n.C,obj.n.d50,obj.n.u,obj.n.h,width,obj.T_C);
		obj.n.qs.yang = total_transport_yang(obj.n.C,obj.n.d50,obj.n.u,obj.n.h,width,obj.T_C);
		obj.n.qs.bagnold = total_transport_bagnold(obj.n.C,obj.n.d50,obj.n.u,obj.n.h,width,obj.T_C);

		% reshape fields : rows : days, columns : cross-section
		obj.n = structfun_deep(@(x) reshape_conditional(x,[nt,nn]), obj.n);

end % compute_sediment_transport

