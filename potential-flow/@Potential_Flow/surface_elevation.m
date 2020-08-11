% Wed 14 Mar 10:15:13 CET 2018
%% compute surface elevation according to Bernoulli's law
%% note : this is likely very different from the true surface elevation,
%%        as streamline curvature causes a transverse pressure gradient
function zs = potential_surface_elevation(obj,idir,varargin)
	g = obj.g;
if (1)
	umag = obj.mag('u',varargin{:});
	% 1/2 u^2 + g zs + c
	zs = -1/(2*g)*umag.^2
else
	% note: this is not quite correct to mix the potential flow with the secondary flow
	n = obj.mesh.n();

	Us             = obj.mag('u');
	[ds_dx, ds_dy] = obj.dir('u');

	% gradient of the velocity magnitude
	[dUs_dx, dUs_dy] = obj.grad_mag('u');

	% gradient along streamline (gradient projected to streamline direction)
	% gradient orthogonal to streamline
	[dUs_ds, dUs_dn] = obj.xy2sn(dUs_dx,dUs_dy,'u');

	R = obj.streamline_radius_of_curvature();

	%dz_ds = 2*h*u/C.^2*du_ds;
	dz_ds = -2/g*Us.*dUs_ds;
	dz_dn = -1./g*Us.^2./R;
	dz_dn = 1.*dz_dn;
	dz_ds = 1.*dz_ds;

	if (1)
		dz_ds = reshape(dz_ds,n);
		dz_dn = reshape(dz_dn,n);
		dz_ds(:,1) = dz_ds(:,2);
		dz_dn(:,1) = dz_dn(:,2);
		dz_ds = dz_ds(:);
		dz_dn = dz_dn(:);
	end
	%dz_dn(:,2:end) = 0;
	% TODO - what about friction?

	% rotate back to xy for plotting
	[dz_dx, dz_dy] = obj.sn2xy(dz_ds, dz_dn,'u');
	
	% rotate to mesh-SN coordinates
	[dz_dS, dz_dN] = obj.mesh.xy2sn(dz_dx,dz_dy);

	dz_dS = reshape(dz_dS,n);
	dz_dN = reshape(dz_dN,n);
	
	% TODO do not just select left right but also top-bottom by idir
	if (nargin()> 1 && idir)
		% determine integration constants at inflow boundary
		% set level at grid point 1,1 to zero
		z1 = cumintL(dz_dN(:,1), hN(:,1));

		% integrate along mesh-SN coordinates
		z  = cumintL(dz_dS.', hS.').';
	else
		% determine integration constants at inflow boundary
		% set level at grid point 1,1 to zero
		z1 = cumintL(dz_dS(1,:)',1)';
		z  = cumintL(dz_dN,1);
	end

	% apply integration constants
	z  = bsxfun(@plus,z1,z);
end

end % potential_surface_elevation

