% Sat 19 May 12:37:44 CEST 2018
%
%% objective function for determining the bed level
function f = objective_bed_level(obj,h)
	obj.h = h;
	obj.assemble_potential_matrix();
	obj.solve_potential();

	[qsx,qsy] = obj.sediment_transport(h);
%	qsx = obj.mesh.smooth(qsx);
%	qsy = obj.mesh.smooth(qsy);

	% 0 = d/dx qsx + d/dy qsy
	% this actually does not need to be squared
%	f  = (obj.mesh.Dx*qsx + obj.mesh.Dy*qsy).^2;
	f  = obj.mesh.Dx*qsx + obj.mesh.Dy*qsy;
%	f  = obj.mesh.smooth(f);
%	figure()
%	obj.plot(f)
%pause
end % objective_bed_level

