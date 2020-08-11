% Sun 17 May 21:35:39 +08 2020
function [mesh,fun] = test_case(id,n)

if (nargin()<2)
	n = 10;
end

n = n*[1,1];

mesh = StructuredMesh();

switch (id)
case {1}
	mesh.generate_rectangle([0,1],[0,1],n);
	% flow along y
	fun.u = @(x,y) ones(size(x));
	fun.v = @(x,y) zeros(size(x));
%	x = linspace(0,1,n)';
%	y = linspace(0,1,n)';
%	XX = repmat(x,m);
%	YY = repmat(y,m);
%	%s = hypot(cdiff(x),cdiff(y));
case {2}
	% flow along curve
	mesh.generate_disk([0.5,1.5],[0,2*pi/4],n);
	fun.u = @(x,y)  y.*1./(hypot(x,y)).^2;%.*sin(atan(y,x));
	fun.v = @(x,y) -x.*1./(hypot(x,y)).^2;%.*sin(atan(y,x));
	% flow along x
%	y = sc*linspace(0,1);
%	x = sc*0*y; %zeros(size(x));
%	s = hypot(cdiff(x),cdiff(y));
%case {3}
%case {4}
%	t = linspace(0,1/4);
%	x = sc*sin(2*pi*t);
%	y = sc*cos(2*pi*t);
%	s = hypot(cdiff(x),cdiff(y));
% semicircle
%case {5}
	%w = sc*0.5*ones(size(x));
	%ds = sc*0.05;
	%nn = 7;
	%pf.mesh.generate_from_centreline(x,y,w,ds,nn);
end
	mesh.X = mesh.X';
	mesh.Y = mesh.Y';

	bc      = struct();
	bc(1).y = [mesh.X(  1,1), mesh.X(1,end)];
	bc(2).y = [mesh.X(end,1), mesh.X(end,end)];
	bc(1).x = [  mesh.Y(1,1),mesh.Y(1,end);]
	bc(2).x = [mesh.Y(end,1),mesh.Y(end,end)];
	bc(1).p = 0; % 0 set flow : 1 upstream
        bc(2).p = 1;  % potential (level) : 0 downstream
        bc(1).v = 1;  % unit inflow velocity 
        bc(2).v = 0;   % zero level
	bc(1).d_tol = sqrt(eps);
	bc(2).d_tol = sqrt(eps);
	mesh.boundary_condition_s = bc;
end

