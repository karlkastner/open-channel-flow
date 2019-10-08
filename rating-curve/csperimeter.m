% 2015-03-16 19:16:15.139306769 +0100
% Karl Kastner, Berlin
%
%% compute wetted perimeter
%
% function P = csperimeter(z_s, z_b, dn)
% 
% z_s : [1,nt]  : water surface level
% z_b : [nn,nt] : bed level
%
function P = csperimeter(z_s, z_b, dn)
	ntz = length(z_s);

%	bottom = [-2*max(h); -max(h); bottom; -max(h); -2*max(h)];
%	% exten boundary by 45deg element
%	bottom = [bottom(1)-dx; bottom; bottom(end)-dx];
	z_s = rvec(z_s);
	if (isvector(z_b))
		z_b = cvec(z_b);
	end
	ntb = size(z_b,2);
	%nn = size(z_b,1);
	if (isscalar(z_s))
		z_s = repmat(z_s,1,ntb);
	end
	if (isvector(z_b))
		z_b = repmat(z_b,1,ntz);
	end

	% extend banks by rectangular walls
	% TODO extend slanting
	%max_zs = max(z_s);
	z_b    = [rvec(z_s);
                  z_b;
                  rvec(z_s)];

	% depth across section
	D = bsxfun(@minus,z_s,z_b);

	% zero depth in dry parts
	D_ = max(0,D);

if (0)
	% this is less accurate, as it does not trim the depth at 0
	dz = cdiff(D);
	n = length(z_b);
	dP = hypot(D,dn*ones(n,1));
	P = sum(dP,D>0);
else
	% change of width
	dn  = dn*ones(size(z_b,1)-1,ntz); % was +1

	% change of bed elevation
	dz  = diff(D);

	% change of bed elevation below water
	dz_ = diff(D_);

	% fraction of the perimeter that is submerged
	f = dz_./dz;

	% if the slope is zero,
	% the fraction is 1, if the segment is submerged, 0 otherwise
	fdx    = abs(dz) < sqrt(eps);
	Dm     = 0.5*(D(1:end-1,:)+D(2:end,:));
	f(fdx) = double(Dm(fdx) > 0);

	% change of width below water
	% dn_ = dn.*(dz_./dz);

	% perimeter of segment below water
	dP = f.*hypot(dz,dn);

	% total wetted perimeter of cross section
	P = sum(dP);
end

%z_s
%close all
%subplot(2,2,1)
%plot([D])
%subplot(2,2,2)
%plot([D_])
%subplot(2,2,3)
%plot(dp)
%pause

if (0)
	B = bsxfun(@plus,cvec(bottom),rvec(h));
	P = zeros(length(h),1);

	% this should better be done by linear interpolation
	for idx=1:size(B,2)
		if (isnan(h(idx)))
		P(idx,1) = NaN;
		else
		b  = B(:,idx);
		f1 = find(b>0,1,'first');
		f2 = find(b>0,1,'last');
		%bottom = [0; B(mask(:,idx),idx); 0];
		%bottom = [0; b(mask(:,idx)); 0];
		%bottom = [b(mask(:,idx))];
		%dy_dx = (1./dx)*diff(b);
		if (1 == f1)
			% extend by 45deg
			dp1 = sqrt(2)*(h(idx)+bottom(f1));
		else
			% first piece
			dh1  = -(h(idx)+bottom(f1-1));
			dx1  = dx*dh1./(bottom(f1-1)-bottom(f1));
			dp1 = sqrt(dx1*dx1 + dh1*dh1);
		end
		% last piece
		if (f2 == length(bottom))
			dp2 = sqrt(2)*(h(idx)+bottom(f2));
		else
			dh2  = -(bottom(f2+1)+h(idx));
			dx2  = dx*dh2./(bottom(f2+1)-bottom(f2));
			dp2 = sqrt(dx2*dx2 + dh2*dh2);
		end
			% wetted perimeter
			dy = diff(b(f1:f2));
			P(idx,1) = sum(sqrt(dx*dx + dy.*dy)) + dp1 + dp2;
		end
	end % for idx
end
end % csperimeter

