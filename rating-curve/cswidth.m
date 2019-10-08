% 2015-03-16 19:16:15.139306769 +0100
% Karl Kastner, Berlin
%
%% determine cross section width
%%
% z_s : water surface level
% z_b : bed level
%
% h      : d-reference (larger value is larger depth)
% bottom : d-reference (larger value is larger depth)
function W = cswidth(z_s, z_b, dn)
%	bottom = [-2*max(h); -max(h); bottom; -max(h); -2*max(h)];
%	% exten boundary by 45deg element
%	bottom = [bottom(1)-dx; bottom; bottom(end)-dx];

	% extend banks by rectangular walls
	% TODO extend slanting
	max_z = max(z_s);
	z_b = [max_z;
               z_b;
               max_z];

	% depth
	D = bsxfun(@minus,rvec(z_s),cvec(z_b));

	% depth only below water
	D_ = max(0,D);

	% change of depth
	dz  = diff(D);
	% change of vertical only below water
	dz_ = diff(D_);

	% change of width
	 %dn  = dn*ones(size(dz));
	% change of width below water
	dn_ = dn.*(dz_./dz);

	% perimeter per segment
	%dp = hypot(dz_,dn_);
	% perimeter per cross section
	%P = sum(dp)';
	W = sum(dn_);

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

