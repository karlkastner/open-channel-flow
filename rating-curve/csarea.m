% 2015-03-16 18:14:47.261958638 +0100
% Karl Kastner, Berlin
%
%% predict cross sectional area from transverse bed level profile
%% and surface elevation
% function A = csarea(zs,zb,dw)
function A = csarea(zs,zb,dw)
	if (isvector(zb))
		zb = cvec(zb);
	end

	% depth across section
	D = bsxfun(@minus,rvec(zs),zb);
%	D = bsxfun(@plus,cvec(bottom),rvec(h));

	% zero parts of cross section that are dry
	D = max(D,0);

	% integrate area
	A = dw*sum(D);

	% make a column vector
	%A = vec(A);

%    % area
%    A = zeros(length(h),1);
%    % for each water level
%    for idx=1:length(h)
%        if (isnan(h(idx)))
%            A(idx,1) = NaN;
%        else
%	    % depth
%            d = D(:,idx);
%	    % part of
%            fdx = d > 0;
%            A(idx,1) = sum(d(fdx))*dw;
%    %        afunc = @(h) cvec(nansum(bsxfun(@plus,cvec(double(calib.bottom)),rvec(h))));
%        end % end
%    end % for idx
end % csarea

