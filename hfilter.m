% 2020-01-04 23:15:58.059737635 +0800
% quasi hydrograph by pt2-filtered rectangular 'rainfall event' (or wet season)
function y = hfilter(mag,offset,width,centre,rr,rl,n,ys)
	mag = mag*1e4;
	offset = offset*1e3;
	width = width*100;
	centre = centre*100;
%	rr = 2*atan(rr)/pi;
%	rl = 2*atan(rl)/pi;
	order = length(rr);
% w=200; r=0.99; r2=0.985; x = 0*linspace(0,1,4e3)'; for i=1:3; x([length(x)*i/4+(-w:w)])=1; end; plot([x,1e1*filter(1-r2,[1,-r2],filter(1-r2,[1,-r2],flipud(filter(1-r,[1,-r],filter(1-r,[1,-r],x)))))]); xlim([1500,2500])
	% rr and rl constants can be split in 2 (or 3)
	x = (1:n)';
	y = zeros(3*n,1);
	% TODO, do this by interpolation
	for idx=1:3
	x_ = linspace(-10,10,1e3);
	y_ = ones(size(x_));
	y_(x_>=-0.5 & x_<=0.5) = mag;
	x_ = x_*width+centre;
	y((idx-1)*n+(1:n)) = interp1(x_,y_,x,'linear');
%[width,centre]
%	y((idx-1)*n+round((centre-width/2:centre+width/2))) = mag;
	end
	y = y+offset;
	% filter right
	if (0)
	for idx=1:order
		y = filter(1-rr(idx),[1,-rr(idx)],y);
	end
	% filter left
	y = flipud(y);
	for idx=1:order
		y = filter(1-rl(idx),[1,-rl(idx)],y);
	end
	y = flipud(y);
	end
	y = filter(0.1,[1,-rr],y);
	y = flipud(y);
	y = filter(0.1,[1,-rl],y);
	y = flipud(y);
%[[mag/1e4,offset/1e3] [width,centre]/100, [rr, rl]]
	y = y(n+1:2*n);
%mag
%width
%centre
%offset
%rr
%rl
%	if (abs(rr)>=1 || abs(rl)>=1)
%		y = y*NaN;
%	end
%figure(1e3)
%clf
%plot(x,[y,ys])
%drawnow
%pause(0.1)
end


