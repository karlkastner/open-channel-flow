% Thu 17 May 13:29:10 CEST 2018
%% cut the boundary from the domain
%% wa : width of inlet to side channel
%% wb : width of side channel
function rdx = cut_boundary(x,y,Lx,Ly,yc,wa,wb,dx)
	   %rdx = (x>=Lx-wb) & (y>=yc+wa);
	   %rdx(:) = 0;
	
	% smooth downstream end of inlet
	   d   = (y-(yc+wa/2+0*pi*wb/2))/(wb*pi);
	   % remove piece on top of right channel
	   rdx = (x>=Lx-wb) & d > 1/4;
	   rdx = rdx ...
                 | ( (d >= -1/4) & (d <= 1/4) ...
	             & (Lx-x < 1/2*wb*(1+sin(2*pi*d))) );

%	   d   = (y-yc)/wa;
%	   rdx = rdx ...
%                 | ( (d >= 0) & (d <= 1) ...
%	             & (Lx-x < 1/2*wb*(1-cos(pi*d))) );
	% cut thin dam
	   d   = x-(Lx-wb);
	   rdx = rdx | ...
		 ( (y<=yc-wa/2) & (d>=0) & (d<dx) );
%		figure(2)
%		imagesc(d)
%		pause
end

