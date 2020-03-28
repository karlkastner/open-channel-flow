% Fri 14 Jul 14:33:04 CEST 2017
%
%% linear wave routing (linearised kinematic wave)
%
function x = linear_wave(c,dt,D,x)
	% make coefficients time step invariant
	% c(2) = (c(2) + (1/((1+c(3)))-1))/dt;
	% c3 is interpreted as diffusion coefficient (negative sign)

	

if (0)
	c(1) = c(1)*(1-c(3));
	c(2) = c(2)/dt - (1/((1-c(3))));
%	c(2) = c(2)/dt + 1/(1+c(3));
	% c(3) = sign(c(3))*abs(c(3))^dt;

	% y_i = c_3/dt*y_i-1 + c_1*(1-c_3/dt)*x(t/dt-c_2)
	% TODO use zi
	% TODO simply fill in zeros before c(3) instead of using mpoweri

	% scale and diffuse
	% needs to skip startup time
	a = [1 -c(3)];
	b = [c(1)];
	x = filter(b,a,x);
else
	%
	% zero-phase filtering
	%
	c10 = c(1);

	[b a] = digital_low_pass_filter(c(3),dt);

	% amplification (filter is run twice, so take square root)
	% c(1) = c(1)^(0.5)*(1-c(3));
	b(1) = c(1)^0.5*b(1);

	% time shift
	c(2) = c(2)/dt;

	if (1)
		% initial value
		zi = c10.^0.5*x(1);
		% filter fw in time
		x = filter(b,a,x,zi);
		% filter bw in time
		% initial value
		zi = c10.^0.5*x(end);
		x = flipud(filter(b,a,flipud(x)));
	else
		x = filtfilt(b,a,x);
	end
	
end

	% shift in time
	switch (2)
	case {0}
		x = circshift_fractional(x,-c(2));
	case {1}
		if (c(2) > 0)
		x = mpoweri(D{1},c(2),x);
		else
		x = mpoweri(D{2},-c(2),x);
		end
	case {2}
		n = length(x);
		dn_max = 3/24*1/dt;	
		% c2 is interpreted as delay, so at out(i) is in(i-c2)
		x = interp1_limited((1:n),x,(1:n)-c(2),dn_max,'linear');
	end
end

