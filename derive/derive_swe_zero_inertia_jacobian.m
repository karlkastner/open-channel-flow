% Fri  4 Oct 14:57:13 CEST 2024
	syms z zl  bl br zr
	syms C h hl hr dx positive

	Sl = (h +z  - hl-zl)/dx;
	Sr = (hr+zr - h -z )/dx;

	% velocity, this is the coefficient
	ul = -C*sign(Sl)*sqrt(0.5*(hl + h ).*abs(Sl)) 
	ur = -C*sign(Sr)*sqrt(0.5*(h  + hr).*abs(Sr)) 

	dh_dt = -(     0.5*ur.*( (1-br).*h   + (1+br).*hr) ...
	            - 0.5*ul.*( (1-bl).*hl  + (1+bl).*h) ...
		 )/dx;

Jl = diff(dh_dt,hl)
Jr = diff(dh_dt,hr)
Jc = diff(dh_dt,h)

		

	

