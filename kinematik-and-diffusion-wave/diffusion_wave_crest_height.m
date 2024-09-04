% Sun 23 May 12:00:51 CEST 2021
% TODO qmax, qmin
% crest height of a sinusoidal diffusion wave in a rectangular channel of arbitrary width
function [h,dh_dt] = diffusion_wave_crest_height(h0,C,S0,T,t)
	h = h0./(1+(sqrt(h0./S0.^3).*t.*pi.^2)./(2*T.^2.*C)).^2;
	dh_dt = -((h./S0).^(3./2).*1./T.^2.*pi.^2)./C;
end

