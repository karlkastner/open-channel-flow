% Tue 19 May 16:13:02 +08 2020
function cs_div_cin = sediment_division_meijer_ksiazek_1(Qs_div_Qin)
	if (1)
		syms ci cs Qi Qs k 
		Si = ci*Qi;
		Ss = cs*Qs;
		solve(Ss/(Si-Ss) == k*Qs/(Qi - Qs),cs)
	end
	k = 2.63;
	cs_div_cin = k./(1 + Qs_div_Qin*(k-1));
end

