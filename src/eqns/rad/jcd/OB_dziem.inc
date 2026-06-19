
if (trans_) then
B(1, 1) = 1 - (alpha_omg*c_1*omega_c**2 + 4)/V
B(2, 1) = -alpha_omg*c_1*omega_c**2
else
B(1, 1) = 1 - (alpha_omg*c_1*omega_c**2 + 4)/V
B(1, 2) = -alpha_omg*c_1*omega_c**2
endif
if (PRESENT(s)) then
s(1) = 0
endif
