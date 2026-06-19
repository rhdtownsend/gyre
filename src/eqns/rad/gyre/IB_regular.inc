
if (trans_) then
B(1, 1) = alpha_omg*c_1*omega_c**2
B(2, 1) = 0
else
B(1, 1) = alpha_omg*c_1*omega_c**2
B(1, 2) = 0
endif
if (PRESENT(s)) then
s(1) = 0
endif
