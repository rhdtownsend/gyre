
if (trans_) then
B(1, 1) = a_11 - chi
B(2, 1) = a_12*alpha_omg*c_1*omega_c**2
else
B(1, 1) = a_11 - chi
B(1, 2) = a_12*alpha_omg*c_1*omega_c**2
endif
if (PRESENT(s)) then
s(1) = 0
endif
