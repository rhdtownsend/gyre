
if (trans_) then
B(1, 1) = 0
B(2, 1) = alpha_omg*c_1*omega_c**2/lambda
B(3, 1) = 1 - alpha_grv
B(4, 1) = 0
B(1, 2) = 0
B(2, 2) = 0
B(3, 2) = -2
B(4, 2) = -1
else
B(1, 1) = 0
B(2, 1) = 0
B(1, 2) = alpha_omg*c_1*omega_c**2/lambda
B(2, 2) = 0
B(1, 3) = 1 - alpha_grv
B(2, 3) = -2
B(1, 4) = 0
B(2, 4) = -1
endif
if (PRESENT(s)) then
s(1) = -c_1*y_T_1
s(2) = 0
endif
