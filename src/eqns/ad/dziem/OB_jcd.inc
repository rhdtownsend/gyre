
if (trans_) then
B(1, 1) = -b_11 + chi
B(2, 1) = -b_12
B(3, 1) = alpha_grv*b_12 + alpha_grv*b_12*(-l_e - 1 + lambda/(alpha_omg* &
      c_1*omega_c**2))/(As + V_g)
B(4, 1) = 0
B(1, 2) = U*alpha_gbc
B(2, 2) = 0
B(3, 2) = l_e + 1
B(4, 2) = 1
else
B(1, 1) = -b_11 + chi
B(2, 1) = U*alpha_gbc
B(1, 2) = -b_12
B(2, 2) = 0
B(1, 3) = alpha_grv*b_12 + alpha_grv*b_12*(-l_e - 1 + lambda/(alpha_omg* &
      c_1*omega_c**2))/(As + V_g)
B(2, 3) = l_e + 1
B(1, 4) = 0
B(2, 4) = 1
endif
if (PRESENT(s)) then
s(1) = -b_12*c_1*y_T_1*(-l_e - 1 + lambda/(alpha_omg*c_1*omega_c**2))/( &
      As + V_g)
s(2) = 0
endif
