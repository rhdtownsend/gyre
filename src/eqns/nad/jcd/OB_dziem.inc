
if (trans_) then
B(1, 1) = 1 + (-alpha_omg*c_1*omega_c**2 - 4 + lambda/(alpha_omg*c_1* &
      omega_c**2))/V
B(2, 1) = -alpha_omg*c_1*omega_c**2/lambda
B(3, 1) = -1 - alpha_grv*(-l_e - 1 + lambda/(alpha_omg*c_1*omega_c**2))/ &
      V
B(4, 1) = 0
B(5, 1) = 0
B(6, 1) = 0
B(1, 2) = U*alpha_gbc*alpha_grv
B(2, 2) = 0
B(3, 2) = -alpha_grv*l_e + alpha_grv*(1 - U) - 1
B(4, 2) = -alpha_grv
B(5, 2) = 0
B(6, 2) = 0
B(1, 3) = -4*V*nabla_ad + 2
B(2, 3) = 4*V*alpha_omg*c_1*nabla_ad*omega_c**2/lambda
B(3, 3) = 4*V*nabla_ad
B(4, 3) = 0
B(5, 3) = 4*f_rht
B(6, 3) = -1
else
B(1, 1) = 1 + (-alpha_omg*c_1*omega_c**2 - 4 + lambda/(alpha_omg*c_1* &
      omega_c**2))/V
B(2, 1) = U*alpha_gbc*alpha_grv
B(3, 1) = -4*V*nabla_ad + 2
B(1, 2) = -alpha_omg*c_1*omega_c**2/lambda
B(2, 2) = 0
B(3, 2) = 4*V*alpha_omg*c_1*nabla_ad*omega_c**2/lambda
B(1, 3) = -1 - alpha_grv*(-l_e - 1 + lambda/(alpha_omg*c_1*omega_c**2))/ &
      V
B(2, 3) = -alpha_grv*l_e + alpha_grv*(1 - U) - 1
B(3, 3) = 4*V*nabla_ad
B(1, 4) = 0
B(2, 4) = -alpha_grv
B(3, 4) = 0
B(1, 5) = 0
B(2, 5) = 0
B(3, 5) = 4*f_rht
B(1, 6) = 0
B(2, 6) = 0
B(3, 6) = -1
endif
if (PRESENT(s)) then
s(1) = -c_1*y_T_1*(-l_e - 1 + lambda/(alpha_omg*c_1*omega_c**2))/V
s(2) = 0
s(3) = 0
endif
