B(1, 1) = 1 + (-alpha_omg*omega_c**2*c_1 - 4 + lambda/(alpha_omg*omega_c &
      **2*c_1))/V
B(2, 1) = alpha_gbc*alpha_grv*U
B(3, 1) = -4*V*nabla + 2
B(1, 2) = -alpha_omg*omega_c**2*c_1
B(2, 2) = 0
B(3, 2) = 4*V*alpha_omg*nabla*omega_c**2*c_1
B(1, 3) = -1 - alpha_grv*(-l_e - 1 + lambda/(alpha_omg*omega_c**2*c_1))/ &
      V
B(2, 3) = -alpha_grv*l_e + alpha_grv*(1 - U) - 1
B(3, 3) = 4*V*nabla
B(1, 4) = 0
B(2, 4) = -alpha_grv
B(3, 4) = 0
B(1, 5) = 0
B(2, 5) = 0
B(3, 5) = 4*f_rh
B(1, 6) = 0
B(2, 6) = 0
B(3, 6) = -1
