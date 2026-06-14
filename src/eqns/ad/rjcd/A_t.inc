A_t(1, 1) = del_x*(V_g - l_i - 1)
A_t(2, 1) = del_x*om2*(-V_g*alpha_gam + eta)
A_t(3, 1) = del_x*(-V_g*alpha_gam - alpha_grv*eta + eta)
A_t(4, 1) = 0
A_t(1, 2) = del_x*(-As_iso + om2)/om2
A_t(2, 2) = del_x*(-om2*(3 - U) + om2*(As - U - l_i + 3))/om2
A_t(3, 2) = del_x*((alpha_grv*(1 - U) - alpha_grv*(-U - l_i + 3))/om2 + &
      (As - U - alpha_grv*(1 - U) - l_i + 3)/om2)
A_t(4, 2) = 0
A_t(1, 3) = 0
A_t(2, 3) = 0
A_t(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A_t(4, 3) = alpha_grv*del_x
A_t(1, 4) = -As*U*alpha_grv*del_x
A_t(2, 4) = -U*V_g*alpha_grv*del_x*om2
A_t(3, 4) = del_x*(-U*V_g*alpha_grv - U*(-As - U - V_g + 3) + alpha_grv* &
      lambda - alpha_grv*(1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 &
      - U) - alpha_grv*(-U - l_i + 3)))
A_t(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
if (PRESENT(f)) then
f(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
f(2) = -del_x*l_e*y_T_1/(alpha_omg*omega_c**2)
f(3) = 0
f(4) = 0
endif
