A_t(1, 1) = del_x*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A_t(2, 1) = del_x*(-V_g*alpha_gam + eta)/V_2
A_t(3, 1) = alpha_grv*del_x*eta
A_t(4, 1) = 0
A_t(1, 2) = del_x*(-V_2*(-V_g*alpha_gam + V_g + eta - l_i - 1) + V_2*(As &
      - As_iso - U - l_i + om2 + 3))
A_t(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U - l_i + 3)/V_2) + V_g* &
      alpha_gam - eta)
A_t(3, 2) = -V_2*alpha_grv*del_x*eta
A_t(4, 2) = -V_2*alpha_grv*del_x
A_t(1, 3) = 0
A_t(2, 3) = 0
A_t(3, 3) = alpha_grv*del_x*(-U - l_i + 3)
A_t(4, 3) = alpha_grv*del_x
A_t(1, 4) = del_x*(As*U*alpha_grv + U*V_g*alpha_grv)
A_t(2, 4) = U*V_g*alpha_grv*del_x/V_2
A_t(3, 4) = alpha_grv*del_x*lambda
A_t(4, 4) = alpha_grv*del_x*(-U - l_i + 2)
if (PRESENT(f)) then
f(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
f(2) = del_x*(-V_2*c_1*l_e*y_T_1 - V_2*lambda*y_T_1/(alpha_omg*omega_c** &
      2))
f(3) = 0
f(4) = 0
endif
