
if (trans_) then
A(1, 1) = del_x*(V_g - l_i - 1)
A(2, 1) = del_x*(-V_g*alpha_gam + eta)
A(3, 1) = -alpha_grv*del_x*eta
A(4, 1) = 0
A(1, 2) = del_x*(-As_iso + om2)
A(2, 2) = del_x*(As - U - l_i + 3)
A(3, 2) = -alpha_grv*del_x*(1 - U)
A(4, 2) = alpha_grv*del_x
A(1, 3) = 0
A(2, 3) = 0
A(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A(4, 3) = alpha_grv*del_x
A(1, 4) = -As*U*alpha_grv*del_x
A(2, 4) = -U*V_g*alpha_grv*del_x
A(3, 4) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv*( &
      1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv*( &
      -U - l_i + 3)))
A(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
else
A(1, 1) = del_x*(V_g - l_i - 1)
A(2, 1) = del_x*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = -As*U*alpha_grv*del_x
A(1, 2) = del_x*(-V_g*alpha_gam + eta)
A(2, 2) = del_x*(As - U - l_i + 3)
A(3, 2) = 0
A(4, 2) = -U*V_g*alpha_grv*del_x
A(1, 3) = -alpha_grv*del_x*eta
A(2, 3) = -alpha_grv*del_x*(1 - U)
A(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A(4, 3) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv*( &
      1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv*( &
      -U - l_i + 3)))
A(1, 4) = 0
A(2, 4) = alpha_grv*del_x
A(3, 4) = alpha_grv*del_x
A(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
endif
if (PRESENT(r)) then
r(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
r(2) = -c_1*del_x*l_e*y_T_1
r(3) = 0
r(4) = 0
endif
