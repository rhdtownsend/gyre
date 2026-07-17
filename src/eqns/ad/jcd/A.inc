
if (trans_) then
A(1, 1) = del_x*(V_g + beta_pi*(As - As_iso) - l_i - 1)
A(2, 1) = del_x*om2*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)/ &
      lambda
A(3, 1) = del_x*(-alpha_grv*eta + alpha_grv*(-alpha_gam*(V_g + beta_pi*( &
      As - As_iso)) + eta))
A(4, 1) = 0
A(1, 2) = del_x*eta*(-As_iso + om2)
A(2, 2) = del_x*eta*(-om2*(3 - U)/lambda + om2*(As - U - beta_pi*(As - &
      As_iso) - l_i + 3)/lambda)
A(3, 2) = del_x*(alpha_grv*eta*(l_i - 2) + eta*(-alpha_grv*(1 - U) + &
      alpha_grv*(As - U - beta_pi*(As - As_iso) - l_i + 3)))
A(4, 2) = 0
A(1, 3) = 0
A(2, 3) = 0
A(3, 3) = del_x*(2 - l_i)
A(4, 3) = del_x
A(1, 4) = -U*del_x*(As - beta_pi*(As - As_iso))
A(2, 4) = -U*del_x*om2*(V_g + beta_pi*(As - As_iso))/lambda
A(3, 4) = del_x*(-U*alpha_grv*(V_g + beta_pi*(As - As_iso)) - U*(-As - U &
      - V_g + 3) + lambda - (1 - U)*(-U - l_i + 2) + (U - 1)*(l_i - 2))
A(4, 4) = del_x*(-2*U - l_i + 3)
else
A(1, 1) = del_x*(V_g + beta_pi*(As - As_iso) - l_i - 1)
A(2, 1) = del_x*eta*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = -U*del_x*(As - beta_pi*(As - As_iso))
A(1, 2) = del_x*om2*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)/ &
      lambda
A(2, 2) = del_x*eta*(-om2*(3 - U)/lambda + om2*(As - U - beta_pi*(As - &
      As_iso) - l_i + 3)/lambda)
A(3, 2) = 0
A(4, 2) = -U*del_x*om2*(V_g + beta_pi*(As - As_iso))/lambda
A(1, 3) = del_x*(-alpha_grv*eta + alpha_grv*(-alpha_gam*(V_g + beta_pi*( &
      As - As_iso)) + eta))
A(2, 3) = del_x*(alpha_grv*eta*(l_i - 2) + eta*(-alpha_grv*(1 - U) + &
      alpha_grv*(As - U - beta_pi*(As - As_iso) - l_i + 3)))
A(3, 3) = del_x*(2 - l_i)
A(4, 3) = del_x*(-U*alpha_grv*(V_g + beta_pi*(As - As_iso)) - U*(-As - U &
      - V_g + 3) + lambda - (1 - U)*(-U - l_i + 2) + (U - 1)*(l_i - 2))
A(1, 4) = 0
A(2, 4) = 0
A(3, 4) = del_x
A(4, 4) = del_x*(-2*U - l_i + 3)
endif
if (PRESENT(r)) then
r(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
r(2) = -del_x*l_e*lambda*y_T_1/(alpha_omg*omega_c**2)
r(3) = 0
r(4) = 0
endif
