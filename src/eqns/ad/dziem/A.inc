
if (trans_) then
A(1, 1) = del_x*(V_g + beta_pi*(As - As_iso) - l_i - 1)
A(2, 1) = del_x*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)
A(3, 1) = del_x*(alpha_grv*eta - alpha_grv*(-alpha_gam*(V_g + beta_pi*( &
      As - As_iso)) + eta))
A(4, 1) = 0
A(1, 2) = del_x*(-As_iso + om2)
A(2, 2) = del_x*(As - U - beta_pi*(As - As_iso) - l_i + 3)
A(3, 2) = del_x*(alpha_grv*(-U - l_i + 3) - alpha_grv*(As - U - beta_pi* &
      (As - As_iso) - l_i + 3))
A(4, 2) = 0
A(1, 3) = 0
A(2, 3) = 0
A(3, 3) = del_x*(-U - l_i + 3)
A(4, 3) = del_x
A(1, 4) = U*del_x*(As - beta_pi*(As - As_iso))
A(2, 4) = U*del_x*(V_g + beta_pi*(As - As_iso))
A(3, 4) = del_x*(-U*alpha_grv*(V_g + beta_pi*(As - As_iso)) + lambda)
A(4, 4) = del_x*(-U - l_i + 2)
else
A(1, 1) = del_x*(V_g + beta_pi*(As - As_iso) - l_i - 1)
A(2, 1) = del_x*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = U*del_x*(As - beta_pi*(As - As_iso))
A(1, 2) = del_x*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)
A(2, 2) = del_x*(As - U - beta_pi*(As - As_iso) - l_i + 3)
A(3, 2) = 0
A(4, 2) = U*del_x*(V_g + beta_pi*(As - As_iso))
A(1, 3) = del_x*(alpha_grv*eta - alpha_grv*(-alpha_gam*(V_g + beta_pi*( &
      As - As_iso)) + eta))
A(2, 3) = del_x*(alpha_grv*(-U - l_i + 3) - alpha_grv*(As - U - beta_pi* &
      (As - As_iso) - l_i + 3))
A(3, 3) = del_x*(-U - l_i + 3)
A(4, 3) = del_x*(-U*alpha_grv*(V_g + beta_pi*(As - As_iso)) + lambda)
A(1, 4) = 0
A(2, 4) = 0
A(3, 4) = del_x
A(4, 4) = del_x*(-U - l_i + 2)
endif
if (PRESENT(r)) then
r(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
r(2) = -c_1*del_x*l_e*y_T_1
r(3) = 0
r(4) = 0
endif
