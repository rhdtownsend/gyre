
if (trans_) then
A(1, 1) = del_x*(V_g - alpha_gam*(V_g + beta_pi*(As - As_iso)) + beta_pi &
      *(As - As_iso) + eta - l_i - 1)
A(2, 1) = del_x*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)/V_2
A(3, 1) = alpha_grv*del_x*eta
A(4, 1) = 0
A(1, 2) = del_x*(-V_2*(V_g - alpha_gam*(V_g + beta_pi*(As - As_iso)) + &
      beta_pi*(As - As_iso) + eta - l_i - 1) + V_2*(As - As_iso - U - &
      beta_pi*(As - As_iso) - l_i + om2 + 3))
A(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U - beta_pi*(As - As_iso) - l_i + &
      3)/V_2) + alpha_gam*(V_g + beta_pi*(As - As_iso)) - eta)
A(3, 2) = -V_2*alpha_grv*del_x*eta
A(4, 2) = -V_2*alpha_grv*del_x
A(1, 3) = 0
A(2, 3) = 0
A(3, 3) = del_x*(-U - l_i + 3)
A(4, 3) = del_x
A(1, 4) = del_x*(U*(As - beta_pi*(As - As_iso)) + U*(V_g + beta_pi*(As - &
      As_iso)))
A(2, 4) = U*del_x*(V_g + beta_pi*(As - As_iso))/V_2
A(3, 4) = del_x*lambda
A(4, 4) = del_x*(-U - l_i + 2)
else
A(1, 1) = del_x*(V_g - alpha_gam*(V_g + beta_pi*(As - As_iso)) + beta_pi &
      *(As - As_iso) + eta - l_i - 1)
A(2, 1) = del_x*(-V_2*(V_g - alpha_gam*(V_g + beta_pi*(As - As_iso)) + &
      beta_pi*(As - As_iso) + eta - l_i - 1) + V_2*(As - As_iso - U - &
      beta_pi*(As - As_iso) - l_i + om2 + 3))
A(3, 1) = 0
A(4, 1) = del_x*(U*(As - beta_pi*(As - As_iso)) + U*(V_g + beta_pi*(As - &
      As_iso)))
A(1, 2) = del_x*(-alpha_gam*(V_g + beta_pi*(As - As_iso)) + eta)/V_2
A(2, 2) = del_x*(V_2*(dV_2/V_2 + (As - U - beta_pi*(As - As_iso) - l_i + &
      3)/V_2) + alpha_gam*(V_g + beta_pi*(As - As_iso)) - eta)
A(3, 2) = 0
A(4, 2) = U*del_x*(V_g + beta_pi*(As - As_iso))/V_2
A(1, 3) = alpha_grv*del_x*eta
A(2, 3) = -V_2*alpha_grv*del_x*eta
A(3, 3) = del_x*(-U - l_i + 3)
A(4, 3) = del_x*lambda
A(1, 4) = 0
A(2, 4) = -V_2*alpha_grv*del_x
A(3, 4) = del_x
A(4, 4) = del_x*(-U - l_i + 2)
endif
if (PRESENT(r)) then
r(1) = del_x*eta*z
r(2) = del_x*(-V_2*dz - V_2*eta*z)
r(3) = 0
r(4) = 0
endif
