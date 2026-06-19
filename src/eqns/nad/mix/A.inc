
if (trans_) then
A(1, 1) = del_x*(V_g - l_i - 1)
A(2, 1) = del_x*(-V_g*alpha_gam + eta)
A(3, 1) = -alpha_grv*del_x*eta
A(4, 1) = 0
A(5, 1) = del_x*ups_T
A(6, 1) = 0
A(1, 2) = del_x*(-As_iso + om2)
A(2, 2) = del_x*(As - U - l_i + 3)
A(3, 2) = -alpha_grv*del_x*(1 - U)
A(4, 2) = alpha_grv*del_x
A(5, 2) = del_x*ups_T
A(6, 2) = 0
A(1, 3) = 0
A(2, 3) = 0
A(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A(4, 3) = alpha_grv*del_x
A(5, 3) = 0
A(6, 3) = 0
A(1, 4) = -As*U*alpha_grv*del_x
A(2, 4) = -U*V_g*alpha_grv*del_x
A(3, 4) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv*( &
      1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv*( &
      -U - l_i + 3)))
A(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
A(5, 4) = U*alpha_grv*del_x*ups_T
A(6, 4) = 0
A(1, 5) = V*del_x*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - &
      om2) - 4*nabla_ad)/f_rht
A(2, 5) = V*del_x*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/ &
      f_rht
A(3, 5) = del_x*(-V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht + V* &
      alpha_grv*nabla_ad*(1 - U)/f_rht)
A(4, 5) = -V*alpha_grv*del_x*nabla_ad/f_rht
A(5, 5) = del_x*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A(6, 5) = -V*del_x*nabla/(c_rad*f_rht)
A(1, 6) = del_x*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl* &
      c_rad*lambda*(-1 + nabla_ad/nabla))
A(2, 6) = del_x*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A(3, 6) = -alpha_grv*del_x*eta*f_conv
A(4, 6) = 0
A(5, 6) = del_x*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S - &
      alpha_hfl*c_rad*lambda/(V*nabla))
A(6, 6) = del_x*(-l_i - 1)
else
A(1, 1) = del_x*(V_g - l_i - 1)
A(2, 1) = del_x*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = -As*U*alpha_grv*del_x
A(5, 1) = V*del_x*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - &
      om2) - 4*nabla_ad)/f_rht
A(6, 1) = del_x*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl* &
      c_rad*lambda*(-1 + nabla_ad/nabla))
A(1, 2) = del_x*(-V_g*alpha_gam + eta)
A(2, 2) = del_x*(As - U - l_i + 3)
A(3, 2) = 0
A(4, 2) = -U*V_g*alpha_grv*del_x
A(5, 2) = V*del_x*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/ &
      f_rht
A(6, 2) = del_x*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A(1, 3) = -alpha_grv*del_x*eta
A(2, 3) = -alpha_grv*del_x*(1 - U)
A(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A(4, 3) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv*( &
      1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv*( &
      -U - l_i + 3)))
A(5, 3) = del_x*(-V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht + V* &
      alpha_grv*nabla_ad*(1 - U)/f_rht)
A(6, 3) = -alpha_grv*del_x*eta*f_conv
A(1, 4) = 0
A(2, 4) = alpha_grv*del_x
A(3, 4) = alpha_grv*del_x
A(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
A(5, 4) = -V*alpha_grv*del_x*nabla_ad/f_rht
A(6, 4) = 0
A(1, 5) = del_x*ups_T
A(2, 5) = del_x*ups_T
A(3, 5) = 0
A(4, 5) = U*alpha_grv*del_x*ups_T
A(5, 5) = del_x*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A(6, 5) = del_x*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S - &
      alpha_hfl*c_rad*lambda/(V*nabla))
A(1, 6) = 0
A(2, 6) = 0
A(3, 6) = 0
A(4, 6) = 0
A(5, 6) = -V*del_x*nabla/(c_rad*f_rht)
A(6, 6) = del_x*(-l_i - 1)
endif
if (PRESENT(r)) then
r(1) = del_x*lambda*y_T_1/(alpha_omg*omega_c**2)
r(2) = -c_1*del_x*l_e*y_T_1
r(3) = 0
r(4) = 0
r(5) = c_1*del_x*y_T_1*(V*eta*(-nabla + nabla_ad)/f_rht + V*l_e*nabla_ad &
      /f_rht)
r(6) = del_x*f_conv*lambda*y_T_1/(alpha_omg*omega_c**2)
endif
