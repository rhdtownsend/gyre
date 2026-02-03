A(1, 1) = invx*(V_g - l_i - 1)
A(2, 1) = invx*(-As_iso + om2)
A(3, 1) = 0
A(4, 1) = As*U*alpha_grv*invx
A(5, 1) = V*invx*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - om2 &
      ) - 4*nabla_ad)/f_rht
A(6, 1) = invx*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl* &
      c_rad*lambda*(-1 + nabla_ad/nabla))
A(1, 2) = invx*(-V_g*alpha_gam + eta)
A(2, 2) = invx*(As - U - l_i + 3)
A(3, 2) = 0
A(4, 2) = U*V_g*alpha_grv*invx
A(5, 2) = V*invx*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/ &
      f_rht
A(6, 2) = invx*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A(1, 3) = invx*(V_g*alpha_gam + alpha_grv*eta - eta)
A(2, 3) = invx*(-As + U + alpha_grv*(-U - l_i + 3) + l_i - 3)
A(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A(4, 3) = invx*(-U*V_g*alpha_grv + alpha_grv*lambda)
A(5, 3) = invx*(V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht - V*(-V* &
      c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/f_rht)
A(6, 3) = invx*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_grv*eta &
      *f_conv + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)
A(1, 4) = 0
A(2, 4) = 0
A(3, 4) = alpha_grv*invx
A(4, 4) = alpha_grv*invx*(-U - l_i + 2)
A(5, 4) = V*alpha_grv*invx*nabla_ad/f_rht
A(6, 4) = 0
A(1, 5) = invx*ups_T
A(2, 5) = invx*ups_T
A(3, 5) = 0
A(4, 5) = -U*alpha_grv*invx*ups_T
A(5, 5) = invx*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A(6, 5) = invx*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S - &
      alpha_hfl*c_rad*lambda/(V*nabla))
A(1, 6) = 0
A(2, 6) = 0
A(3, 6) = 0
A(4, 6) = 0
A(5, 6) = -V*invx*nabla/(c_rad*f_rht)
A(6, 6) = invx*(-l_i - 1)
