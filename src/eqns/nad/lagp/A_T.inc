A_T(1, 1) = invx*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A_T(2, 1) = invx*(-V_g*alpha_gam + eta)/V_2
A_T(3, 1) = alpha_grv*eta*invx
A_T(4, 1) = 0
A_T(5, 1) = invx*ups_T
A_T(6, 1) = 0
A_T(1, 2) = invx*(-V_2*(-V_g*alpha_gam + V_g + eta - l_i - 1) + V_2*(As &
      - As_iso - U - l_i + om2 + 3))
A_T(2, 2) = invx*(V_2*(dV_2/V_2 + (As - U - l_i + 3)/V_2) + V_g* &
      alpha_gam - eta)
A_T(3, 2) = -V_2*alpha_grv*eta*invx
A_T(4, 2) = -V_2*alpha_grv*invx
A_T(5, 2) = 0
A_T(6, 2) = 0
A_T(1, 3) = 0
A_T(2, 3) = 0
A_T(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A_T(4, 3) = alpha_grv*invx
A_T(5, 3) = 0
A_T(6, 3) = 0
A_T(1, 4) = invx*(As*U*alpha_grv + U*V_g*alpha_grv)
A_T(2, 4) = U*V_g*alpha_grv*invx/V_2
A_T(3, 4) = alpha_grv*invx*lambda
A_T(4, 4) = alpha_grv*invx*(-U - l_i + 2)
A_T(5, 4) = -U*alpha_grv*invx*ups_T
A_T(6, 4) = 0
A_T(1, 5) = invx*(V*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad &
      ))/f_rht + V*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - &
      om2) - 4*nabla_ad)/f_rht)
A_T(2, 5) = V*invx*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad)) &
      /(V_2*f_rht)
A_T(3, 5) = V*alpha_grv*eta*invx*(-nabla + nabla_ad)/f_rht
A_T(4, 5) = V*alpha_grv*invx*nabla_ad/f_rht
A_T(5, 5) = invx*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A_T(6, 5) = -V*invx*nabla/(c_rad*f_rht)
A_T(1, 6) = invx*(alpha_hfl*c_rad*lambda*(-1 + nabla_ad/nabla) - &
      alpha_hfl*c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A_T(2, 6) = invx*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)/V_2
A_T(3, 6) = alpha_grv*eta*f_conv*invx
A_T(4, 6) = 0
A_T(5, 6) = invx*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S &
      - alpha_hfl*c_rad*lambda/(V*nabla))
A_T(6, 6) = invx*(-l_i - 1)
