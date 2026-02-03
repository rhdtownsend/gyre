A(1, 1) = invx*(-V_g*alpha_gam + V_g + eta - l_i - 1)
A(2, 1) = invx*(-V_2*(-V_g*alpha_gam + V_g + eta - l_i - 1) + V_2*(As - &
      As_iso - U - l_i + om2 + 3))
A(3, 1) = 0
A(4, 1) = invx*(As*U*alpha_grv + U*V_g*alpha_grv)
A(5, 1) = invx*(V*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/ &
      f_rht + V*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - om2 &
      ) - 4*nabla_ad)/f_rht)
A(6, 1) = invx*(alpha_hfl*c_rad*lambda*(-1 + nabla_ad/nabla) - alpha_hfl &
      *c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A(1, 2) = invx*(-V_g*alpha_gam + eta)/V_2
A(2, 2) = invx*(V_2*(dV_2/V_2 + (As - U - l_i + 3)/V_2) + V_g*alpha_gam &
      - eta)
A(3, 2) = 0
A(4, 2) = U*V_g*alpha_grv*invx/V_2
A(5, 2) = V*invx*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/( &
      V_2*f_rht)
A(6, 2) = invx*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)/V_2
A(1, 3) = alpha_grv*eta*invx
A(2, 3) = -V_2*alpha_grv*eta*invx
A(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A(4, 3) = alpha_grv*invx*lambda
A(5, 3) = V*alpha_grv*eta*invx*(-nabla + nabla_ad)/f_rht
A(6, 3) = alpha_grv*eta*f_conv*invx
A(1, 4) = 0
A(2, 4) = -V_2*alpha_grv*invx
A(3, 4) = alpha_grv*invx
A(4, 4) = alpha_grv*invx*(-U - l_i + 2)
A(5, 4) = V*alpha_grv*invx*nabla_ad/f_rht
A(6, 4) = 0
A(1, 5) = invx*ups_T
A(2, 5) = 0
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
