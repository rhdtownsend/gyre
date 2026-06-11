        A_t(1, 1) = del_x*(V_g - l_i - 1)
A_t(2, 1) = del_x*(-V_g*alpha_gam + eta)
A_t(3, 1) = -alpha_grv*del_x*eta
A_t(4, 1) = 0
A_t(5, 1) = del_x*ups_T
A_t(6, 1) = 0
A_t(1, 2) = del_x*(-As_iso + om2)
A_t(2, 2) = del_x*(As - U - l_i + 3)
A_t(3, 2) = -alpha_grv*del_x*(1 - U)
A_t(4, 2) = alpha_grv*del_x
A_t(5, 2) = del_x*ups_T
A_t(6, 2) = 0
A_t(1, 3) = 0
A_t(2, 3) = 0
A_t(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A_t(4, 3) = alpha_grv*del_x
A_t(5, 3) = 0
A_t(6, 3) = 0
A_t(1, 4) = -As*U*alpha_grv*del_x
A_t(2, 4) = -U*V_g*alpha_grv*del_x
A_t(3, 4) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv &
      *(1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv* &
      (-U - l_i + 3)))
A_t(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
A_t(5, 4) = U*alpha_grv*del_x*ups_T
A_t(6, 4) = 0
A_t(1, 5) = V*del_x*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - &
      om2) - 4*nabla_ad)/f_rht
A_t(2, 5) = V*del_x*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad &
      ))/f_rht
A_t(3, 5) = del_x*(-V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht + V* &
      alpha_grv*nabla_ad*(1 - U)/f_rht)
A_t(4, 5) = -V*alpha_grv*del_x*nabla_ad/f_rht
A_t(5, 5) = del_x*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2 &
      )
A_t(6, 5) = -V*del_x*nabla/(c_rad*f_rht)
A_t(1, 6) = del_x*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl* &
      c_rad*lambda*(-1 + nabla_ad/nabla))
A_t(2, 6) = del_x*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A_t(3, 6) = -alpha_grv*del_x*eta*f_conv
A_t(4, 6) = 0
A_t(5, 6) = del_x*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S &
      - alpha_hfl*c_rad*lambda/(V*nabla))
A_t(6, 6) = del_x*(-l_i - 1)
        if (PRESENT(f)) then
        A_t(1, 1) = del_x*(V_g - l_i - 1)
A_t(2, 1) = del_x*(-V_g*alpha_gam + eta)
A_t(3, 1) = -alpha_grv*del_x*eta
A_t(4, 1) = 0
A_t(5, 1) = del_x*ups_T
A_t(6, 1) = 0
A_t(1, 2) = del_x*(-As_iso + om2)
A_t(2, 2) = del_x*(As - U - l_i + 3)
A_t(3, 2) = -alpha_grv*del_x*(1 - U)
A_t(4, 2) = alpha_grv*del_x
A_t(5, 2) = del_x*ups_T
A_t(6, 2) = 0
A_t(1, 3) = 0
A_t(2, 3) = 0
A_t(3, 3) = del_x*(-alpha_grv*(1 - U) + alpha_grv*(-U - l_i + 3))
A_t(4, 3) = alpha_grv*del_x
A_t(5, 3) = 0
A_t(6, 3) = 0
A_t(1, 4) = -As*U*alpha_grv*del_x
A_t(2, 4) = -U*V_g*alpha_grv*del_x
A_t(3, 4) = del_x*(-U*(-As - U - V_g + 3) + alpha_grv*lambda - alpha_grv &
      *(1 - U)*(-U - l_i + 2) + (U - 1)*(alpha_grv*(1 - U) - alpha_grv* &
      (-U - l_i + 3)))
A_t(4, 4) = del_x*(-alpha_grv*(U - 1) + alpha_grv*(-U - l_i + 2))
A_t(5, 4) = U*alpha_grv*del_x*ups_T
A_t(6, 4) = 0
A_t(1, 5) = V*del_x*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - &
      om2) - 4*nabla_ad)/f_rht
A_t(2, 5) = V*del_x*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad &
      ))/f_rht
A_t(3, 5) = del_x*(-V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht + V* &
      alpha_grv*nabla_ad*(1 - U)/f_rht)
A_t(4, 5) = -V*alpha_grv*del_x*nabla_ad/f_rht
A_t(5, 5) = del_x*(V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2 &
      )
A_t(6, 5) = -V*del_x*nabla/(c_rad*f_rht)
A_t(1, 6) = del_x*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl* &
      c_rad*lambda*(-1 + nabla_ad/nabla))
A_t(2, 6) = del_x*(V*alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl* &
      c_rad*lambda*nabla_ad/nabla + eta*f_conv)
A_t(3, 6) = -alpha_grv*del_x*eta*f_conv
A_t(4, 6) = 0
A_t(5, 6) = del_x*(alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S &
      - alpha_hfl*c_rad*lambda/(V*nabla))
A_t(6, 6) = del_x*(-l_i - 1)
        endif

