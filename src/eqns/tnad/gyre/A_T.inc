A_T(1, 1) = invx*(F_trb*(V_g - 2)*(V_g*alpha_gam - eta)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + V_g - l_i - 1)
A_T(2, 1) = invx*(F_trb*(-V_g*alpha_gam + eta)*(V_g*alpha_gam - eta)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) - V_g*alpha_gam + eta)
A_T(3, 1) = invx*(F_trb*alpha_grv*eta*(V_g*alpha_gam - eta)/(F_trb*(-V_g &
      *alpha_gam + eta) + 1) + alpha_grv*eta)
A_T(4, 1) = 0
A_T(5, 1) = invx*(F_trb*ups_T*(V_g*alpha_gam - eta)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + ups_T)
A_T(6, 1) = 0
A_T(1, 2) = invx*(-As_iso + F_trb*(V_g - 2)**2/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) + om2)
A_T(2, 2) = invx*(As + F_trb*(V_g - 2)*(-V_g*alpha_gam + eta)/(F_trb*( &
      -V_g*alpha_gam + eta) + 1) - U - l_i + 3)
A_T(3, 2) = F_trb*alpha_grv*eta*invx*(V_g - 2)/(F_trb*(-V_g*alpha_gam + &
      eta) + 1)
A_T(4, 2) = -alpha_grv*invx
A_T(5, 2) = invx*(F_trb*ups_T*(V_g - 2)/(F_trb*(-V_g*alpha_gam + eta) + &
      1) + ups_T)
A_T(6, 2) = 0
A_T(1, 3) = 0
A_T(2, 3) = 0
A_T(3, 3) = alpha_grv*invx*(-U - l_i + 3)
A_T(4, 3) = alpha_grv*invx
A_T(5, 3) = 0
A_T(6, 3) = 0
A_T(1, 4) = invx*(As*U*alpha_grv - F_trb*U*V_g*alpha_grv*(V_g - 2)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1))
A_T(2, 4) = invx*(-F_trb*U*V_g*alpha_grv*(-V_g*alpha_gam + eta)/(F_trb*( &
      -V_g*alpha_gam + eta) + 1) + U*V_g*alpha_grv)
A_T(3, 4) = invx*(-F_trb*U*V_g*alpha_grv**2*eta/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) + alpha_grv*lambda)
A_T(4, 4) = alpha_grv*invx*(-U - l_i + 2)
A_T(5, 4) = invx*(-F_trb*U*V_g*alpha_grv*ups_T/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) - U*alpha_grv*ups_T)
A_T(6, 4) = 0
A_T(1, 5) = invx*(-F_trb*V*(V_g - 2)*(-V*c_kap_ad*nabla - c_dif + eta*( &
      -nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1)) + &
      V*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - om2) - 4* &
      nabla_ad)/f_rht)
A_T(2, 5) = invx*(-F_trb*V*(-V_g*alpha_gam + eta)*(-V*c_kap_ad*nabla - &
      c_dif + eta*(-nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + &
      eta) + 1)) + V*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + &
      nabla_ad))/f_rht)
A_T(3, 5) = invx*(-F_trb*V*alpha_grv*eta*(-V*c_kap_ad*nabla - c_dif + &
      eta*(-nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1 &
      )) + V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht)
A_T(4, 5) = V*alpha_grv*invx*nabla_ad/f_rht
A_T(5, 5) = invx*(-F_trb*V*ups_T*(-V*c_kap_ad*nabla - c_dif + eta*( &
      -nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1)) + &
      V*nabla*(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A_T(6, 5) = -V*invx*nabla/(c_rad*f_rht)
A_T(1, 6) = invx*(F_trb*(V_g - 2)*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) - V*alpha_egv*c_egv*nabla_ad - &
      V*c_eps_ad + alpha_hfl*c_rad*lambda*(-1 + nabla_ad/nabla))
A_T(2, 6) = invx*(F_trb*(-V_g*alpha_gam + eta)*(-V*alpha_egv*c_egv* &
      nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - &
      eta*f_conv)/(F_trb*(-V_g*alpha_gam + eta) + 1) + V*alpha_egv* &
      c_egv*nabla_ad + V*c_eps_ad - alpha_hfl*c_rad*lambda*nabla_ad/ &
      nabla + eta*f_conv)
A_T(3, 6) = invx*(F_trb*alpha_grv*eta*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) + alpha_grv*eta*f_conv)
A_T(4, 6) = 0
A_T(5, 6) = invx*(F_trb*ups_T*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad &
      + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/(F_trb*( &
      -V_g*alpha_gam + eta) + 1) + alpha_egv*c_egv + alpha_thm*c_thk* &
      i_omega_c + c_eps_S - alpha_hfl*c_rad*lambda/(V*nabla))
A_T(6, 6) = invx*(-l_i - 1)
