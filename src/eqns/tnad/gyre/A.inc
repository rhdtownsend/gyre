A(1, 1) = del_x*(F_trb*(V_g - 2)*(V_g*alpha_gam - eta)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + V_g - l_i - 1)
A(2, 1) = del_x*(-As_iso + F_trb*(V_g - 2)**2/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) + om2)
A(3, 1) = 0
A(4, 1) = del_x*(As*U*alpha_grv - F_trb*U*V_g*alpha_grv*(V_g - 2)/(F_trb &
      *(-V_g*alpha_gam + eta) + 1))
A(5, 1) = del_x*(-F_trb*V*(V_g - 2)*(-V*c_kap_ad*nabla - c_dif + eta*( &
      -nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1)) + &
      V*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - om2) - 4* &
      nabla_ad)/f_rht)
A(6, 1) = del_x*(F_trb*(V_g - 2)*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) - V*alpha_egv*c_egv*nabla_ad - &
      V*c_eps_ad + alpha_hfl*c_rad*lambda*(-1 + nabla_ad/nabla))
A(1, 2) = del_x*(F_trb*(-V_g*alpha_gam + eta)*(V_g*alpha_gam - eta)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) - V_g*alpha_gam + eta)
A(2, 2) = del_x*(As + F_trb*(V_g - 2)*(-V_g*alpha_gam + eta)/(F_trb*( &
      -V_g*alpha_gam + eta) + 1) - U - l_i + 3)
A(3, 2) = 0
A(4, 2) = del_x*(-F_trb*U*V_g*alpha_grv*(-V_g*alpha_gam + eta)/(F_trb*( &
      -V_g*alpha_gam + eta) + 1) + U*V_g*alpha_grv)
A(5, 2) = del_x*(-F_trb*V*(-V_g*alpha_gam + eta)*(-V*c_kap_ad*nabla - &
      c_dif + eta*(-nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + &
      eta) + 1)) + V*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + &
      nabla_ad))/f_rht)
A(6, 2) = del_x*(F_trb*(-V_g*alpha_gam + eta)*(-V*alpha_egv*c_egv* &
      nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - &
      eta*f_conv)/(F_trb*(-V_g*alpha_gam + eta) + 1) + V*alpha_egv* &
      c_egv*nabla_ad + V*c_eps_ad - alpha_hfl*c_rad*lambda*nabla_ad/ &
      nabla + eta*f_conv)
A(1, 3) = del_x*(F_trb*alpha_grv*eta*(V_g*alpha_gam - eta)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + alpha_grv*eta)
A(2, 3) = F_trb*alpha_grv*del_x*eta*(V_g - 2)/(F_trb*(-V_g*alpha_gam + &
      eta) + 1)
A(3, 3) = alpha_grv*del_x*(-U - l_i + 3)
A(4, 3) = del_x*(-F_trb*U*V_g*alpha_grv**2*eta/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) + alpha_grv*lambda)
A(5, 3) = del_x*(-F_trb*V*alpha_grv*eta*(-V*c_kap_ad*nabla - c_dif + eta &
      *(-nabla + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1)) &
      + V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht)
A(6, 3) = del_x*(F_trb*alpha_grv*eta*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/( &
      F_trb*(-V_g*alpha_gam + eta) + 1) + alpha_grv*eta*f_conv)
A(1, 4) = 0
A(2, 4) = -alpha_grv*del_x
A(3, 4) = alpha_grv*del_x
A(4, 4) = alpha_grv*del_x*(-U - l_i + 2)
A(5, 4) = V*alpha_grv*del_x*nabla_ad/f_rht
A(6, 4) = 0
A(1, 5) = del_x*(F_trb*ups_T*(V_g*alpha_gam - eta)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + ups_T)
A(2, 5) = del_x*(F_trb*ups_T*(V_g - 2)/(F_trb*(-V_g*alpha_gam + eta) + 1 &
      ) + ups_T)
A(3, 5) = 0
A(4, 5) = del_x*(-F_trb*U*V_g*alpha_grv*ups_T/(F_trb*(-V_g*alpha_gam + &
      eta) + 1) - U*alpha_grv*ups_T)
A(5, 5) = del_x*(-F_trb*V*ups_T*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla &
      + nabla_ad))/(f_rht*(F_trb*(-V_g*alpha_gam + eta) + 1)) + V*nabla &
      *(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A(6, 5) = del_x*(F_trb*ups_T*(-V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + &
      alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv)/(F_trb*(-V_g* &
      alpha_gam + eta) + 1) + alpha_egv*c_egv + alpha_thm*c_thk* &
      i_omega_c + c_eps_S - alpha_hfl*c_rad*lambda/(V*nabla))
A(1, 6) = 0
A(2, 6) = 0
A(3, 6) = 0
A(4, 6) = 0
A(5, 6) = -V*del_x*nabla/(c_rad*f_rht)
A(6, 6) = del_x*(-l_i - 1)
