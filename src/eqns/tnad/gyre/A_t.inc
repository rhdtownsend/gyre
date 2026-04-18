A_t(1, 1) = del_x*(F_trb*(V_g - 2)*(V_g*alpha_gam - eta) + F_trb*(-V_g* &
      alpha_gam + eta)*(V_g*alpha_gam - eta)*(F_trb*V_g*om2 - 2*F_trb* &
      om2)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + V_g - l_i - &
      1)
A_t(2, 1) = del_x*(-F_trb*om2*(-V_g*alpha_gam + eta)*(V_g*alpha_gam - &
      eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) - V_g* &
      alpha_gam + eta)
A_t(3, 1) = del_x*(F_trb**2*alpha_grv*lambda*(-V_g*alpha_gam + eta)*(V_g &
      *alpha_gam - eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) &
      + F_trb*alpha_grv*eta*(V_g*alpha_gam - eta) + alpha_grv*eta)
A_t(4, 1) = 0
A_t(5, 1) = del_x*(F_trb**2*om2*ups_T*(-V_g*alpha_gam + eta)*(V_g* &
      alpha_gam - eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + &
      F_trb*ups_T*(V_g*alpha_gam - eta) + ups_T)
A_t(6, 1) = 0
A_t(1, 2) = del_x*(-As_iso + F_trb*(V_g - 2)**2 + F_trb*(V_g - 2)*(-V_g* &
      alpha_gam + eta)*(F_trb*V_g*om2 - 2*F_trb*om2)/(F_trb*V_g* &
      alpha_gam*om2 - F_trb*lambda - om2) + om2)
A_t(2, 2) = del_x*(As - F_trb*om2*(V_g - 2)*(-V_g*alpha_gam + eta)/( &
      F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) - U - l_i + 3)
A_t(3, 2) = del_x*(F_trb**2*alpha_grv*lambda*(V_g - 2)*(-V_g*alpha_gam + &
      eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + F_trb* &
      alpha_grv*eta*(V_g - 2))
A_t(4, 2) = -alpha_grv*del_x
A_t(5, 2) = del_x*(F_trb**2*om2*ups_T*(V_g - 2)*(-V_g*alpha_gam + eta)/( &
      F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + F_trb*ups_T*(V_g &
      - 2) + ups_T)
A_t(6, 2) = 0
A_t(1, 3) = 0
A_t(2, 3) = 0
A_t(3, 3) = alpha_grv*del_x*(-U - l_i + 3)
A_t(4, 3) = alpha_grv*del_x
A_t(5, 3) = 0
A_t(6, 3) = 0
A_t(1, 4) = del_x*(As*U*alpha_grv - F_trb*U*V_g*alpha_grv*(V_g - 2) - &
      F_trb*U*V_g*alpha_grv*(-V_g*alpha_gam + eta)*(F_trb*V_g*om2 - 2* &
      F_trb*om2)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2))
A_t(2, 4) = del_x*(F_trb*U*V_g*alpha_grv*om2*(-V_g*alpha_gam + eta)/( &
      F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + U*V_g*alpha_grv)
A_t(3, 4) = del_x*(-F_trb**2*U*V_g*alpha_grv**2*lambda*(-V_g*alpha_gam + &
      eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) - F_trb*U*V_g &
      *alpha_grv**2*eta + alpha_grv*lambda)
A_t(4, 4) = alpha_grv*del_x*(-U - l_i + 2)
A_t(5, 4) = del_x*(-F_trb**2*U*V_g*alpha_grv*om2*ups_T*(-V_g*alpha_gam + &
      eta)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) - F_trb*U*V_g &
      *alpha_grv*ups_T - U*alpha_grv*ups_T)
A_t(6, 4) = 0
A_t(1, 5) = del_x*(-F_trb*V*(V_g - 2)*(-V*c_kap_ad*nabla - c_dif + eta*( &
      -nabla + nabla_ad))/f_rht - F_trb*V*(-V_g*alpha_gam + eta)*(F_trb &
      *V_g*om2 - 2*F_trb*om2)*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla &
      + nabla_ad))/(f_rht*(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2 &
      )) + V*(V*c_kap_ad*nabla + c_dif + 4*nabla + nabla_ad*(U - om2) - &
      4*nabla_ad)/f_rht)
A_t(2, 5) = del_x*(F_trb*V*om2*(-V_g*alpha_gam + eta)*(-V*c_kap_ad*nabla &
      - c_dif + eta*(-nabla + nabla_ad))/(f_rht*(F_trb*V_g*alpha_gam* &
      om2 - F_trb*lambda - om2)) + V*(-V*c_kap_ad*nabla - c_dif + eta*( &
      -nabla + nabla_ad))/f_rht)
A_t(3, 5) = del_x*(-F_trb**2*V*alpha_grv*lambda*(-V_g*alpha_gam + eta)*( &
      -V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/(f_rht*( &
      F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2)) - F_trb*V* &
      alpha_grv*eta*(-V*c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad &
      ))/f_rht + V*alpha_grv*eta*(-nabla + nabla_ad)/f_rht)
A_t(4, 5) = V*alpha_grv*del_x*nabla_ad/f_rht
A_t(5, 5) = del_x*(-F_trb**2*V*om2*ups_T*(-V_g*alpha_gam + eta)*(-V* &
      c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/(f_rht*(F_trb* &
      V_g*alpha_gam*om2 - F_trb*lambda - om2)) - F_trb*V*ups_T*(-V* &
      c_kap_ad*nabla - c_dif + eta*(-nabla + nabla_ad))/f_rht + V*nabla &
      *(-c_kap_S + 4*f_rht)/f_rht - df_rht - l_i + 2)
A_t(6, 5) = -V*del_x*nabla/(c_rad*f_rht)
A_t(1, 6) = del_x*(F_trb*(V_g - 2)*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv) + &
      F_trb*(-V_g*alpha_gam + eta)*(F_trb*V_g*om2 - 2*F_trb*om2)*(-V* &
      alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda* &
      nabla_ad/nabla - eta*f_conv)/(F_trb*V_g*alpha_gam*om2 - F_trb* &
      lambda - om2) - V*alpha_egv*c_egv*nabla_ad - V*c_eps_ad + &
      alpha_hfl*c_rad*lambda*(-1 + nabla_ad/nabla))
A_t(2, 6) = del_x*(-F_trb*om2*(-V_g*alpha_gam + eta)*(-V*alpha_egv*c_egv &
      *nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - &
      eta*f_conv)/(F_trb*V_g*alpha_gam*om2 - F_trb*lambda - om2) + V* &
      alpha_egv*c_egv*nabla_ad + V*c_eps_ad - alpha_hfl*c_rad*lambda* &
      nabla_ad/nabla + eta*f_conv)
A_t(3, 6) = del_x*(F_trb**2*alpha_grv*lambda*(-V_g*alpha_gam + eta)*(-V* &
      alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda* &
      nabla_ad/nabla - eta*f_conv)/(F_trb*V_g*alpha_gam*om2 - F_trb* &
      lambda - om2) + F_trb*alpha_grv*eta*(-V*alpha_egv*c_egv*nabla_ad &
      - V*c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv &
      ) + alpha_grv*eta*f_conv)
A_t(4, 6) = 0
A_t(5, 6) = del_x*(F_trb**2*om2*ups_T*(-V_g*alpha_gam + eta)*(-V* &
      alpha_egv*c_egv*nabla_ad - V*c_eps_ad + alpha_hfl*c_rad*lambda* &
      nabla_ad/nabla - eta*f_conv)/(F_trb*V_g*alpha_gam*om2 - F_trb* &
      lambda - om2) + F_trb*ups_T*(-V*alpha_egv*c_egv*nabla_ad - V* &
      c_eps_ad + alpha_hfl*c_rad*lambda*nabla_ad/nabla - eta*f_conv) + &
      alpha_egv*c_egv + alpha_thm*c_thk*i_omega_c + c_eps_S - alpha_hfl &
      *c_rad*lambda/(V*nabla))
A_t(6, 6) = del_x*(-l_i - 1)
