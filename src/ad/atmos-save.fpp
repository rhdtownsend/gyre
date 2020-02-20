      ! Set up elements of the Jacobian matrix

      A(1,1) = V_g - 3._WP
      A(1,2) = lambda/(c_1*alpha_om*omega_c**2) - V_g
      A(1,3) = alpha_gr*(lambda/(c_1*alpha_om*omega_c**2))
      A(1,4) = alpha_gr*(0._WP)

      A(2,1) = c_1*alpha_om*omega_c**2 - As
      A(2,2) = As + 1._WP
      A(2,3) = alpha_gr*(0._WP)
      A(2,4) = alpha_gr*(-1._WP)

      A(3,1) = alpha_gr*(0._WP)
      A(3,2) = alpha_gr*(0._WP)
      A(3,3) = alpha_gr*(1._WP)
      A(3,4) = alpha_gr*(1._WP)

      A(4,1) = alpha_gr*(0._WP)
      A(4,2) = alpha_gr*(0._WP)
      A(4,3) = alpha_gr*(lambda)
      A(4,4) = alpha_gr*(0._WP)

      ! Evaluate the eigenvalues

      psi2 = (A(2,2) - A(1,1))**2 + 4._WP*A(1,2)*A(2,1)

      if (psi2 < 0._WP) then

         if (.NOT. warned) then
            $WARN(WARNING: Discarding imaginary part of atmospheric radial wavenumber)
            warned = .TRUE.
         endif

      endif

      beta(1) = 0.5_WP*(A(1,1) + A(2,2) + sqrt(psi2))
      beta(2) = 0.5_WP*(A(1,1) + A(2,2) - sqrt(psi2))
      beta(3) = -l_e
      beta(4) = l_e + 1._WP

      ! Evaluate the associated right eigenvectors

      U(1,1) = A(1,2)
      U(2,1) = beta(1) - A(1,1)
      U(3,1) = 0._WP
      U(4,1) = 0._WP

      U(1,2) = A(1,2)
      U(2,2) = beta(2) - A(1,1)
      U(3,2) = 0._WP
      U(4,2) = 0._WP

      U(1,3) = ((A(2,2)+l_e)*(-A(1,3) + A(1,4)*(l_e+1._WP)) - A(1,2)*(-A(2,3) + A(2,4)*(l_e+1))) * &
               ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))           
      U(2,3) = ((A(1,1)+l_e)*(-A(2,3) + A(2,4)*(l_e+1._WP)) - A(2,1)*(-A(1,3) + A(1,4)*(l_e+1))) * &
               ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      U(3,3) = 1._WP * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      U(4,3) = -(l_e-1._WP) * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))

      U(1,4) = (-(A(2,2)-l_e-1._WP)*(A(1,3) + A(1,4)*l_e) + A(1,2)*(A(2,3) + A(2,4)*l_e)) * &
               ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1))
      U(2,4) = (-(A(1,1)-l_e-1._WP)*(A(2,3) + A(2,4)*l_e) + A(2,1)*(A(1,3) + A(1,4)*l_e)) * &
               ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1))
      U(3,4) = 1._WP * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      U(4,4) = l_e * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))


      ! Evaluate the associated left eigenvectors

      W(1,1) = -(beta(2) - A(1,1)) * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      W(1,2) = A(1,2) * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      W(1,3) = ((beta(2) - A(1,1))*(U(1,3)*l_e + U(1,4)*(l_e+1._WP)) - A(1,2)*(U(2,3)*l_e + U(2,4)*(l_e+1._WP)))/ &
               (2._WP*l_e+1._WP)
      W(1,4) = (-(beta(2) - A(1,1))*(U(1,3) - U(1,4)) + A(1,2)*(U(2,3) - U(2,4)))/ &
               (2._WP*l_e+1._WP)

      ! W(1,1) = beta(2) - A(1,1)
      ! W(1,2) = -A(1,2)
      ! W(1,3) = (-(beta(2) - A(1,1))*(U(1,3)*l_e + U(1,4)*(l_e+1._WP)) + A(1,2)*(U(2,3)*l_e + U(2,4)*(l_e+1._WP)))/ &
      !          (2._WP*l_e+1._WP)
      ! W(1,4) = ((beta(2) - A(1,1))*(U(1,3) - U(1,4)) - A(1,2)*(U(2,3) - U(2,4)))/ &
      !          (2._WP*l_e+1._WP)

      W(2,1) = -(beta(1) - A(1,1)) * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      W(2,2) = A(1,2) * ((A(1,1)+l_e)*(A(2,2)+l_e) - A(1,2)*A(2,1)) * ((A(1,1)-l_e-1._WP)*(A(2,2)-l_e-1._WP) - A(1,2)*A(2,1))
      W(2,3) = ((beta(1) - A(1,1))*(U(1,3)*l_e + U(1,4)*(l_e+1._WP)) - A(1,2)*(U(2,3)*l_e + U(2,4)*(l_e+1._WP)))/ &
               (2._WP*l_e+1._WP)
      W(2,4) = (-(beta(1) - A(1,1))*(U(1,3) - U(1,4)) + A(1,2)*(U(2,3) - U(2,4)))/ &
               (2._WP*l_e+1._WP)

      W(3,1) = 0._WP
      W(3,2) = 0._WP
      W(3,3) = l_e/(2._WP*l_e+1._WP)
      W(3,4) = -1._WP/(2._WP*l_e+1._WP)

      W(4,1) = 0._WP
      W(4,2) = 0._WP
      W(4,3) = (l_e+1._WP)/(2._WP*l_e+1._WP)
      W(4,4) = 1._WP/(2._WP*l_e+1._WP)

      ! Perform validations on the eigenvectors

      I = identity_matrix(4)

      print *,'Chk right:', &
           NORM2(MATMUL(A-beta(1)*I, U(:,1)))/NORM2(U(:,1)), &
           NORM2(MATMUL(A-beta(2)*I, U(:,2)))/NORM2(U(:,2)), &
           NORM2(MATMUL(A-beta(3)*I, U(:,3)))/NORM2(U(:,3)), &
           NORM2(MATMUL(A-beta(4)*I, U(:,4)))/NORM2(U(:,4))

      print *,'Chk left:', &
           NORM2(MATMUL(W(1,:), A-beta(1)*I))/NORM2(W(1,:)), &
           NORM2(MATMUL(W(2,:), A-beta(2)*I))/NORM2(W(2,:)), &
           NORM2(MATMUL(W(3,:), A-beta(3)*I))/NORM2(W(3,:)), &
           NORM2(MATMUL(W(4,:), A-beta(4)*I))/NORM2(W(4,:))

