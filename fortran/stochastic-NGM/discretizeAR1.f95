module discretizeAR1
use constants
implicit none
contains


     !Discretizes an AR(1) of the following form:
     !z' = rho*z + eps
     !eps ~ N(mu,sigma^2)
     !with an evenly-spaced grid over
     ![mu_z-m*sigma_z,mu_z+m*sigma_z] and
     !a PTM that is based on the areas
     !under the conditional densities.
     !
     !zvect (1 x N) is a vector of the grid points.
     !Pmat (N x N) is the transition probability matrix.
     !Pmat(i,j) is the probability of z'=zvect(j) given z=zvect(i).
     !
      subroutine Tauchen(rho, mu_eps, sigma_eps, N, m, zvect, Pmat)
         use mathutils
          use distributions
        implicit none

          real(dbl), intent(in):: rho, mu_eps, sigma_eps, m
          integer, intent(in):: N
          real(dbl), intent(out):: zvect(N)
          real(dbl), intent(out):: Pmat(N,N)

         real(dbl) mu_z, sigma_z, sigma_zcond
         real(dbl) w
         integer i,j

         mu_z = mu_eps/(1-rho)
         sigma_z = sigma_eps/sqrt(1-rho**2)
          sigma_zcond = sigma_eps

         zvect = linspace(mu_z-m*sigma_z,mu_z+m*sigma_z,N)
        w = (zvect(2)-zvect(1))/2


          do i=1,N
              Pmat(i,1) = enordf((zvect(1) + w - rho*zvect(i) - mu_eps)/sigma_eps)
             do j=2,N-1
               Pmat(i,j) = enordf((zvect(j) + w- rho*zvect(i) - mu_eps)/sigma_eps) - &
                  enordf((zvect(j) - w- rho*zvect(i) - mu_eps)/sigma_eps)
             end do
              Pmat(i,N) = 1 - enordf((zvect(N) - w - rho*zvect(i) - mu_eps)/sigma_eps)
         end do
     end subroutine Tauchen


     subroutine Rouwenhurst(rho, mu_eps, sigma_eps, N, zvect, Pmat)
         use mathutils
         implicit none

         real(dbl), intent(in):: rho, mu_eps, sigma_eps
         integer, intent(in):: N
         real(dbl), intent(out):: zvect(N)
         real(dbl), intent(out):: Pmat(N,N)

         real(dbl) mu_z, sigma_z, q, eps
         real(dbl), allocatable, dimension(:,:):: P1, P2
         integer status, i, j


         mu_z = mu_eps/(1-rho)
         sigma_z = sigma_eps/sqrt(1-rho**2)

         q = (rho+1)/2
         eps = sqrt(dble(N-1)) * sigma_z

         if (N == 1) then
             Pmat = 1.0d0
             zvect = mu_z
             return
         else if (N == 2) then
             Pmat = reshape((/q, 1-q, 1-q, q/),(/2,2/))
             zvect = (/mu_z-eps,mu_z+eps/)
             return
         end if

         allocate(P1(2,2),stat=status)
         P1 = reshape((/q, 1-q, 1-q, q/),(/2,2/))

         do i=2,N-1
             allocate(P2(i+1,i+1),stat=status)
            P2 = q * reshape( (/  (/(P1(:,j),0.0d0 ,j=1,i)/) ,  (/(0.0d0,j=1,i+1)/)    /), (/i+1,i+1/) ) + &
                  (1-q) * reshape( (/  (/(0.0d0,j=1,i+1)/), (/ (P1(:,j),0.0d0 ,j=1,i)/)   /) ,   (/i+1,i+1/) ) + &
                  (1-q) * reshape( (/  (/ (0.0d0,P1(:,j) ,j=1,i) /) ,  (/(0.0d0,j=1,i+1)/)  /), (/i+1,i+1/) ) + &
                  q * reshape( (/ (/(0.0d0,j=1,i+1)/), (/(0.0d0,P1(:,j) ,j=1,i)/)   /) ,   (/i+1,i+1/) )

            P2(2:i,:) = P2(2:i,:)/2

             deallocate(P1,stat=status)

             if (i==N-1) then
                 Pmat = P2
             else
                 allocate(P1(i+1,i+1), stat=status)
                 P1 = P2
             end if

             deallocate(P2,stat=status)
         end do

         zvect = linspace(mu_z-eps,mu_z+eps,N)

    end subroutine Rouwenhurst


end module discretizeAR1
