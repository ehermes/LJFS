module analyze2
    
    implicit none

    contains

    subroutine e_lj(ep1,ep2,sig1,sig2,rvec,elj)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: elj
        real*8 :: e_lj12, e_lj6

        call elj12(ep1,ep2,sig1,sig2,rvec,e_lj12)
        call elj6(ep1,ep2,sig1,sig2,rvec,e_lj6)

        elj = e_lj12 + e_lj6

        return
        end subroutine e_lj

    subroutine f_lj(ep1,ep2,sig1,sig2,rvec,flj)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: flj(3)
        real*8 :: f_lj12(3), f_lj6(3)

        call flj12(ep1,ep2,sig1,sig2,rvec,f_lj12)
        call flj6(ep1,ep2,sig1,sig2,rvec,f_lj6)
        
        flj = f_lj12 + f_lj6

        return
        end subroutine f_lj

    subroutine dedf_lj(ep1,ep2,sig1,sig2,rvec,dedep,dedsig,dfdep,dfdsig)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dedep, dedsig, dfdep(3), dfdsig(3)
        real*8 :: dedep_lj12, dedep_lj6, dedsig_lj12, dedsig_lj6
        real*8 :: dfdep_lj12(3), dfdep_lj6(3), dfdsig_lj12(3), dfdsig_lj6(3)

        call ddepelj12(ep1,ep2,sig1,sig2,rvec,dedep_lj12)
        call ddepelj6(ep1,ep2,sig1,sig2,rvec,dedep_lj6)

        dedep = dedep_lj12 + dedep_lj6

        call ddsigelj12(ep1,ep2,sig1,sig2,rvec,dedsig_lj12)
        call ddsigelj6(ep1,ep2,sig1,sig2,rvec,dedsig_lj6)

        dedsig = dedsig_lj12 + dedsig_lj6

        call ddepflj12(ep1,ep2,sig1,sig2,rvec,dfdep_lj12)
        call ddepflj6(ep1,ep2,sig1,sig2,rvec,dfdep_lj6)

        dfdep = dfdep_lj12 + dfdep_lj6

        call ddsigflj12(ep1,ep2,sig1,sig2,rvec,dfdsig_lj12)
        call ddsigflj6(ep1,ep2,sig1,sig2,rvec,dfdsig_lj6)

        dfdsig = dfdsig_lj12 + dfdsig_lj6

        return
        end subroutine dedf_lj

    subroutine dedf_coul(conv,q1,q2,rvec,dedq_coul,dfdq_coul)

        implicit none

        real*8, intent(in) :: conv, q1, q2, rvec(3)
        real*8, intent(out) :: dedq_coul, dfdq_coul(3)

        call ddqecoul(conv,q1,q2,rvec,dedq_coul)
        call ddqfcoul(conv,q1,q2,rvec,dfdq_coul)

        return
        end subroutine dedf_coul

    subroutine elj12(ep1,ep2,sig1,sig2,rvec,e_lj12)

        implicit none
        
        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: e_lj12
        real*8 :: ep, sigma, rsq 

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        e_lj12 = 4.0*ep*(sigma**12)*(rsq**-6)
        return
        end subroutine elj12

    subroutine elj6(ep1,ep2,sig1,sig2,rvec,e_lj6)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: e_lj6
        real*8 :: ep, sigma, rsq

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        e_lj6 = -4.0*ep*(sigma**6)*(rsq**-3)
        return
        end subroutine elj6

    subroutine ecoul(conv,q1,q2,rvec,e_coul)
        
        implicit none

        real*8, intent(in) :: conv, q1, q2, rvec(3)
        real*8, intent(out) :: e_coul
        real*8 :: rsq, rmag

        rsq = dot_product(rvec,rvec)
        rmag = sqrt(rsq)
        e_coul = conv * q1 * q2 / rmag
        return
        end subroutine ecoul

    subroutine flj12(ep1,ep2,sig1,sig2,rvec,f_lj12)
        
        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: f_lj12(3)
        real*8 :: ep, sigma, rsq
        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        f_lj12 = (48.0*ep*rvec)*(sigma**12)*(rsq**-7)
        return
        end subroutine flj12

    subroutine flj6(ep1,ep2,sig1,sig2,rvec,f_lj6)
        
        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: f_lj6(3)
        real*8 :: ep, sigma, rsq
        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        f_lj6 = (-24.0*ep*rvec)*(sigma**6)*(rsq**-4)
        return
        end subroutine flj6

    subroutine fcoul(conv,q1,q2,rvec,f_coul)

        implicit none

        real*8, intent(in) :: conv, q1, q2, rvec(3)
        real*8, intent(out) :: f_coul(3)
        real*8 :: rsq, rmag

        rsq = dot_product(rvec,rvec)
        rmag = sqrt(rsq)
        f_coul = conv * q1 * q2 * (rvec/rmag) / rsq
        return
        end subroutine fcoul

    subroutine ddepelj12(ep1,ep2,sig1,sig2,rvec,dedep_lj12)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dedep_lj12
        real*8 :: dep, sigma, rsq

        dep = 0.5 * sqrt(ep2/ep1)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        dedep_lj12 = 4.0 * dep * (sigma**12) * (rsq**-6)
        return
        end subroutine ddepelj12

    subroutine ddepelj6(ep1,ep2,sig1,sig2,rvec,dedep_lj6)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dedep_lj6
        real*8 :: dep, sigma, rsq

        dep = 0.5 * sqrt(ep2/ep1)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        dedep_lj6 = -4.0 * dep * (sigma**6) * (rsq**-3)
        return
        end subroutine ddepelj6

    subroutine ddsigelj12(ep1,ep2,sig1,sig2,rvec,dedsig_lj12)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dedsig_lj12
        real*8 :: ep, sigma, dsigma, rsq

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        dsigma = 0.5 * sqrt(sig2/sig1)
        rsq = dot_product(rvec,rvec)
        dedsig_lj12 = 48.0*ep*(sigma**11)*(rsq**-6)*dsigma
        return
        end subroutine ddsigelj12

    subroutine ddsigelj6(ep1,ep2,sig1,sig2,rvec,dedsig_lj6)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dedsig_lj6
        real*8 :: ep, sigma, dsigma, rsq

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        dsigma = 0.5 * sqrt(sig2/sig1)
        rsq = dot_product(rvec,rvec)
        dedsig_lj6 = -24.0*ep*(sigma**5)*(rsq**-3)*dsigma
        return
        end subroutine ddsigelj6

    subroutine ddqecoul(conv,q1,q2,rvec,dedq_coul)

        implicit none

        real*8, intent(in) :: conv, q1, q2, rvec(3)
        real*8, intent(out) :: dedq_coul
        real*8 :: rsq, rmag
        rsq = dot_product(rvec,rvec)
        rmag = sqrt(rsq)
        dedq_coul = conv * q2 / rmag
        return
        end subroutine ddqecoul

    subroutine ddepflj12(ep1,ep2,sig1,sig2,rvec,dfdep_lj12)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dfdep_lj12(3)
        real*8 :: dep, sigma, rsq

        dep = 0.5 * sqrt(ep2/ep1)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        dfdep_lj12 = (48.0 * dep * rvec) * (sigma**12) * (rsq**-7)
        return
        end subroutine ddepflj12

    subroutine ddepflj6(ep1,ep2,sig1,sig2,rvec,dfdep_lj6)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dfdep_lj6(3)
        real*8 :: dep, sigma, rsq

        dep = 0.5 * sqrt(ep2/ep1)
        sigma = sqrt(sig1*sig2)
        rsq = dot_product(rvec,rvec)
        dfdep_lj6 = (-24.0 * dep * rvec) * (sigma**6) * (rsq**-4)
        return
        end subroutine ddepflj6

    subroutine ddsigflj12(ep1,ep2,sig1,sig2,rvec,dfdsig_lj12)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dfdsig_lj12(3)
        real*8 :: ep, sigma, dsigma, rsq

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        dsigma = 0.5 * sqrt(sig2/sig1)
        rsq = dot_product(rvec,rvec)
        dfdsig_lj12 = (576.0 * ep * rvec) * (sigma**11) * (rsq**-7) * dsigma
        return
        end subroutine ddsigflj12

    subroutine ddsigflj6(ep1,ep2,sig1,sig2,rvec,dfdsig_lj6)

        implicit none

        real*8, intent(in) :: ep1, ep2, sig1, sig2, rvec(3)
        real*8, intent(out) :: dfdsig_lj6(3)
        real*8 :: ep, sigma, dsigma, rsq

        ep = sqrt(ep1*ep2)
        sigma = sqrt(sig1*sig2)
        dsigma = 0.5 * sqrt(sig2/sig1)
        rsq = dot_product(rvec,rvec)
        dfdsig_lj6 = (-144.0 * ep * rvec) * (sigma**5) * (rsq**-4) * dsigma
        return
        end subroutine ddsigflj6

    subroutine ddqfcoul(conv,q1,q2,rvec,dfdq_coul)

        implicit none

        real*8, intent(in) :: conv, q1, q2, rvec(3)
        real*8, intent(out) :: dfdq_coul(3)
        real*8 :: rsq, rmag

        rsq = dot_product(rvec,rvec)
        rmag = sqrt(rsq)
        dfdq_coul = conv * q2 * (rvec/rmag) / rsq
        return
        end subroutine ddqfcoul

end module analyze2
