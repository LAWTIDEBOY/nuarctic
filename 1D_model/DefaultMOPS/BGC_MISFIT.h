C$Header: /Users/ikriest/CVS/mops/BGC_MISFIT.h,v 1.2 2018/03/12 06:44:38 ikriest Exp $
C$Name: mops-2_3 $

! arrays for diagnostics
      real*8 m1_out(bgc_ktotal),
     &       m2_out(bgc_ktotal),
     &       m3_out(bgc_ktotal)

      COMMON/COSTVARS/m1_out,m2_out,m3_out

