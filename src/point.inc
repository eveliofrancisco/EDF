      real(kind=8) ::     grad(3),hess(3,3),ehess(3),fielde(3),
     &                    fieldn(3),fplus(3),fminus(3),fieldt(3),
     &                    rho3(10),dsig
      real(kind=8) ::     rho,gradmod,lap,kink,king,kinp,vne,vnn,vee
      real(kind=8) ::     elf,rhospin,rhothis
      logical             domaxwell,dobfunc
      common /datum/      rho,grad,hess,gradmod,lap,ehess,
     &                    king,kink,vne,vnn,vee,fielde,fieldn, 
     &                    dsig(3),elf,kinp,fplus,fminus,rho3, 
     &                    rhospin,rhothis
      integer(kind=4) ::  nindex,ntype
      common /datui/      nindex,ntype
      common /datul/      domaxwell,dobfunc
      real(kind=8) ::     fxc,dernuc,derele,derho2,divfxc,hessel
      common/ehrdat/      fxc(3),dernuc(3),derele(3),derho2(3),
     &                    divfxc(3,3),hessel(3,3)
c
c$omp threadprivate(/datum/,/datui/,/ehrdat/)
