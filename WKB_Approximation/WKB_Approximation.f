!     Numerical analysis of the ordinary differential equations
!     with Runge-Kutta method.
!===============================================================================
!     微分方程式 y"+e^x*y=0 を数値的に解くプログラム。
!     初期条件はy(x=0)=1,y'(x=0)=0
!     出力は(x,y)を出力する。
!     この出力結果はWKB近似で求めた解と比べる。それによりWKB近似について理解深め、
!     特に有効条件についての感覚的理解を得ることが目的である。
!===============================================================================
!-----文字および定数の定義---
      program WKB_Approximation
      implicit none
      character FNAME(4)*20,ERROR*20
      real*8 s(4,50000),x(50000),r(4,50000)
      real*8 x0,x1,d,x2
      integer n,ni,l,i,j,nfile,ix2,ix3
      data n,x0,x1,d,ni/2, 0.0d0, 1.0d1, 1.0d-3, 10/
      data nfile,x2/1, -9.998d0/
!================================================================================
!     積分範囲は[x2,x1]
!     区間[x2,x0]での積分をまず行い、次に
!     区間[x0,x1]での積分を行う。
!     初期値を与える点はx0 (x2<x0<x1)
!     nは方程式の数
!     dはサンプリングの間隔
!     nfileは出力したい(x,y(i))の数(つくるファイルの数)
!     ix2はx(ix2)=x2となるところ。区間[x2,x0]の積分で
!     最後に値を入れるところ。
!     niは積分の計算のときの分割数をd/niとする。
!
!     ni=10のときと比較して、ni=2とすると結果は最大10^-6程度変わる。(特に激しく振動するところで)
!     ni=4のときは差は最大10^-8程度であるが
!     5<=ni<10ときは少なくとも差は10^-8以下である。ni=20でも同じ。したがって、
!     この計算(ni=10)は少なくとも小数点7桁まで信用できると思われる。
!===============================================================================
!-----初期値および条件の書き出し---
      write(*,*) ' *** Numerical Analysis of ODE ***'
      write(*,*) '      (With Runge-Kutta Method)'
      write(*,*)
      write(*,200) n,x2,x1,d,ni  
      l=1
      x(1)=x0
      call INITIAL(s)
      call INITIAL(r)
      write(*,*) 'Initial condition'
      write(*,210) x(1),(s(i,1),i=1,n)
!-----RungeKutta法で微分方程式を計算-------------------------------------------
      call DRUNGE(l,x,n,s,d,ni,x1,ix2)
      call DRUNGE(l,x,n,r,d,7,x1,ix3)        !誤差見積もり用(niを変えたもの)
!-----微分方程式の解の画面への書き出し-----------------------------------------
!      write(*,*) ' * solution *'
!       i=1
! 10   write(*,210) x(ix2+1-i),(s(j,i),j=1,n)
!      if (x(i).ge.x2) then
!           i=i+1
!            goto 10
!      endif
!      i=i+2
!      
!      
! 30   write(*,210) x(i),(s(j,i),j=1,n)
!      if (x(i).le.x1) then
!            i=i+1
!           goto 30
!      endif     
!-----常微分方程式の解のファイルへの書き出し--------------------------
      write(6,*)'Input the output filename.'
      do j=1,nfile
       write(6,230) j,' th'
       read(5,*) FNAME(j)
      end do
      do j=1,nfile
       open(unit=j,file=FNAME(j))        !!!ファイルを開く
      end do
      
       write(6,*) 'Input the filename of error.'  !誤差見積もりようのファイルをつくる。
       read(5,*) ERROR
       open(unit=7,file=ERROR) 
       
!-----ファイルにデータを書き込む------------------------------------ 
      i=1
!===================================================================
!     まず区間[x2,x0]における解
!     x(i)にはx(1)=x0から始めてx2までiとは"逆向きに"値を入れて
!     いる。なのでx(i)の値の小さい順に出力されるようにしている。
!     つまりx(ix2)からx(1)までiの降順に書き込んでいる。
!     ただし判定条件はiの昇順に見ていることに注意。
!===================================================================
 20   do j=1,nfile
       write(j,220) x(ix2+1-i),s(j,ix2+1-i)
      end do
      
       write(7,220) x(ix2+1-i),s(1,ix2+1-i)-r(1,ix2+1-i)  !誤差評価用
       
      if (x(i).ge.x2) then
            i=i+1
            goto 20
      endif
      write(6,*) ix2
      write(6,*) ix3
      i=i+2    ! x(ix2+1)はx0なのでこれはx(1)と被る。新たに必要なのはx(ix2+2)以降
      
!-----区間(x0,x1]における解-----------------------------------------------------

 40   do j=1,nfile
      write(j,220) x(i),s(j,i)
      end do
      
      write(7,220) x(i),s(1,i)-r(1,i)     !誤差評価用
      
      if (x(i).le.x1) then
            i=i+1
            goto 40
      endif     
      
      do j=1,nfile       !ファイルを閉じる。
      close(j)
      end do
      
      close(7)    !誤差評価用のファイルを閉じる。
      
!-----書式の設定-----------------------------------------------------
 200   format(1h ,'number of equations',i2/' upper limit=',f5.1
     1/' lower limit=',f5.1/' sampling interval=',f6.3
     2/' division number=',i2)
 210   format(1h ,2x,'x=',f4.1,4x,'y=(',f8.5,','
     &,f8.5,')')
 220   format(f8.5,' ',f10.7)            !ファイル出力のフォーマット
 230   format(' ',i1,a)
      end
!------------------------------------------------
!
!-----初期値代入---------------------------------
      subroutine INITIAL(s)
      real*8 s(4,50000)
      s(1,1)=1.0d0
      s(2,1)=0.0d0
      s(3,1)=1.0d0
      s(4,1)=1.0d0
      return
      end
!------------------------------------------------
!
!-----微分方程式---------------------------------
      subroutine EQATION(x,y,f)
      real*8 x,y(4),f(4)
      f(1)=y(2)
      f(2)=-dexp(x)*y(1)
      f(3)=y(2)
      f(4)=y(3)
      return
      end
!------------------------------------------------
!====================================================================
!-----Runge-Kutta-Methodのドライバールーチン-----
!     x(i)は独立変数x。s(j,i)は求めたい関数y_j
!====================================================================
      subroutine DRUNGE(l,x,n,s,d,ni,x1,ix2)  !出力はx,s,l
      real*8 s(4,50000),y(4)
      real*8 x(50000),h,d,x1
      integer ni,l,i,ix2
      h=d/real(ni)
!-----まず区間[x2,x0]での積分-----------------------------------------
      y(1)=s(1,1)
      y(2)=s(2,1)
      y(3)=s(3,1)
      y(4)=s(4,1)      !y(i)をもとにy(i+1)を計算するため
      l=1
 10   x(l+1)=x(l)
        do i=1,ni
           call RUNGE(-h,n,x(l+1),y)    !x(l+1)におけるy_j,つまりs(j,l+1)になる量を計算。
        end do
        l=l+1
        s(1,l)=y(1)
        s(2,l)=y(2)
        s(3,l)=y(3)
        s(4,l)=y(4)
        if(x(l).le.-9.998d0) goto 40
      goto 10
 40   ix2=l      !前半の区間[x2,x0]の積分でどこまでx(i)を使ったかをメモ。
!-----次に区間[x0,x1]での積分--------------------------------------------
      l=l+1
      x(l)=0.0d0
      y(1)=s(1,1)
      y(2)=s(2,1)
      y(3)=s(3,1)
      y(4)=s(4,1)        !再び初期値を入れてx0から積分してゆく。
 20   x(l+1)=x(l)
        do i=1,ni
           call RUNGE(h,n,x(l+1),y)    !x(l+1)におけるy_j,つまりs(j,l+1)になる量を計算。
        end do
        l=l+1
        s(1,l)=y(1)
        s(2,l)=y(2)
        s(3,l)=y(3)
        s(4,l)=y(4)
        if(x(l).ge.x1) goto 30
      goto 20

 30   return
      end
!------------------------------------------------------------------------
!
!-----Runge-Kutta-Method-------------------------------------------------
      subroutine RUNGE(h,n,x,y)   ! 出力はx,y
      real*8 f(4),y(4),s1(4),s2(4),s3(4),s4(4),ysav(4)
      real*8 x,xsav,h
      xsav=x
      do 10 j=1,n
            ysav(j)=y(j)
 10   continue
      call EQATION(x,y,f)
      do 20 j=1,n
            s1(j)=h*f(j)
 20   continue
 
      x=xsav+0.5d0*h                 ! x(i+1)=x(i)+0.5h
      
      do 30 j=1,n
            y(j)=ysav(j)+0.5d0*s1(j)
 30   continue
      call EQATION(x,y,f)
      do 40 j=1,n
            s2(j)=h*f(j)
            y(j)=ysav(j)
 40   continue
      x=xsav+0.5d0*h
      do 50 j=1,n
            y(j)=ysav(j)+0.5d0*s2(j)
 50   continue
      call EQATION(x,y,f)
      do 60 j=1,n
            s3(j)=h*f(j)
            y(j)=ysav(j)
 60   continue
      x=xsav+h
      do 70 j=1,n
            y(j)=ysav(j)+s3(j)
 70   continue
      call EQATION(x,y,f)
      do 80 j=1,n
            s4(j)=h*f(j)
            y(j)=ysav(j)
 80   continue
      do 90 j=1,n
            y(j)=ysav(j)+(s1(j)+2.0d0*s2(j)+2.0d0*s3(j)+s4(j))/6.0d0
 90   continue
      return
      end
!------------------------------------------------
