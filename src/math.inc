c math parameters

      real pi,rad
      parameter(pi=3.1415926535898,rad = pi/180.0)
      real tsin(-65000:65000),tcos(-65000:65000)
      real tacos(-10000:10000),texp(0:10000)
      real dlt(6, 6)
      common /math/tsin,tcos,tacos,texp,dlt
