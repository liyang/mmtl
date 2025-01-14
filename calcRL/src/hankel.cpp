/************************* HANKEL FUNCTIONS **********************/

//#include "stdafx.h"
#include <math.h>
#include "complex.h"
#include "hankel.h"


#define abs(a) ((a)>=0 ? (a) : -(a))

double recurrence( int, int, double, double* );


/***** Function Hank20 - to calculate H^(2)_0(x), x - double  ****/

Complex Hank20( double x ) {

   double vj0, vy0;

   hank01( vj0, vy0, x, 1 );

   return( Complex( vj0, -vy0 ) );
}

/***** Function Hank21 - to calculate H^(2)_1(x), x - double  ****/

Complex Hank21( double x ) {

   double vj1, vy1;

   hank11( vj1, vy1, x, 1 );

   return( Complex( vj1, -vy1 ) );
}

/***** Function Hank10 - to calculate H^(1)_0(x), x - double  ****/

Complex Hank10( double x ) {

   double vj0, vy0;

   hank01( vj0, vy0, x, 1 );

   return( Complex( vj0, vy0 ) );
}

/***** Function Hank11 - to calculate H^(1)_1(x), x - double  ****/

Complex Hank11( double x ) {

   double vj1, vy1;

   hank11( vj1, vy1, x, 1 );

   return( Complex( vj1, vy1 ) );
}

int hank01( double& vj0, double& vy0, double xd, int n ) {

 int    n1, n2;
 double x, y, z, fx, x1, x2, x3, x4;
 double xlg = 1.0e+70;

 double a[]={0.0,
     -.17e-18                  , .1222e-16             ,
     -.75885e-15               , .4125321e-13          ,
     -.194383469e-11           , .7848696314e-10       ,
     -.267925353056e-08        , .7608163592419e-07    ,
     -.176194690776215e-05     , .3246032882100508e-04 ,
     -.46062616620627505e-03   , .48191800694676045e-02,
     -.34893769411408885e-01   , .15806710233209726    ,
     -.37009499387264978       , .26517861320333681    ,
     -.87234423528522213e-02   , .31545594294978024    ,
     -.1e-19                   , .39e-18               ,
     -.2698e-16                , .164349e-14           ,
     -.8747341e-13             , .402633082e-11        ,
     -.15837552542e-09         , .524879478733e-08     ,
     -.14407233274019e-06      , .32065325376548e-05   ,
     -.5632079141056987e-04    , .75311359325777423e-03,
     -.72879624795520792e-02   , .47196689595763387e-01,
     -.17730201278114358       , .26156734625504664    ,
      .17903431407718266       ,-.27447430552974527    ,
     -.66292226406569883e-01   ,
     -.1e-19                   , .2e-19                ,
     -.11e-18                  , .55e-18               ,
     -.288e-17                 , .1631e-16             ,
     -.10012e-15               , .67481e-15            ,
     -.506903e-14              , .4326596e-13          ,
     -.43045789e-12            , .516826239e-11        ,
     -.7864091377e-10          , .163064646352e-08     ,
     -.5170594537606e-07       , .307518478751947e-05  ,
     -.53652204681321174e-03   , .19989206986950373e01 ,
      .1e-19                   ,-.3e-19                ,
      .13e-18                  ,-.62e-18               ,
      .311e-17                 ,-.1669e-16             ,
      .9662e-16                ,-.60999e-15            ,
      .425523e-14              ,-.3336328e-13          ,
      .30061451e-12            ,-.320674742e-11        ,
      .4220121905e-10          ,-.72719159369e-09      ,
      .1797245724797e-07       ,-.74144984110606e-06   ,
     .683851994261165e-04      ,-.31111709210674018e-01};

 x = xd;
 y = abs(x);
 z = y * .125;

 if(z == 0.0) {
   vj0 = 1.0;
   vj0 =-xlg;
   return 0;
 }

 if(z > 1.0) {
   z = 1.0 / z;
   x2= 4.0 * z*z - 2.0;
   n1= 38;
   n2= 55;
   x1 = recurrence(n1, n2, x2, a);

   n1 = 56;
   n2 = 73;
   fx = recurrence(n1, n2, x2, a);
   x2 = cos(y - 0.7853981633974483);
   x3 = sin(y - 0.7853981633974483);
   x4 = 0.7978845608028654 / sqrt(y);
   fx = fx * z;
   vj0= x4 * (x1*x2 - fx*x3);
   vy0= x4 * (fx*x2 + x1*x3);
   return 0;
 }
 else {
   x2 = 4.0 * z*z - 2.0;
   n1 = 1;
   n2 = 18;
   vj0= recurrence(n1, n2, x2, a);

   if(n <= 0) return 0;

   n1 = 19;
   n2 = 37;

   fx = recurrence(n1, n2, x2, a);
   vy0= .6366197723675813 * log(y) * vj0 + fx;
   return 0;
 }
}

int hank11( double& vj1, double& vy1, double xd, int n ) {

 int    n1, n2;

 double x, y, z, fx, x1, x2, x3, x4;
 double xlg = 1.0e+70;
 double a[] = {0.0,
     -.4e-19                   , .295e-17              ,
     -.19554e-15               , .1138572e-13          ,
     -.57774042e-12            , .2528123664e-10       ,
     -.94242129816e-9          , .2949707007278e-7     ,
     -.76175878054003e-6       , .1588701923993213e-4  ,
     -.26044438934858068e-3    , .32402701826838575e-2 ,
     -.29175524806154208e-1    , .17770911723972828e0  ,
     -.66144393413454325e0     , .12879940988576776e1  ,
     -.11918011605412169e1     , .12967175412105298e1  ,
      .9e-19                   ,-.658e-17              ,
      .42773e-15               ,-.2440949e-13          ,
      .121143321e-11           ,-.5172121473e-10       ,
      .187547032473e-8         ,-.5688440039919e-7     ,
      .141662436449235e-5      ,-.283046401495148e-4   ,
      .44047862986709951e-3    ,-.51316411610610848e-2 ,
      .42319180353336904e-1    ,-.22662499155675492e0  ,
      .67561578077218767e0     ,-.76729636288664594e0  ,
     -.12869738438135000e0     , .40608211771868508e-1 ,
      .1e-19                   ,-.2e-19                ,
      .12e-18                  ,-.58e-18               ,
      .305e-17                 ,-.1731e-16             ,
      .10668e-15               ,-.72212e-15            ,
      .545267e-14              ,-.4684224e-13          ,
      .46991955e-12            ,-.570486364e-11        ,
      .881689866e-10           ,-.187189074911e-8      ,
      .6177633960644e-7        ,-.398728430048891e-5   ,
      .89898983308594085e-3    , .20018060817200274e1  ,
     -.1e-19                   , .3e-19                ,
     -.14e-18                  , .65e-18               ,
     -.328e-17                 , .1768e-16             ,
     -.10269e-15               , .65083e-15            ,
     -.456125e-14              , .3596777e-13          ,
     -.32643157e-12            , .351521879e-11        ,
     -.4686363688e-10          , .82291933277e-9       ,
     -.2095978138408e-7        , .91386152579555e-6    ,
     -.9627723549157079e-4     , .93555574139070650e-1};

 x = xd;
 y = abs(x);
 z = y * .125;

 if(z <= 0.0) {
   vj1 = 0.0;
   vj1 =-xlg;
   return 0;
 }

 if(z > 1) {
   z = 1.0 / z;
   x2= 4.0 * z*z - 2.0;
   n1= 37;
   n2= 54;
   x1= recurrence(n1, n2, x2, a);

   n1 = 55;
   n2 = 72;
   fx = recurrence(n1, n2, x2, a);
   x2 = cos(y - 2.356194490192345);
   x3 = sin(y - 2.356194490192345);
   x4 = 0.7978845608028654 / sqrt(y);
   fx = fx * z;
   vj1= x4 *(x1 * x2 - fx * x3);
   vy1= x4 *(fx * x2 + x1 * x3);
   return 0;
 }
 else {
   x2 = 4.0 * z*z - 2.0;
   n1 = 1;
   n2 = 18;
   fx = recurrence(n1, n2, x2, a);
   vj1= fx * z;
   if(n <= 0) return 0;

   n1 = 19;
   n2 = 36;
   fx = recurrence(n1, n2, x2, a);
   vy1= 0.6366197723675813 * (log(y) * vj1 - 1.0 / y) + fx * z;
   return 0;
 }
}

double recurrence( int n1, int n2, double x2, double* a ) {

 int i;

 double q1, q2 = 0.0, q3 = 0.0;

 for(i = n1; i <= n2; i++) {
    q1 = q2;
    q2 = q3;
    q3 = x2 * q2 - q1 + a[i];
 }
 
 return (q3 - q1) * 0.5;
}
