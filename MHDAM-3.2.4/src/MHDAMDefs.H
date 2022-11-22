#ifndef  MHDAM_DEFS_H
#define  MHDAM_DEFS_H

#define  FOREVER             for (;;)
#define  THEN                {
#define  ELSEIF              } else {
#define  ENDIF               }
#define  ENDFOR              }
#define  ENDSWITCH           }

#ifndef  ABS
#define  ABS(x)              ((x)<0?-(x):(x))
#endif
#ifndef  MAX
#define  MAX(x,y)            (((x)>(y))?(x):(y))
#endif
#ifndef  MIN
#define  MIN(x,y)            (((x)<=(y))?(x):(y))
#endif
#ifndef  INRANGE
#define  INRANGE(a,b,c)      (((b)<(a))?(a):(((b)<(c))?(b):(c)))
#endif
#ifndef  MINMOD
#define  MINMOD(a,b)         ((a)>0.0?((b)<=0.0?0.0:((a)>(b)?(b):(a))):\
                                      ((b)>=0.0?0.0:((a)<(b)?(b):(a))))
#endif
#ifndef  DIV2
#define  DIV2(x)             (((unsigned)(x))>>1)
#endif
#ifndef  TIMES2
#define  TIMES2(x)           ((x)<<1)
#endif
#ifndef  TOHIBYTE
#define  TOHIBYTE(a)         ((a)<<8)
#endif
#ifndef  TOLOBYTE
#define  TOLOBYTE(a)         (((a)>>8)&0xff)
#endif
#ifndef  ISODD
#define  ISODD(a)            ((a) & 0x01)
#endif
#ifndef  ISEVEN
#define  ISEVEN(a)           (((a) % 2) == 0)
#endif
#ifndef  E
#define  E(i,j)              setw((i)) << setprecision((j))
#endif

                                        // Borland's style constants definitions
#ifndef  d_E
#define  d_E            2.71828182845904523536028747135266249    //            E
#endif
#ifndef  d_LOG2E
#define  d_LOG2E        1.44269504088896340735992468100189213    //    log2( E )
#endif
#ifndef  d_LOG10E
#define  d_LOG10E       0.434294481903251827651128918916605082   //      lg( E )
#endif
#ifndef  d_LN2
#define  d_LN2          0.693147180559945309417232121458176568   //      ln( 2 )
#endif
#ifndef  d_LN10
#define  d_LN10         2.30258509299404568401799145468436420    //     ln( 10 )
#endif
#ifndef  d_PI
#define  d_PI           3.14159265358979323846264338327950288    //           PI
#endif
#ifndef  d_2PI
#define  d_2PI          6.28318530717958647692528676655900576    //         2*PI
#endif
#ifndef  d_3PI
#define  d_3PI          9.42477796076937971538793014983850865    //         3*PI
#endif
#ifndef  d_PI_2
#define  d_PI_2         1.57079632679489661923132169163975144    //         PI/2
#endif
#ifndef  d_3PI_2
#define  d_3PI_2        4.71238898038468985769396507491925432    //       3*PI/2
#endif
#ifndef  d_PI_3
#define  d_PI_3         1.04719755119659774615421446109316762    //         PI/3
#endif
#ifndef  d_PI_4
#define  d_PI_4         0.785398163397448309615660845819875721   //         PI/4
#endif
#ifndef  d_3PI_4
#define  d_3PI_4        2.35619449019234492884698253745962716    //       3*PI/4
#endif
#ifndef  d_PI_180
#define  d_PI_180       0.0174532925199432957692369076848861271  //       PI/180
#endif
#ifndef  d_1_PI
#define  d_1_PI         0.318309886183790671537767526745028724   //         1/PI
#endif
#ifndef  d_2_PI
#define  d_2_PI         0.636619772367581343075535053490057448   //         2/PI
#endif
#ifndef  d_1_4PI
#define  d_1_4PI        0.079577471545947667884441881686257181   //     1/(4*PI)
#endif
#ifndef  d_1_8PI
#define  d_1_8PI        0.039788735772973833942220940843128591   //     1/(4*PI)
#endif
#ifndef  d_180_PI
#define  d_180_PI       57.2957795130823208767981548141051703    //       180/PI
#endif
#ifndef  d_SQRTPI
#define  d_SQRTPI       1.77245385090551602729816748334114518    //   sqrt( PI )
#endif
#ifndef  d_1_SQRTPI
#define  d_1_SQRTPI     0.564189583547756286948079451560772586   // 1/sqrt( PI )
#endif
#ifndef  d_2_SQRTPI
#define  d_2_SQRTPI     1.12837916709551257389615890312154517    // 2/sqrt( PI )
#endif
#ifndef  d_LNPI
#define  d_LNPI         1.14472988584940017414342735135305871    //     ln( PI )
#endif
#ifndef  d_SQRT2
#define  d_SQRT2        1.41421356237309504880168872420969808    //    sqrt( 2 )
#endif
#ifndef  d_SQRT_2
#define  d_SQRT_2       0.707106781186547524400844362104849039   //  sqrt( 2 )/2
#endif
#ifndef  d_GAMMA
#define  d_GAMMA        0.577215664901532860606512090082402431   //        gamma
#endif

#endif // MHDAM_DEFS_H
