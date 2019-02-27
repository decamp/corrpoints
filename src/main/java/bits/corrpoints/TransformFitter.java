/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

import org.netlib.lapack.*;
import org.netlib.util.intW;


/**
 * Computes transform functions from correspondence points.
 *
 * @author decamp
 */
public class TransformFitter {
    
    public static final int ERR_OK                = 0;
    public static final int ERR_SINGULAR          = 1;
    public static final int ERR_DID_NOT_CONVERGE  = 2;
    public static final int ERR_ILLEGAL_VALUE     = 3;
    public static final int ERR_INSUFFICIENT_RANK = 4;
    
    static final double MACHINE_EPSILON_DOUBLE     = 1.1102230246251565e-16;
    static final double LOG_2                      = Math.log( 2.0 );
    

    /**
     * Attempts to find a projective transform that maps a set of 2D 
     * <i>fromPoints</i> to a set of 2D <i>toPoints</i>, referred to 
     * together as the <i>correspondence points</i>.  
     * <p>
     * At least 4 correspondence points are needed to define a projective
     * transform, although in some cases more points may be required if the 
     * points define a degenerate mapping (eg, providing the same correspondence point
     * four times). 
     * <p>
     * If more correspondence points are provided, the transform may be overdefined,
     * in which case the transform will attempt to find the mapping that minimizes
     * the least-squares error.
     * 
     * @param fromPoints      Array countaining set of points: (x0, y0, x1, y1, ... xn, yn). 
     * @param fromPointsOff   Offset into fromPoints array where data begins.
     * @param toPoints        Array containing set of points: (x0, y0, x1, y1, ... xn, yn)
     * @param toPointsOff     Offset into toPoints array where data begins.
     * @param nPoints         Number of correspondence points used for training.  There must be at least 4 points.
     * @param errorOut        <code>errorOut[0]</code> will hold the resulting error code on return.  May be <code>null</code>.
     * @return a ProjectiveTransform object that optimally maps the provided fromPoints to toPoints.  On failure, returns <code>null</code>
     */
    public static ProjectiveTransform corrPointsToProjectiveTransform(
        double[] fromPoints,
        int fromPointsOff,
        double[] toPoints,
        int toPointsOff,
        int nPoints,
        int[] errorOut
    ) {
        //For projection:
        // u = (Ax + By + C) / (Gx + Hy + I);
        // v = (Dx + Ey + F) / (Gx + Hy + I);
        
        //Assume I = 1.0 and multiply both equations by denominator.
        // u = [x y 1 0 0 0 -ux -uy] * [A B C D E F G H]'
        // v = [0 0 0 x y 1 -vx -vy] * [A B C D E F G H]'

        // With 4 or more correspondence points we can combine the u equations and 
        // the v equations for one linear system to solve for [A B C D E F G H]: 
        // 
        // [ u1  ] = [ x1  y1  1  0   0   0  -u1*x1  -u1*y1 ] * [A] 
        // [ u2  ] = [ x2  y2  1  0   0   0  -u2*x2  -u2*y2 ]   [B] 
        // [ u3  ] = [ x3  y3  1  0   0   0  -u3*x3  -u3*y3 ]   [C] 
        // [ u4  ] = [ x4  y4  1  0   0   0  -u4*x4  -u4*y4 ]   [D]
        // [ ... ]   [ ...                                  ]   [E] 
        // [ un  ] = [ xn  yn  1  0   0   0  -un*xn  -un*yn ]   [F] 
        // [ v1  ] = [ 0   0   0  x1  y1  1  -v1*x1  -v1*y1 ]   [G] 
        // [ v2  ] = [ 0   0   0  x2  y2  1  -v2*x2  -v2*y2 ]   [H] 
        // [ v3  ] = [ 0   0   0  x3  y3  1  -v3*x3  -v3*y3 ] 
        // [ v4  ] = [ 0   0   0  x4  y4  1  -v4*x4  -v4*y4 ] 
        // [ ... ]   [ ...                                  ]   
        // [ vn  ] = [ 0   0   0  xn  yn  1  -vn*xn  -vn*yn ] 
        // 
        // Or rewriting the above matrix equation: 
        // U = X * Tvec, where Tvec = [A B C D E F G H]' 
        // so Tvec = X\U. 
        // 
        
        if( nPoints < 4 ) {
            if( errorOut != null ) {
                errorOut[0] = ERR_INSUFFICIENT_RANK;
            }
            return null;
        }

        //Main matrix dimensions.
        final int m     = nPoints * 2;
        final int n     = 8;
        final int dim0  = Math.min(m, n);
        final int dim1  = Math.max(m, n);
        
        final double[] aMat = new double[m*n];
        final double[] bMat = new double[dim1];
        
        for(int i = 0; i < nPoints; i++) {
            double xx = fromPoints[i*2   + fromPointsOff];
            double yy = fromPoints[i*2+1 + fromPointsOff];
            double uu = toPoints  [i*2   + toPointsOff  ];
            double vv = toPoints  [i*2+1 + toPointsOff  ];

            aMat[i+m*0] = xx;
            aMat[i+m*1] = yy;
            aMat[i+m*2] = 1.0;
            aMat[i+m*6] = -uu*xx;
            aMat[i+m*7] = -uu*yy;

            aMat[i+m*3 + nPoints] = xx; 
            aMat[i+m*4 + nPoints] = yy;
            aMat[i+m*5 + nPoints] = 1.0;
            aMat[i+m*6 + nPoints] = -vv*xx;
            aMat[i+m*7 + nPoints] = -vv*yy;
            
            bMat[i  ]       = toPoints[i*2   + toPointsOff];
            bMat[i+nPoints] = toPoints[i*2+1 + toPointsOff];
        }
        
        //Tedious computations of required workspace sizes.
        final int smlsiz   = Ilaenv.ilaenv(9, "DGELSD", "", 0, 0, 0, 0);
        final int nrhs     = 1;
        final int nlvl     = Math.max( 0, (int)( Math.log(dim0 / (smlsiz+1.0)) / LOG_2 )) + 1;
        final int lwork    = 12 * dim0 + 2 * dim0 * smlsiz + 8 * dim0 * nlvl + dim0 * nrhs + (smlsiz + 1) * (smlsiz + 1);
        final int liwork   = 3 * dim0 * nlvl + 11 * dim0;
        
        //Allocate stuff.
        final double[] sMat    = new double[9];
        final double rcond     = MACHINE_EPSILON_DOUBLE;
        final intW rankBox     = new intW(0);
        final double[] workMat = new double[lwork];
        final int[] iWork      = new int[liwork];
        final intW infoBox     = new intW(0);

        //Solve for least squares using Divide-and-Conquer SVD.
        Dgelsd.dgelsd(
            m,
            n,
            nrhs,
            aMat,
            0,
            m,
            bMat,
            0,
            dim1,
            sMat,
            0,
            rcond,
            rankBox,
            workMat,
            0,
            lwork,
            iWork,
            0,
            infoBox
        );
        
        if( infoBox.val != 0 || rankBox.val < 8 ) {
            if( errorOut != null ) {
                if( infoBox.val == 0 ) {
                    errorOut[0] = ERR_INSUFFICIENT_RANK;
                } else {
                    errorOut[0] = infoBox.val < 0 ? ERR_ILLEGAL_VALUE : ERR_DID_NOT_CONVERGE;
                }
            }
            return null;
        }
        
        double[] forwardMat = new double[9];
        System.arraycopy( bMat, 0, forwardMat, 0, Math.min( 8, m ) );
        
        //We always assumed that W = 1.0.
        forwardMat[8] = 1.0;

        double[] backMat = forwardMat.clone();
        boolean hasBack  = false;
        
        //Invert coefficient matrix.
        Dgetrf.dgetrf( 3, 3, backMat, 0, 3, iWork, 0, infoBox );
        if( infoBox.val == 0 ) {
            Dgetri.dgetri( 3, backMat, 0, 3, iWork, 0, workMat, 0, lwork, infoBox );
            hasBack = infoBox.val == 0;
        }

        if( errorOut != null ) {
            errorOut[0] = ERR_OK;
        }
        
        return new ProjectiveTransform( forwardMat, hasBack ? backMat : null );
    }
    

    /**
     * Attempts to find a 2nd-degree polynomial transform that maps a set of 2D 
     * <i>fromPoints</i> to a set of 2D <i>toPoints</i>, referred to 
     * together as the <i>correspondence points</i>.  
     * <p>
     * At least 6 correspondence points are needed to define a 2nd degree polynomial
     * transform, although in some cases more points may be required if the 
     * points define a degenerate mapping.
     * <p>
     * If more correspondence points are provided, the transform may be overdefined,
     * in which case the transform will attempt to find the mapping that minimizes
     * the least-squares error.
     * 
     * @param fromPoints     Array countaining set of points: (x0, y0, x1, y1, ... xn, yn). 
     * @param fromPointsOff  Offset into fromPoints array where data begins.
     * @param toPoints       Array containing set of points: (x0, y0, x1, y1, ... xn, yn)
     * @param toPointsOff    Offset into toPoints array where data begins.
     * @param nPoints        Number of correspondence points used for training.  There must be at least 6 points.
     * @param errorOut       <code>errorOut[0]</code> will hold the resulting error code on return.  May be <code>null</code>.
     * @return a Poly2Transform object that optimally maps the provided fromPoints to toPoints.  On failure, returns <code>null</code>
     */
    public static Poly2Transform corrPointsToPoly2Transform(
        double[] fromPoints,
        int fromPointsOff,
        double[] toPoints,
        int toPointsOff,
        int nPoints,
        int[] errorOut
    ) {
        //For 3rd degree polynomial:
        // u = A0 + B0x + C0y + D0xy + E0xx + F0yy
        // v = A1 + B1x + C1y + D1xy + E1xx + F1yy
        //
        //Separate knows from unknows.
        // [u0 v0] = [1.0 x0 y0 x0*y0 x0*x0 y0y0] * [A0 B0 C0 D0 E0 F0]'
        // [u1 v1] = [1.0 x1 y1 x1*y1 x1*x1 y1y1]   [A1 B1 C1 D1 E1 F1]
        //   ...   =    ...                                                             
        // [un vn] = [1.0 xn yn xn*yn xn*xn ynyn]
        //
        //Restated:
        // postPoints_nx2 = xMat_nx6 * coeffs_6x2
        //
        //Coefficients can then be solved using least squares.
        // coeffs6x2 = xMat_nx6 \ postPoints_nx2
        
        final int k = 6; //Coefficient count. 

        if( nPoints < k ) {
            if( errorOut != null ) {
                errorOut[0] = ERR_INSUFFICIENT_RANK;
            }
            return null;
        }


        final int m    = nPoints;
        final int n    = k;
        final int dim0 = Math.min( m, n );

        final double[] xMat = new double[m * n];
        final double[] uMat = new double[m * 2];
        
        
        for(int i = 0; i < m; i++) {
            double xx = fromPoints[i*2   + fromPointsOff];
            double yy = fromPoints[i*2+1 + fromPointsOff];
            
            xMat[i+m*0] = 1.0;
            xMat[i+m*1] = xx;
            xMat[i+m*2] = xx * xx;
            xMat[i+m*3] = yy;
            xMat[i+m*4] = yy * xx;
            xMat[i+m*5] = yy * yy;
            
            uMat[i  ] = toPoints[i*2   + toPointsOff];
            uMat[i+m] = toPoints[i*2+1 + toPointsOff];
        }

        //Tedious workspace allocations.
        final int smlsiz = Ilaenv.ilaenv(9, "DGELSD", "", 0, 0, 0, 0);
        final int nrhs   = 2;
        final int nlvl   = Math.max( 0, (int)( Math.log(dim0 / (smlsiz+1.0)) / LOG_2 )) + 1;
        final int lwork  = 12 * dim0 + 2 * dim0 * smlsiz + 8 * dim0 * nlvl + dim0 * nrhs + (smlsiz + 1) * (smlsiz + 1);
        final int liwork = 3 * dim0 * nlvl + 11 * dim0;
        
        //Allocate stuff.
        final double[] sMat    = new double[dim0];
        final double rcond     = MACHINE_EPSILON_DOUBLE;
        final intW rankBox     = new intW(0);
        final double[] workMat = new double[lwork];
        final int[] iWork      = new int[liwork];
        final intW infoBox     = new intW(0);
        
        //Solve Least squares using divide-and-conquer SVD.
        //Output will be stored in uMat.
        Dgelsd.dgelsd(
            m,
            n,
            nrhs,
            xMat,
            0,
            m,
            uMat,
            0,
            m,
            sMat,
            0,
            rcond,
            rankBox,
            workMat,
            0,
            workMat.length,
            iWork,
            0,
            infoBox
        );
        
        //Check if error or if input data had insufficient rank.
        if( infoBox.val != 0 || rankBox.val < k ) {
            if( errorOut != null ) {
                if( infoBox.val == 0 ) {
                    errorOut[0] = ERR_INSUFFICIENT_RANK;
                } else {
                    errorOut[0] = infoBox.val < 0 ? ERR_ILLEGAL_VALUE : ERR_DID_NOT_CONVERGE;
                }
            }
            
            return null;
        }
        
        //Place coeffs into different matrix, if necessary.
        double[] coeffs;
        if( m == k ) {
            coeffs = uMat;
        } else {
            coeffs = new double[k * 2];
            System.arraycopy( uMat, 0, coeffs, 0, k );
            System.arraycopy( uMat, m, coeffs, k, k );
        }
        
        //System.out.println(MatFormat.formatToScreen( m, 2, uMat, 0, m ));
        return new Poly2Transform( coeffs, null );
    }
    
    
    /**
     * Attempts to find a 3rd-degree polynomial transform that maps a set of 2D 
     * <i>fromPoints</i> to a set of 2D <i>toPoints</i>, referred to 
     * together as the <i>correspondence points</i>.  
     * <p>
     * At least 10 correspondence points are needed to define a 3rd degree polynomial
     * transform, although in some cases more points may be required if the 
     * points define a degenerate mapping.
     * <p>
     * If more correspondence points are provided, the transform may be overdefined,
     * in which case the transform will attempt to find the mapping that minimizes
     * the least-squares error.
     * 
     * @param fromPoints      Array countaining set of points: (x0, y0, x1, y1, ... xn, yn). 
     * @param fromPointsOff   Offset into fromPoints array where data begins.
     * @param toPoints        Array containing set of points: (x0, y0, x1, y1, ... xn, yn)
     * @param toPointsOff     Offset into toPoints array where data begins.
     * @param nPoints         Number of correspondence points used for training.  There must be at least 4 points.
     * @param errorOut        <code>errorOut[0]</code> will hold the resulting error code on return.  May be <code>null</code>.
     * @return a Poly3Transform object that optimally maps the provided fromPoints to toPoints.  On failure, returns <code>null</code>
     */
    public static Poly3Transform corrPointsToPoly3Transform(
        double[] fromPoints,
        int fromPointsOff,
        double[] toPoints,
        int toPointsOff,
        int nPoints,
        int[] errorOut
    ) {
        // For 3rd degree polynomial:
        //  u = A0 + B0x + C0y + D0xy + E0xx + F0yy +G0xxy + H0yyx + I0xxx + J0yyy
        //  v = A1 + B1x + C1y + D1xy + E1xx + F1yy +G1xxy + H1yyx + I1xxx + J1yyy
        //
        // Separate knows from unknowns.
        //  [u0 v0] = [1.0 x0 y0 x0*y0 x0*x0 y0*y0 x0*x0*y0 y0*y0*x0 x0*x0*x0 y0*y0*y0] * [A0 B0 C0 D0 E0 F0 G0 H0 I0 J0]'
        //  [u1 v1] = [1.0 x1 y1 x1*y1 x1*x1 y1*y1 x1*x1*y1 y1*y1*x1 x1*x1*x1 y1*y1*y1]   [A1 B1 C1 D1 E1 F1 G1 H1 I1 J1]
        //    ...   =    ...
        //  [un vn] = [1.0 xn yn xn*yn xn*xn yn*yn xn*xn*yn yn*yn*xn xn*xn*xn yn*yn*yn]
        //
        // Restated:
        //  postPoints_nx2 = xMat_nx10 * coeffs_10x2
        //
        // Coefficients can then be solved using least squares.
        //  coeffs10x2 = xMat_nx10 \ postPoints_nx2
        
        final int k = 10; //Coefficient count.
        if( nPoints < k ) {
            if( errorOut != null ) {
                errorOut[0] = ERR_INSUFFICIENT_RANK;
            }
            return null;
        }


        final int m    = nPoints;
        final int n    = 10;
        final int dim0 = Math.min(m, n);

        final double[] xMat = new double[m * n];
        final double[] uMat = new double[m * 2];
        
        
        for(int i = 0; i < m; i++) {
            double xx = fromPoints[i*2   + fromPointsOff];
            double yy = fromPoints[i*2+1 + fromPointsOff];
            
            xMat[i+m*0] = 1.0;
            xMat[i+m*1] = xx;
            xMat[i+m*2] = yy;
            xMat[i+m*3] = xx * yy;
            xMat[i+m*4] = xx * xx;
            xMat[i+m*5] = yy * yy;
            xMat[i+m*6] = xx * xx * yy;
            xMat[i+m*7] = yy * yy * xx;
            xMat[i+m*8] = xx * xx * xx;
            xMat[i+m*9] = yy * yy * yy;
            
            uMat[i  ] = toPoints[i*2   + toPointsOff];
            uMat[i+m] = toPoints[i*2+1 + toPointsOff];
        }
        
        //Tedious workspace allocations.
        final int smlsiz = Ilaenv.ilaenv( 9, "DGELSD", "", 0, 0, 0, 0 );
        final int nrhs   = 2;
        final int nlvl   = Math.max( 0, (int)( Math.log(dim0 / (smlsiz+1.0)) / LOG_2 )) + 1;
        final int lwork  = 12 * dim0 + 2 * dim0 * smlsiz + 8 * dim0 * nlvl + dim0 * nrhs + (smlsiz + 1) * (smlsiz + 1);
        final int liwork = 3 * dim0 * nlvl + 11 * dim0;
        
        //Allocate stuff.
        final double[] sMat    = new double[dim0];
        final double rcond     = MACHINE_EPSILON_DOUBLE;
        final intW rankBox     = new intW(0);
        final double[] workMat = new double[lwork];
        final int[] iWork      = new int[liwork];
        final intW infoBox     = new intW(0);
        
        //Solve Least squares using divide-and-conquer SVD.
        //Output will be stored in uMat.
        Dgelsd.dgelsd(
            m,
            n,
            nrhs,
            xMat,
            0,
            m,
            uMat,
            0,
            m,
            sMat,
            0,
            rcond,
            rankBox,
            workMat,
            0,
            workMat.length,
            iWork,
            0,
            infoBox
        );
        
        //Check if error or if input data had insufficient rank.
        if( infoBox.val != 0 || rankBox.val < k ) {
            if( errorOut != null ) {
                if( infoBox.val == 0 ) {
                    errorOut[0] = ERR_INSUFFICIENT_RANK;
                } else {
                    errorOut[0] = infoBox.val < 0 ? ERR_ILLEGAL_VALUE : ERR_DID_NOT_CONVERGE;
                }
            }
            return null;
        }

        //Place coeffs into different matrix, if necessary.
        double[] coeffs;
        if( m == k ) {
            coeffs = uMat;
        } else {
            coeffs = new double[20];
            System.arraycopy( uMat, 0, coeffs, 0, 10 );
            System.arraycopy( uMat, m, coeffs, 10, 10 );
        }

        //System.out.println(MatFormat.formatToScreen( m, 2, uMat, 0, m ));
        return new Poly3Transform( coeffs, null );
    }
    

    /**
     * Attempts to find a N-degree polynomial transform that maps a set of 2D 
     * <i>fromPoints</i> to a set of 2D <i>toPoints</i>, referred to 
     * together as the <i>correspondence points</i>.  
     * <p>
     * At least (N+2)(N+1) correspondence points are needed to define a N-degree polynomial
     * transform, although in some cases more points may be required if the 
     * points define a degenerate mapping.
     * <p>
     * If more correspondence points are provided, the transform may be overdefined,
     * in which case the transform will attempt to find the mapping that minimizes
     * the least-squares error.
     * 
     * @param fromPoints     Array countaining set of points: (x0, y0, x1, y1, ... xn, yn). 
     * @param fromPointsOff  Offset into fromPoints array where data begins.
     * @param toPoints       Array containing set of points: (x0, y0, x1, y1, ... xn, yn)
     * @param toPointsOff    Offset into toPoints array where data begins.
     * @param nPoints        Number of correspondence points used for training.  There must be at least 6 points.
     * @param errorOut       <code>errorOut[0]</code> will hold the resulting error code on return.  May be <code>null</code>.
     * @return a PolyTransform object that optimally maps the provided fromPoints to toPoints.  On failure, returns <code>null</code>
     */
    public static PolyTransform corrPointsToPolyTransform(
        double[] fromPoints,
        int fromPointsOff,
        double[] toPoints,
        int toPointsOff,
        int nPoints,
        int degree,
        int[] errorOut
    ) {
        //For N-degree polynomial:
        // u = A0 + B0x + C0x^2 ... K0y^n
        // v = A1 + B1x + C1x^2 ... K1y^n
        // k = number of coefficients / 2 = (N+1)(N+2) / 2
        //
        //Separate knows from unknows.
        // [u0 v0] = [1.0 x0 x0*x0 ...] * [A0 B0 C0 ...]'
        // [u1 v1] = [1.0 x1 x1*x1 ...]   [A1 B1 C1 ...]
        //   ...   =    ...                                                             
        // [un vn] = [1.0 xn xn*xn ...]
        //
        //Restated:
        // postPoints_nx2 = xMat_nxk * coeffs_kx2
        //
        //Coefficients can then be solved using least squares.
        // coeffs_kx2 = xMat_nxk \ postPoints_nx2

        if( degree < 1 ) {
            if( errorOut != null ) {
                errorOut[0] = ERR_ILLEGAL_VALUE;
            }
            return null;
        }

        final int k = (degree + 1) * (degree + 2) / 2; //Coefficient count / 2 

        if( nPoints < k ) {
            if( errorOut != null ) {
                errorOut[0] = ERR_INSUFFICIENT_RANK;
            }
            return null;
        }
        
        final int m    = nPoints;
        final int n    = k;
        final int dim0 = Math.min( m, n );

        final double[] xMat = new double[ m * n ];
        final double[] uMat = new double[ m * 2 ];


        for( int i = 0; i < m; i++ ) {
            double x   = fromPoints[i * 2 + fromPointsOff];
            double y   = fromPoints[i * 2 + 1 + fromPointsOff];
            double yy  = 1.0;
            int    col = 0;

            for( int py = 0; py <= degree; py++ ) {
                double xx = 1.0;
                for( int px = degree - py; px >= 0; px-- ) {
                    xMat[i + m * col++] = xx * yy;
                    xx *= x;
                }
                yy *= y;
            }

            uMat[i    ] = toPoints[i * 2 +     toPointsOff];
            uMat[i + m] = toPoints[i * 2 + 1 + toPointsOff];
        }
        
        //Tedious workspace allocations.
        final int smlsiz = Ilaenv.ilaenv( 9, "DGELSD", "", 0, 0, 0, 0 );
        final int nrhs   = 2;
        final int nlvl   = Math.max( 0, (int)(Math.log( dim0 / (smlsiz + 1.0) ) / LOG_2) ) + 1;
        final int lwork  = 12 * dim0 + 2 * dim0 * smlsiz + 8 * dim0 * nlvl + dim0 * nrhs + (smlsiz + 1) * (smlsiz + 1);
        final int liwork = 3 * dim0 * nlvl + 11 * dim0;
        
        //Allocate stuff.
        final double[] sMat    = new double[ dim0 ];
        final double rcond     = MACHINE_EPSILON_DOUBLE;
        final intW rankBox     = new intW( 0 );
        final double[] workMat = new double[ lwork ];
        final int[] iWork      = new int[ liwork ];
        final intW infoBox     = new intW( 0 );

        //Solve Least squares using divide-and-conquer SVD.
        //Output will be stored in uMat.
        Dgelsd.dgelsd(
            m,
            n,
            nrhs,
            xMat,
            0,
            m,
            uMat,
            0,
            m,
            sMat,
            0,
            rcond,
            rankBox,
            workMat,
            0,
            workMat.length,
            iWork,
            0,
            infoBox
        );
        
        //Check if error or if input data had insufficient rank.
        if( infoBox.val != 0 || rankBox.val < k ) {
            if( errorOut != null ) {
                if( infoBox.val == 0 ) {
                    errorOut[0] = ERR_INSUFFICIENT_RANK;
                } else {
                    errorOut[0] = infoBox.val < 0 ? ERR_ILLEGAL_VALUE : ERR_DID_NOT_CONVERGE;
                }
            }
            return null;
        }
        
        
        switch( degree ) {
        case 1:
        {
            return new LinearTransform(
                uMat[0  ],
                uMat[1  ],
                uMat[2  ],
                uMat[m  ],
                uMat[m+1],
                uMat[m+2]
            );
        }
        
        case 2:
        {
            double[] cf;
            
            if( m == k ) {
                cf = uMat;
            } else {
                cf = new double[k*2];
                System.arraycopy( uMat, 0, cf, 0, k );
                System.arraycopy( uMat, m, cf, k, k );
            }
            
            return new Poly2Transform( cf, null );
        }
        
        case 3:
        {
            // This coefficient reordering is a little ridiculous,
            // but I didn't want to make any changes to Poly3Transform.
            double[] cf ={
                uMat[0  ],
                uMat[1  ],
                uMat[4  ],
                uMat[5  ],
                uMat[2  ],
                uMat[7  ],
                uMat[6  ],
                uMat[8  ],
                uMat[3  ],
                uMat[9  ],
                uMat[  m],
                uMat[1+m],
                uMat[4+m],
                uMat[5+m],
                uMat[2+m],
                uMat[7+m],
                uMat[6+m],
                uMat[8+m],
                uMat[3+m],
                uMat[9+m]
            };
            
            return new Poly3Transform( cf, null );
        }
        
        default:
        {
            //Place coeffs into different matrix and transpose.
            double[] cf = new double[k*2];
            for( int i = 0; i < k; i++ ) {
                cf[i*2  ] = uMat[i  ];
                cf[i*2+1] = uMat[i+m];
            }
            return new PolyNTransform( degree, cf, null );
        }}
        
    }

}
