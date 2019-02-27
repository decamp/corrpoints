/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

/**
 * 2nd Degree polynomial transform.
 * 
 * @author decamp
 */
public class Poly2Transform implements PolyTransform {
    
    
    private final double[] mForward;
    private final double[] mBackward;
    
    
    public Poly2Transform( double[] forward6x2Ref, double[] backward6x2Ref ) {
        mForward  = forward6x2Ref;
        mBackward = backward6x2Ref;
    }
    
    
    public boolean isApplicable() {
        return mForward != null;
    }


    public boolean isInvertible() {
        return mBackward != null;
    }    


    public void apply(double x, double y, double[] out2x1, int outOff) {
        apply(mForward, x, y, out2x1, outOff);
    }


    public void invert(double x, double y, double[] out2x1, int outOff) {
        apply(mBackward, x, y, out2x1, outOff);
    }


    public int degree() {
        return 2;
    }
    
    
    /**
     * Package private debugging.
     */
    double[] coeffsRef() {
        return mForward;
    }


    private static void apply( double[] coeffs, double x, double y, double[] out, int outOff ) {
        if( coeffs == null ) {
            throw new UnsupportedOperationException();
        }
        
        double v;
        
        v  = coeffs[ 0];
        v += coeffs[ 1] * x;
        v += coeffs[ 2] * x * x;
        v += coeffs[ 3] * y;
        v += coeffs[ 4] * y * x;
        v += coeffs[ 5] * y * y;
        out[outOff  ] = v;
        
        v  = coeffs[ 6];
        v += coeffs[ 7] * x;
        v += coeffs[ 8] * x * x;
        v += coeffs[ 9] * y;
        v += coeffs[10] * y * x;
        v += coeffs[11] * y * y;
        out[outOff+1] = v;
    }

}
