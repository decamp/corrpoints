/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

/**
 * 3rd Degree polynomial transform.
 * 
 * @author decamp
 */
public class Poly3Transform implements PolyTransform {
    
    
    private final double[] mForward;
    private final double[] mBackward;
    
    
    public Poly3Transform( double[] forward10x2Ref, double[] backward10x2Ref ) {
        mForward  = forward10x2Ref;
        mBackward = backward10x2Ref;
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
        return 3;
    }
    
    
    /**
     * Package private debugging.
     */
    double[] coeffsRef() {
        return mForward;
    }
    
    private static void apply(double[] coeffs, double x, double y, double[] out, int outOff) {
        if(coeffs == null)
            throw new UnsupportedOperationException();
        
        double v;
        
        v  = coeffs[ 0];
        v += coeffs[ 1] * x;
        v += coeffs[ 2] * y;
        v += coeffs[ 3] * x * y;
        v += coeffs[ 4] * x * x;
        v += coeffs[ 5] * y * y;
        v += coeffs[ 6] * x * x * y;
        v += coeffs[ 7] * y * y * x;
        v += coeffs[ 8] * x * x * x;
        v += coeffs[ 9] * y * y * y;
        out[outOff  ] = v;
        
        v  = coeffs[10];
        v += coeffs[11] * x;
        v += coeffs[12] * y;
        v += coeffs[13] * x * y;
        v += coeffs[14] * x * x;
        v += coeffs[15] * y * y;
        v += coeffs[16] * x * x * y;
        v += coeffs[17] * y * y * x;
        v += coeffs[18] * x * x * x;
        v += coeffs[19] * y * y * y;
        out[outOff+1] = v;
    }

}
