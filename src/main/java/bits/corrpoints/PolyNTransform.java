/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

/**
 * Arbitrary degree transform.
 * Due to looping, this will provide worse performance than fixed-size transforms.
 * Worst performance of all transforms.
 *
 * @author decamp
 */
public class PolyNTransform implements PolyTransform {

    private final int      mDegree;
    private final double[] mForward;
    private final double[] mBackward;


    public PolyNTransform( int degree, double[] forwardRef, double[] backwardRef ) {
        mDegree   = degree;
        mForward  = forwardRef;
        mBackward = backwardRef;
    }


    public boolean isApplicable() {
        return mForward != null;
    }


    public boolean isInvertible() {
        return mBackward != null;
    }


    public void apply( double x, double y, double[] out2x1, int outOff ) {
        if( mForward == null ) {
            throw new UnsupportedOperationException();
        }
        apply( mDegree, mForward, x, y, out2x1, outOff );
    }


    public void invert( double x, double y, double[] out2x1, int outOff ) {
        if( mBackward == null ) {
            throw new UnsupportedOperationException();
        }
        apply( mDegree, mBackward, x, y, out2x1, outOff );
    }


    public int degree() {
        return mDegree;
    }



    static int coeffCount( int degree ) {
        return (degree + 1) * (degree + 2);
    }


    static void apply( int degree, double[] coeffs, double x, double y, double[] out, int outOff ) {
        double yy  = 1.0;
        double u   = 0.0;
        double v   = 0.0;
        int    idx = 0;

        for( int py = 0; py <= degree; py++ ) {
            double xx = 1.0;
            for( int px = degree - py; px >= 0; px-- ) {
                u += xx * yy * coeffs[idx++];
                v += xx * yy * coeffs[idx++];
                xx *= x;
            }
            yy *= y;
        }

        out[outOff    ] = u;
        out[outOff + 1] = v;
    }

}
