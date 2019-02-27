/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

/**
 * 2D projective transform.
 *
 * @author decamp
 */
public class ProjectiveTransform implements Transform22 {

    private final double[] mForward;
    private final double[] mBackward;


    public ProjectiveTransform( double[] forwardMatRef, double[] backwardMatRef ) {
        mForward  = forwardMatRef;
        mBackward = backwardMatRef;
    }


    private static void apply( double[] coeffs, double x, double y, double[] out2x1, int outOff ) {
        if( coeffs == null ) {
            throw new UnsupportedOperationException();
        }

        double a = coeffs[0] * x + coeffs[1] * y + coeffs[2];
        double b = coeffs[3] * x + coeffs[4] * y + coeffs[5];
        double c = coeffs[6] * x + coeffs[7] * y + coeffs[8];

        out2x1[outOff    ] = a / c;
        out2x1[outOff + 1] = b / c;
    }


    public boolean isApplicable() {
        return mForward != null;
    }


    public boolean isInvertible() {
        return mBackward != null;
    }


    public void apply( double x, double y, double[] out2x1, int outOff ) {
        apply( mForward, x, y, out2x1, outOff );
    }


    public void invert( double x, double y, double[] out2x1, int outOff ) {
        apply( mBackward, x, y, out2x1, outOff );
    }

    /**
     * @return 3x3 column-major homography matrix that performs backward mapping.
     */
    public double[] getBackwardMatrix() {
        return coeffsToMatrix( mForward );
    }

    /**
     * @return 3x3 column-major homography matrix that performs forward mapping.
     */
    public double[] getForward3x3Matrix() {
        return coeffsToMatrix( mForward );
    }


    /**
     * For debugging.
     */
    double[] coeffsRef() {
        return mForward;
    }


    double[] inverseCoeffsRef() {
        return mBackward;
    }


    private static double[] coeffsToMatrix( double[] coeffs ) {
        return new double[]{
            coeffs[0],
            coeffs[3],
            coeffs[6],
            coeffs[1],
            coeffs[4],
            coeffs[7],
            coeffs[2],
            coeffs[5],
            coeffs[8]
        };
    }


}
