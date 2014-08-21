package bits.corrpoints;

/**
 * Linear or 1st degree polynomial transform.
 *
 * @author decamp
 */
public class LinearTransform implements PolyTransform {

    private final double mTx;
    private final double mSxx;
    private final double mSyx;
    private final double mTy;
    private final double mSxy;
    private final double mSyy;
    private final double mDet;

    public LinearTransform( double tx,
                            double sxx,
                            double syx,
                            double ty,
                            double sxy,
                            double syy )
    {
        mTx  = tx;
        mSxx = sxx;
        mSyx = syx;
        mTy  = ty;
        mSxy = sxy;
        mSyy = syy;
        double det = sxx * syy - sxy * syx;
        mDet = Math.abs( det ) > TransformFitter.MACHINE_EPSILON_DOUBLE ? det : 0.0;
    }


    public boolean isApplicable() {
        return true;
    }

    public boolean isInvertible() {
        return mDet != 0.0;
    }

    public void apply( double x, double y, double[] out2x1, int outOff ) {
        out2x1[outOff] = mSxx * x + mSyx * y + mTx;
        out2x1[outOff + 1] = mSxy * x + mSyy * y + mTy;
    }

    public void invert( double x, double y, double[] out2x1, int outOff ) {
        if( mDet == 0.0 ) {
            throw new UnsupportedOperationException();
        }

        x -= mTx;
        y -= mTy;

        out2x1[outOff    ] = ( mSyy * x - mSyx * y ) / mDet;
        out2x1[outOff + 1] = (-mSxy * x + mSxx * y ) / mDet;
    }

    public int degree() {
        return 1;
    }

}
