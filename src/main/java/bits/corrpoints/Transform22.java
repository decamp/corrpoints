package bits.corrpoints;

/**
 * A Transform22 maps 2D input points to 2D output points.
 * A given instance may also be able to invert the mapping.
 *
 * @author decamp
 */
public interface Transform22 {

    /**
     * @return true iff this Transform22 has a forward transform and supports calls to <code>apply()</code>.
     */
    public boolean isApplicable();

    /**
     * @return true iff this Transform22 has an inverse transform, and supports calls to <code>invert()</code>.
     */
    public boolean isInvertible();

    /**
     * Maps a 2D point to a 2D point.
     *
     * @param x      Input x-coordinate
     * @param y      Input y-coordinate
     * @param out2x1 Array where output is stored.
     * @param outOff Position in output array where values will be stored.
     * @throws UnsupportedOperationException if <code>!isApplicable()</code>.
     */
    public void apply( double x, double y, double[] out2x1, int outOff );

    /**
     * Maps a 2D point to a 2D point inversely to that of <code>apply()</code>.
     *
     * @param x      Input x-coordinate
     * @param y      Input y-coordinate
     * @param out2x1 Array where output is stored.
     * @param outOff Position in output array where values will be stored.
     * @throws UnsupportedOperationException if <code>!isInvertible()</code>.
     */
    public void invert( double x, double y, double[] out2x1, int outOff );

}
