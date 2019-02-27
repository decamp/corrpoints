/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

/**
 * @author decamp
 */
public class MatFormatter {


    public static String formatToMatlab( int m, int n, double[] mat, int off, int stride ) {
        StringBuilder sb = new StringBuilder();
        formatToMatlab( m, n, mat, off, stride, sb );
        return sb.toString();
    }


    public static String formatToScreen( int m, int n, double[] mat, int off, int stride ) {
        StringBuilder sb = new StringBuilder();
        formatToScreen( m, n, mat, off, stride, sb, "" );
        return sb.toString();
    }


    public static void formatToMatlab( int m, int n, double[] mat, int off, int stride, StringBuilder s ) {
        s.append( "[" );
        for( int y = 0; y < m; y++ ) {
            for( int x = 0; x < n; x++ ) {
                if( x > 0 ) {
                    s.append( "  " );
                }
                s.append( String.format( "% 6.6f", mat[y + x * stride + off] ) );
            }
            s.append( "; " );
        }
        s.append( "]" );
    }


    public static void formatToScreen(
        int m,
        int n,
        double[] mat,
        int off,
        int stride,
        StringBuilder s,
        String prefix
    ) {
        s.append( prefix ).append( "Matrix (" ).append( m ).append( "x" ).append( n ).append( ")\n" );
        for( int y = 0; y < m; y++ ) {
            s.append( prefix );
            for( int x = 0; x < n; x++ ) {
                s.append( String.format( " % 8.4f", mat[y + x * stride + off] ) );
            }
            s.append( "\n" );
        }
    }

}
