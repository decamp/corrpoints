/*
 * Copyright (c) 2014. Philip DeCamp
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.corrpoints;

import java.util.Random;
import org.junit.*;

/**
 * @author decamp
 */
public class TransformFitterTest {
      

    /**
     * Precisely defined case, where there are exactly enough correspondence
     * points to define transform.  In this case, transform should provide
     * exact mapping between correspondence points.
     */
    @Test
    public void testProjectiveExact() {
        double[] xy = {0, 0.78,  1.6, 22,   4.5, 3.54,  3.23, 8.64 };  //Input points
        double[] uv = {0, 0.51, 2.2, 4.52, 5.11, 6.51,  6.44, 12.5 };  //Output points.
        
        //System.out.println(MatFormatter.formatToMatlab(2, 4, xy, 0, 2));
        //System.out.println(MatFormatter.formatToMatlab(2, 4, uv, 0, 2));
        
        int[] err = {0};
        ProjectiveTransform trans = TransformFitter.corrPointsToProjectiveTransform(xy, 0, uv, 0, xy.length / 2, err);
        
        if(trans == null) {
            Assert.assertTrue("Error code: " + err[0], false);
        }
        
        double[] arr = {0,0};
        
        for(int i = 0; i < xy.length; i += 2) {
            trans.apply(xy[i], xy[i+1], arr, 0);
            assertNear(uv[i  ], arr[0]);
            assertNear(uv[i+1], arr[1]);
            
            trans.invert(uv[i], uv[i+1], arr, 0);
            assertNear(xy[i  ], arr[0]);
            assertNear(xy[i+1], arr[1]);
        }
    }
    

    /**
     * Overdefined case, where there are more correspondence points than 
     * required to define parameters.  In this case, resulting transform
     * should find mapping between provided correspondence points that 
     * minimizes the mean squared error. 
     */
    @Test
    public void testProjectiveOverdefined() {
        
        final int n       = 8;
        
        final double[] xy = { 22.514499,  20.674984,  16.123578,  48.870477, -20.001173, -28.072239,  
                              23.120150,   4.651651, -29.437780, -20.782450,  26.536082,  6.718686, 
                             -44.802199,  10.357676,  43.936695,  38.989043 };
        
        final double[] uv = {  2.184435,  29.727904, -35.672099, -37.721562,  84.078352,  17.851000,  
                              36.839035,  21.534735, -14.164054,  30.117624, -18.619282, -33.589262,  
                              58.498715,  117.724447, -49.322516,  105.021209 };
        
        //Coeffs computed by Matlab for inverse transform using same correspondence points.
        final double[] matlabBackward= {   -1.776058168572404,
                                           -0.443822571897106,
                                           60.854375608255914,
                                           -1.867026090213108,
                                            0.708240985664850,
                                           36.750801105746348,
                                           -0.034209255179965,
                                            0.016499728627515,
                                            1.591756726315765 };
        
        //Coeffs computed by Matlab for transform using same set of correspondence points.
        final double[] matlabForward= {   -0.249684990254034,
                                          -0.819810079584899,
                                          28.473635833081730,
                                          -0.821772958067725,
                                           0.357185576033685,
                                          23.170390036717887,
                                           0.003152173427251,
                                          -0.021321447382477,
                                           1.000000000000000 };
      
        //System.out.println(MatFormatter.formatToMatlab(2, n, xy, 0, 2));
        //System.out.println(MatFormatter.formatToMatlab(2, n, uv, 0, 2));
        
        int[] err = {0};
        ProjectiveTransform trans = TransformFitter.corrPointsToProjectiveTransform(xy, 0, uv, 0, xy.length / 2, err);
        
        if(trans == null) {
            Assert.assertTrue("Error code: " + err[0], false);
        }
        
        assertNear(trans.coeffsRef(), 0, matlabForward, 0, 9);
        assertNear(trans.inverseCoeffrsRef(), 0, matlabBackward, 0, 9);
    }


    /**
     * Precisely defined case, where there are exactly enough correspondence
     * points to define transform.  In this case, transform should provide
     * exact mapping between correspondence points.
     */
    @Test
    public void testPoly3Exact() {
        Random rand = new Random(1);
        
        double[] xy = new double[20]; //Input points
        double[] uv = new double[20]; //Output points
        
        for(int i = 0; i < xy.length; i++) {
            xy[i] = rand.nextDouble() * (i + 1);
            uv[i] = xy[i] * (i + 1) + rand.nextDouble() * 0.5;
        }
        
        int m = xy.length / 2;
        
        //System.out.println(MatFormatter.formatToMatlab(1, m * 2, xy, 0, 1));
        //System.out.println(MatFormatter.formatToMatlab(1, m * 2, uv, 0, 1));
        int[] err = {0};
        
        PolyTransform trans = TransformFitter.corrPointsToPolyTransform(xy, 0, uv, 0, m, 3, err);
        
        if(trans == null) {
            Assert.assertTrue("Error code: " + err[0], false);
        }
        
        double[] arr = {0,0};
        
        for(int i = 0; i < xy.length; i += 2) {
            trans.apply(xy[i], xy[i+1], arr, 0);
            assertNear(uv[i  ], arr[0]);
            assertNear(uv[i+1], arr[1]);
        }
    }


    /**
     * Overdefined case, where there are more correspondence points than 
     * required to define parameters.  In this case, resulting transform
     * should find mapping between provided correspondence points that 
     * minimizes the mean squared error. 
     */
    @Test
    public void testPoly3Overdefined() {
        
        //Input points.
        double[] xy = {  0.725235,   0.272905,   0.873075,   1.640757,   1.159879,   0.582343,  
                         4.086223,   7.818314,   5.851982,   3.247697,   0.672492,   3.700837,  
                         1.628778,   9.079670,   7.443098,   3.921615,  10.826424,   4.125007,  
                        15.598678,  18.136991,   9.260569,   7.954107,   4.722733,   1.446441,  
                        10.098565,  17.413205,  21.466626,  12.035879,  23.931892,  25.639706,  
                        25.787209,   7.530278 };
        
        //Output points.
        double[] uv = {   0.770923,    0.853989,    2.818597,    6.616211,    6.136676,    3.554378,   29.079034,  
                         62.822321,   52.991624,   32.495535,    7.816865,   44.556080,   21.293793,  127.212082,  
                        112.087238,   62.902219,  184.191153,   74.401537,  296.753294,  363.133361,  194.661940,  
                        175.078684,  108.672078,   34.786610,  252.630061,  452.903119,  579.833255,  337.077597,  
                        694.337186,  769.553441,  799.447151,  241.427708 };
        
        //Results returned by matlab for same set of correspondence points.
        double[] matlab = { -10.703812978841524,
                             18.055914539503792,
                              2.196046209246088,
                              2.229817025911375,
                             -1.053956394097533,
                             -1.952577693092243,
                              0.154354809921372,
                             -0.329979732102789,
                              0.019207013965372,
                              0.201447835115973,
                             -8.870258409069743,
                             -1.894057161001981,
                             21.056191717907875,
                              4.199148300439201,
                             -0.902943088216324,
                             -4.016250431954839,
                              0.246009153422614,
                             -0.550126374596808,
                             -0.025304959803081,
                              0.370923854631815 };
                
//        for(int i = 0; i < xy.length; i++) {
//            xy[i] = rand.nextDouble() * (i + 1);
//            uv[i] = xy[i] * (i + 1) + rand.nextDouble() * 0.5;
//        }
//        
        int m = xy.length / 2;
        
        //System.out.println(MatFormatter.formatToMatlab(m*2, 1, xy, 0, 1));
        //System.out.println(MatFormatter.formatToMatlab(m*2, 1, uv, 0, 1));
        
        int[] err = {0};
        
        Poly3Transform trans = TransformFitter.corrPointsToPoly3Transform(xy, 0, uv, 0, m, err);
        
        if(trans == null) {
            Assert.assertTrue("Error code: " + err[0], false);
        }
        
        double[] arr = {0,0};
        
        //System.out.println(MatFormatter.formatToScreen(10, 2, trans.coeffsRef(), 0, 10));
        assertNear(trans.coeffsRef(), 0, matlab, 0, 20);
    }

    
    
    @Test
    public void testPoly2Exact() {
        Random rand  = new Random(3);
        double[][] p = generatePolyCorrPoints(rand, 6, 2);
        int[] err    = {0};
        Poly2Transform trans = TransformFitter.corrPointsToPoly2Transform(p[0], 0, p[1], 0, p[0].length / 2, err);
        
        if(trans == null) {
            Assert.assertTrue("Error code: " + err[0], false);
        }
        
        double[] arr = {0,0};
        
        for(int i = 0; i < p[0].length; i+= 2) {
            trans.apply(p[0][i], p[0][i+1], arr, 0);
            assertNear(p[1][i  ], arr[0]);
            assertNear(p[1][i+1], arr[1]);
        }
    }
    
    
    
    @Test
    public void testPolyNExact() {
        Random rand  = new Random(4);
        
        for(int degree = 1; degree < 11; degree++) {
            //System.out.println("Degree: " + degree);
            
            double[][] p = generatePolyCorrPoints(rand, -1, degree);
            int[] err    = {0};
        
            PolyTransform trans = TransformFitter.corrPointsToPolyTransform(p[0], 0, p[1], 0, p[0].length / 2, degree, err);
            
            if(trans == null) {
                Assert.assertTrue("Error code: " + err[0], false);
            }
            
            double[] arr = {0,0};

            for(int i = 0; i < p[0].length; i+= 2) {
                trans.apply(p[0][i], p[0][i+1], arr, 0);
                assertNear(p[1][i  ], arr[0]);
                assertNear(p[1][i+1], arr[1]);
            }
        }
    }
    
    

    public void testPoly3Symmetry() {
        Random rand = new Random(1);
        
        double[] xy = new double[30]; //Input points
        double[] uv = new double[30]; //Output points
        
        for(int i = 0; i < xy.length; i++) {
            xy[i] = rand.nextDouble() * (i + 1);
            uv[i] = 4 + xy[i] * 2 + xy[i] * xy[i] * 0.02; // + rand.nextDouble() * 0.02;
        }
        
        int m = xy.length / 2;
        
        //System.out.println(MatFormatter.formatToMatlab(1, m * 2, xy, 0, 1));
        //System.out.println(MatFormatter.formatToMatlab(1, m * 2, uv, 0, 1));
        int[] err = {0};
        
        Poly3Transform t1 = TransformFitter.corrPointsToPoly3Transform(xy, 0, uv, 0, m, err);
        Poly3Transform t2 = TransformFitter.corrPointsToPoly3Transform(uv, 0, xy, 0, m, err);
        
        double[] arr1 = {0,0};
        double[] arr2 = {0,0};
        
        for(int i = 0; i < xy.length; i += 2) {
            t1.apply(xy[i], xy[i+1], arr1, 0);
            t2.apply(arr1[0], arr1[1], arr2, 0);
            
            //double dx = Math.abs(xy[i  ] - arr2[0]);
            //double dy = Math.abs(xy[i+1] - arr2[1]);
            
            double dx = Math.abs( (xy[i  ] - arr2[0]) / xy[i  ] );
            double dy = Math.abs( (xy[i+1] - arr2[1]) / xy[i+1] );
            
            System.out.println( dx + "\t" + dy );
            //System.out.println(arr1[1] + "\t" + uv[i+1]);
            
        }
    }
    
    
    
    
    static double[][] generatePolyCorrPoints(Random rand, int n, int degree) {
        int k = PolyNTransform.coeffCount(degree);
        
        if(n < 0)
            n = k / 2;
        
        double[][] ret = new double[2][n*2];
        double[] c     = new double[PolyNTransform.coeffCount(degree)];  
        
        for(int i = 0; i < c.length; i++) {
            c[i] = rand.nextDouble() * 3.0 - 1.5;
        }
        
        for(int i = 0; i < n; i++) {
            ret[0][i*2  ] = rand.nextDouble() * 5.0 - 2.5;
            ret[0][i*2+1] = rand.nextDouble() * 5.0 - 2.5;
            
            PolyNTransform.apply(degree, c, ret[0][i*2], ret[0][i*2+1], ret[1], i*2);
        }
        
        return ret;
    }
    
    
    
    static void assertNear(double[] a, int aOff, double[] b, int bOff, int len) {
        for(int i = 0; i < len; i++) {
            double diff = Math.abs(a[aOff] - b[bOff]);
            Assert.assertTrue(diff < 1E-4);
        }
    }
    
    

    static void assertNear(double a, double b) {
        Assert.assertTrue(Math.abs(a - b) < 1E-5);
    }
    
}
