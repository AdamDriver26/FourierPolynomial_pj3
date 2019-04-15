/**
 * Used to transform RealFunction objects in to FourierPolynomial objects.
 */
public class FourierTransformer
{
    /**
     * Approximates a real function on the interval [0,2π] using the composite trapezium rule with 10,000 subintervals.
     * 
     * @param g the RealFunction object to be approximated.
     * @param n the degree of the resulting Fourier polynomial, assumed to be non-negative.
     * @return a FourierPolynomial object which approximates g.
     * @throws IllegalArgumentException if n is negative.
     * @see FourierPolynomial
     */
    public static FourierPolynomial approximate (RealFunction g, int n)
    {
        if (n < 0)
        {
            throw new java.lang.IllegalArgumentException("n cannot be negative");
        }
        int subInt = 10000;
        double h = (2.0*Math.PI)/ (double) subInt;
        double a0 = h*( g.valueAt(0.0) + g.valueAt(2.0*Math.PI) )/(2.0*Math.PI); // Start and end points.
        double[] aj = new double[n];
        double[] bj = new double[n];

        for (int i=1; i<subInt; i++)
        {
                a0 += h*( g.valueAt(i*h) )/Math.PI;
        }
        if (n == 0) // If n=0 returns a0 with aj,bj empty.
        {
            return new FourierPolynomial(a0,aj,bj);
        }
        else // Otherwise calculates the values of aj,bj up to degree n.
        {
            for (int j=1; j<=n; j++)
            {
                // Adds the g(a) and g(b) terms for each value of j.
                aj[j-1] += h*( g.valueAt(0.0)*Math.cos(0.0) + g.valueAt(2.0*Math.PI)*Math.cos(j*2.0*Math.PI) )/(2.0*Math.PI); // Start and end points.
                bj[j-1] += h*( g.valueAt(0.0)*Math.sin(0.0) + g.valueAt(2.0*Math.PI)*Math.sin(j*2.0*Math.PI) )/(2.0*Math.PI); // Start and end points.
                
                for (int i=1; i<subInt; i++)
                {
                    aj[j-1] += h*( g.valueAt(i*h)*Math.cos(j*i*h) )/Math.PI; // h*g(x)*cos(j*x)/PI value of the cosine integral at x=i*h
                    bj[j-1] += h*( g.valueAt(i*h)*Math.sin(j*i*h) )/Math.PI; // h*g(x)*sin(j*x)/PI value of the sine integral at x=i*h
                }
    
            }
            
            return new FourierPolynomial(a0,aj,bj);
        }
        
    } 
}
