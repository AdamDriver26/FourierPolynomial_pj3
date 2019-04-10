
public class HeatEquation
{
    final double alpha;
    final RealFunction g;
    final int n;
    
    /**
     * Constructor for the HeatEquation class. 
     * 
     * @param alpha the diffusivity of a material/space.
     * @param g the initial conditions as a 2π-periodic function.
     * @param n the degree of the associated Fourier polynomial.
     */
    HeatEquation (final double diffusivity, final RealFunction initialCondition, final int degree)
    {
        alpha = diffusivity;
        g = initialCondition;
        n = degree;
    }
    
    /**
     * Approximates the function g with a Fourier polynomial using the approximate method in the FourierTranformer class.
     * Uses the solution to the heat equation given in equation (4) to give the value of u at a point x at time t.
     * 
     * @param x the point in space.
     * @param t time passed.
     * @return the temperature of x at time t.
     * @see FourierPolynomial
     * @see approximate
     */
    public double evaluateSolution (double x, double t)
    {
        FourierPolynomial f = FourierTransformer.approximate(g, n);
        double u = f.a0/2.0;
        
        for (int j=1; j<n; j++)
        {
            u += f.aj[j-1]*Math.exp(-alpha*j*j*t) + f.bj[j-1]*Math.exp(-alpha*j*j*t); // From Equation 4 on the project description.
        }
        return u;
    }
    
    /**
     * Approximates the function g with a Fourier polynomial using the approximate method in the FourierTranformer class.
     * Uses the solution to the heat equation given in equation (4) to represent the temperature profile after time t with a FourierPolynomial object.
     * 
     * @param t  time passed.
     * @return a FourierPolynomial representing the temperature profile.
     */
    public FourierPolynomial getSolution (double t)
    {
        FourierPolynomial f = FourierTransformer.approximate(g, n);
        
        for (int j=1; j<n; j++)
        {
            f.aj[j-1] *= Math.exp(-alpha*j*j*t);
            f.bj[j-1] *= Math.exp(-alpha*j*j*t);
        }
        return f;
    }
}
