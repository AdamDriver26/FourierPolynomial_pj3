/**
 * Used to estimate the solution of the heat equation of a given RealFunction object at a point or as a FourierPolynomial object.
 */
public class HeatEquation
{
    /**
     * Diffusivity of the medium.
     */
    final double alpha;
    
    /**
     * Initial conditions as a 2π-periodic RealFunction object.
     */
    final RealFunction g;
    
    /**
     * Degree of the resulting FourierPolynomial.
     */
    final int n;
    
    /**
     * Constructor for the HeatEquation class. 
     * 
     * @throws IllegalArgumentException if diffusivity is not greater than zero.
     * @throws IllegalArgumentException if degree is negative.
     */
    public HeatEquation (double diffusivity, RealFunction initialCondition, int degree)
    {
        if (diffusivity <= 0)
        {
            throw new java.lang.IllegalArgumentException("alpha must be greater than zero");
        }
        alpha = diffusivity;
        
        g = initialCondition;
        
        if (degree < 0)
        {
            throw new java.lang.IllegalArgumentException("n cannot be negative");
        }
        n = degree;
    }
    
    /**
     * Approximates the function g with a Fourier polynomial using the approximate method in the FourierTranformer class.
     * Uses the solution to the heat equation given in equation (4) to give the value of u at a point x at time t.
     * 
     * @param x the point in space.
     * @param t time passed.
     * @return the value of x at time t.
     * @throws IllegalArgumentException if t is negative.
     * @see FourierPolynomial
     * @see approximate
     */
    public double evaluateSolution (double x, double t)
    {
        if (t < 0)
        {
            throw new java.lang.IllegalArgumentException("t must be non-negative");
        }
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
     * @param t time passed.
     * @return a FourierPolynomial representing the change after time t.
     * @throws IllegalArgumentException if t is negative.
     */
    public FourierPolynomial getSolution (double t)
    {
        if (t < 0)
        {
            throw new java.lang.IllegalArgumentException("t must be non-negative");
        }
        FourierPolynomial f = FourierTransformer.approximate(g, n);
        
        for (int j=1; j<n; j++)
        {
            f.aj[j-1] *= Math.exp(-alpha*j*j*t);
            f.bj[j-1] *= Math.exp(-alpha*j*j*t);
        }
        return f;
    }
}
