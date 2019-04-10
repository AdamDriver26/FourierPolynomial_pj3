
/**
 * Write a description of class Test here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class Test extends RealFunction
{
    public double valueAt(double x)
    {
        return Math.exp(-x);
    }

    static FourierPolynomial testTransformer ()
    {
        Test f = new Test();
        
        return FourierTransformer.approximate(f,4);
    }
    
    public FourierPolynomial testHeat ()
    {
        Test f = new Test();
        HeatEquation h = new HeatEquation(0.1,f,5);
        
        FourierPolynomial ft = h.getSolution(12);
        return ft;
    }
}
