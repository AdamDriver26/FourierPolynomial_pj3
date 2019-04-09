
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
        return Math.pow(1-x*x,0.5)+Math.pow(x,0.6);
    }

    static FourierPolynomial testTransformer ()
    {
        Test f = new Test();
        
        return FourierTransformer.approximate(f,10);
    }
    
    static double result()
    {
        Test f = new Test();
        return FourierTransformer.compositeTrapezium(f,0.0,2*Math.PI,10000)/Math.PI;
    }
}
