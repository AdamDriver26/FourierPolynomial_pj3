import org.apache.commons.math3.*;

public class FourierPolynomial extends RealFunction 
{
    /**
     * The zeroth Fourier coefficient.
     */
    double a0;
    
    /**
     * An array of the cosine Fourier coefficients.
     */
    double[] aj;
    
    /**
     * An array of the sine Fourier coefficients.
     */
    double[] bj;
    
    /**
     * Constructor for the class. 
     * 
     * @param zeroCoefficient value of the constant a0.
     * @param aArray value of the aj array.
     * @param bArray value of the bj array.
     */
    FourierPolynomial(double zeroCoefficient, double[] aArray, double[] bArray)
    {
        a0 = zeroCoefficient;
        // Uses the lengthEqualiser method to ensure aj and bj are the same length.
        aj = lengthEqualiser(aArray, bArray.length);
        bj = lengthEqualiser(bArray, aArray.length);
    }
    
    /**
     * Extends an array with zeros if it is shorter than the given value.
     * 
     * @param the array to be extended.
     * @param the length the array should be extended to.
     * @return an array extended with zeros to length n. 
     */
    private static double[] lengthEqualiser(double[] shortArray, int n)
    {
        // Simply returns the original array if it is not shorter than n.
        if (shortArray.length >= n)
        {
            return shortArray;
        }
        // Otherwise creates a new array of length n which includes the values from the original.
        else
        {
            double[] extendedArray = new double[n];
            for (int i=0; i<shortArray.length; i++)
            {
                extendedArray[i] = shortArray[i];
            }
            return extendedArray;
        }
    }
    
    /**
     * 
     */
    public double getCoefficient(int j, boolean wantB)
    {
        int n = aj.length; //Could use either aj or bj since they are the same length when constructed.
        if (j >= n)
        {
            return 0.0;
        }
        else 
        {
            if (!wantB && j==0)
            {
                return a0;
            }
            else if (!wantB)
            {
                return aj[j+1];
            }
            else
            {
                return bj[j+1];
            }
        }
    }
    
    public double valueAt (double x)
    {
        int n = aj.length;
        double value = a0/2.0;
        for (int j=1; j<n; j++)
        {
            value += aj[j]*Math.cos(j*x) + bj[j]*Math.sin(j*x);
        }
        return value;
    }
    
    public double derivativeValueAt (double x)
    {
        int n = aj.length;
        double value = 0.0;
        for (int j=1; j<n; j++)
        {
            value += bj[j]*Math.cos(j*x) - aj[j]*Math.sin(j*x);
        }
        return value;
    }
    
    public FourierPolynomial add (FourierPolynomial f)
    {
        int n = this.aj.length;
        f.aj = lengthEqualiser(f.aj, n);
        f.bj = lengthEqualiser(f.bj, n);
        
        f.a0 += this.a0;
        for (int i=0; i<n; i++)
        {
            f.aj[i] += this.aj[i];
            f.bj[i] += this.bj[i];
        }
        return f;
    }
    
    public FourierPolynomial multiply (FourierPolynomial f)
    {
        int k = f.aj.length; // order of the given FourierPolynomial
        int l = this.aj.length; // order of the original FourierPolynomial
        int n = k+l;
        
        double newA0 = f.a0; // Wrong
        double[] newA = new double[n];
        double[] newB = new double[n];
        //...
        FourierPolynomial newFourier = new FourierPolynomial(newA0,newA,newB);
        return newFourier;
    }
    
    private static boolean checkInRange (double x)
    {
        if (Math.abs(x) <= Math.pow(10,-10))
        {
            return true;
        }
        else 
        {
            return false;
        }
    }
    
    public FourierPolynomial derivative (FourierPolynomial f)
    {
        //...
        return f;
    }
    
    public FourierPolynomial antiderivative (FourierPolynomial f)
    {
        //...
        return f;
    }
    
    
    // REMOVE THIS BEFORE SUBMITTING
    public static FourierPolynomial test ()
    {
        double Xa0 = 2.0;
        double[] Xaj = {2.5,3.0,3.5};
        double[] Xbj = {4.0,4.5};
        FourierPolynomial X = new FourierPolynomial(Xa0,Xaj,Xbj);
        
        double t1a0 = 1.0;
        double[] t1aj = {0,0,0,1};
        double[] t1bj = {0,0,0,0};
        FourierPolynomial test1 = new FourierPolynomial(t1a0,t1aj,t1bj);
        
        return X.add(test1);
    }
}


