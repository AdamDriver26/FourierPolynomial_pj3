/**
 * An extension to the RealFunction class which represents Fourier polynomials and allows
 * for their addition and multiplication as well as finding derivatives and antiderivatives. 
 */
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
     * Constructor for the FourierPolynomial class. 
     */
    public FourierPolynomial(double zeroCoefficient, double[] aArray, double[] bArray)
    {
        a0 = zeroCoefficient;
        // Uses the lengthEqualiser method to ensure aj and bj are the same length.
        aj = lengthEqualiser(aArray, bArray.length);
        bj = lengthEqualiser(bArray, aArray.length);
    }
    
    /**
     * Extends an array with zeros if it is shorter than the given value.
     * 
     * @param shortArray the array to be extended.
     * @param n the length the array should be extended to.
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
     * Returns the the jth coefficient of the Fourier polynomial.
     * 
     * @param j the term of the coefficient.
     * @param wantB determines whether to return the coefficient of the cosine or sine value.
     * @return the jth coefficient.
     * @throws IllegalArgumentException if j is negative.
     */
    public double getCoefficient(int j, boolean wantB)
    {
        if (j < 0)
        {
            throw new java.lang.IllegalArgumentException("j must be non-negative");
        }

        int n = aj.length; // Could use either aj or bj since they are the same length by construction. This convention is used throughout.
        if (j > n)
        {
            return 0.0;
        }
        else 
        {
            if (j==0)
            {
                return a0;
            }
            else if (!wantB)
            {
                return aj[j-1];
            }
            else
            {
                return bj[j-1];
            }        
        }
    }
    
    /**
     * Calculates the value of the Fourier polynomial at a point.
     * 
     * @param x the point to consider.
     * @return the value of the Fourier polynomial at the point x.
     */
    public double valueAt (double x)
    {
        int n = aj.length;
        double value = a0/2.0; // Immediatly calculates the constant value which is unaffected by x.
        for (int j=1; j<=n; j++)
        {
            value += aj[j-1]*Math.cos(j*x) + bj[j-1]*Math.sin(j*x);
        }
        return value;
    }
    
    /**
     * Calculates the value of the Fourier polynomials derivative at a point.
     * 
     * @param x the point to consider.
     * @return the value of the Fourier polynomials derivative at the point x.
     */
    public double derivativeValueAt (double x)
    {
        int n = aj.length;
        double value = 0.0;
        for (int j=1; j<=n; j++)
        {
            value += j*(bj[j-1]*Math.cos(j*x) - aj[j-1]*Math.sin(j*x));
        }
        return value;
    }
    
    /**
     * Adds a Fourier polynomial.
     * 
     * @param f the Fourier polynomial to add.
     * @return the coefficients of the sum as a FourierPolynomial object.
     * @see FourierPolynomial
     */
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
    
    /**
     * Multiplies the Fourier polynomial with another given Fourier polynomial.
     * 
     * @param f the FourierPolynimal to multiply by.
     * @return the coefficients of the product as a FourierPolynomial object.
     */
    public FourierPolynomial multiply (FourierPolynomial f)
    {
        int l = this.aj.length; // order of the original FourierPolynomial
        int m = f.aj.length; // order of the given FourierPolynomial
        int n = l+m;
        
        double newA0 = (f.a0*this.a0)/2.0; // True for all Fourier polynomials
        // Includes constant values derived from other coefficients 
        for (int i=0; i<l; i++)
        {
            newA0 = this.aj[i]*f.aj[i] + this.bj[i]*f.bj[i]; 
        }
        

        double[] newA = new double[n];
        double[] newB = new double[n];
        
        for (int i=1; i<=n; i++)
        {
            if (i<l)
            {
                newA[i] += this.a0*f.aj[i-1]/2.0; // Adds values of the original a0 times new aj
                newB[i] += this.a0*f.bj[i-1]/2.0; // Adds values of the original a0 times new bj
            }
            if (i<m)
            {
                newA[i] += f.a0*this.aj[i-1]/2.0; // Adds values of the new a0 times original aj
                newB[i] += f.a0*this.bj[i-1]/2.0; // Adds values of the new a0 times original bj
            }
        }
        
        for (int i=1; i<=n; i++) // Cycles through newFourier coefficients to be found
        {
            for (int p=1; p<l; p++) // Cycles through indexes of the original FourierPolynomial used in calculation
            {
                for (int q=1; q<m; q++) // Cycles through indexes of the FourierPolynomial f used in calculation
                {
                    if (p+q == i)
                    {
                        newA[i-1] += (this.aj[p-1]*f.aj[q-1] - this.bj[p-1]*f.bj[q-1])/2.0;
                        newB[i-1] += (this.aj[p-1]*f.bj[q-1] + this.bj[p-1]*f.aj[q-1])/2.0;
                    }
                    
                    if (i<l && i<m)
                    {
                        newA[i-1] += (this.aj[p-1]*f.aj[p+i-1] + this.bj[p-1]*f.bj[p+i-1])/2.0;
                        newB[i-1] += (this.aj[p-1]*f.bj[p+i-1] - this.bj[p-1]*f.aj[p+i-1])/2.0;
                        if (p-i>0)
                        {
                            newA[i-1] += (this.aj[p-1]*f.aj[p-i-1] + this.bj[p-1]*f.bj[p-i-1])/2.0;
                            newB[i-1] += (-this.aj[p-1]*f.bj[p+i-1] + this.bj[p-1]*f.aj[p+i-1])/2.0;
                        }
                    }
                }
            }
        }
        
        FourierPolynomial newFourier = new FourierPolynomial(newA0,newA,newB);
        return newFourier;
    }
    
    private static boolean insignificant (double x)
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
    
    /**
     * Finds the derivative of the Fourier polynomial.
     * 
     * @return the coefficients of the derivative as a FourierPolynomial object.
     * @see FourierPolynomial
     */
    public FourierPolynomial derivative ()
    {
        int n = this.aj.length;
        double[] dAj = new double[n];
        double[] dBj = new double[n];
        
        for (int i=1; i<=n; i++)
        {
            dAj[i-1] = i*this.bj[i-1];
            dBj[i-1] = -i*this.aj[i-1];
        }
        return new FourierPolynomial(0,dAj,dBj);
    }
    
    /**
     * Finds the antiderivative of the Fourier polynomial.
     * 
     * @return the coefficients of the antiderivative as a FourierPolynomial object.
     * @throws IllegalArgumentException if |a0| > 10^(-10). 
     * @see FourierPolynomial
     */
    public FourierPolynomial antiderivative ()
    {
        if (!insignificant(this.a0))
        {
            throw new java.lang.IllegalArgumentException("The antiderivative cannot be found for non-zero a0");
        }
        
        int n = this.aj.length;
        double[] aAj = new double[n];
        double[] aBj = new double[n];
        
        for (int i=1; i<=n; i++)
        {
            aAj[i-1] = -this.bj[i-1]/i;
            aBj[i-1] = this.aj[i-1]/i;
        }
        return new FourierPolynomial(0,aAj,aBj);
    }
    
    
    // REMOVE THIS BEFORE SUBMITTING
    public static FourierPolynomial test ()
    {
        double Xa0 = 0;
        double[] Xaj = {0};
        double[] Xbj = {0};
        FourierPolynomial X = new FourierPolynomial(Xa0,Xaj,Xbj);
        
        double t1a0 = -1;
        double[] t1aj = {-2,-3};
        double[] t1bj = {-4,-5};
        FourierPolynomial test1 = new FourierPolynomial(t1a0,t1aj,t1bj);
        
        return X.antiderivative();
    }
}