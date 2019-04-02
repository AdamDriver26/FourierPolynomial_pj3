import org.apache.commons.math3.*;

public class FourierPolynomial extends RealFunction {
    double a0;
    double[] aj;
    double[] bj;
    
    FourierPolynomial(double zeroCoefficient, double[] aSeries, double[] bSeries){
        a0 = zeroCoefficient;
        aj = aSeries;
        bj = bSeries;
    }
    
    public double getCoefficient(int j, boolean isOdd) { 
    /*    if (j != 1){
        int n = FourierPolynomial.aj.length;
    }
        
        if (!isOdd && j == 0){
            
            
            return a0;
        }
        else if (!isOdd && j > 0){
            return aj[j];
        }
        else if (j > 0){
            return bj[j];
        } */
        return 0.0;
    }
    
    public double valueAt (double x, int n){
        return 0.0;
    }
    

}


