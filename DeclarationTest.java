import static org.junit.Assert.*;
import org.junit.Test;
import java.util.List;
import java.lang.reflect.*;

/**
 * Mathematical Skills II, Programming (2017/18);
 * Project 3.
 * 
 * Non-assessed unit tests.
 * 
 * These tests check whether some of the expected functions have the correct name and signature.
 */
public class DeclarationTest
{
    private Class<?> doubleArray;

    /**
     * Constructor for this test class.
     */
    public DeclarationTest()
    {
        doubleArray = (new double[0]).getClass();
    }

    /**
     * helper method to assert that a class is located in the default package
     */
    private void assertDefaultPackage(Class<?> targetClass)
    {
        String msg = targetClass.getName()+" must be in default package";
        assertNull(msg, targetClass.getPackage());
    }

    /**
     * helper method to verify the type parameter of a return value
     */
    private void assertReturnTypeParameter(Class<?> targetClass, String methodName, Class typePar, Class... paramType)
    {        
        String fullname = targetClass.getName()+"."+methodName; 
        String tpname = typePar.getName();
        String msg = "return type of "+fullname+" must have type parameter <"+ tpname +">; ";

        try
        {
            Method method = targetClass.getDeclaredMethod(methodName, paramType);
            Type returnType = method.getGenericReturnType();

            if(returnType instanceof ParameterizedType){
                ParameterizedType ptype = (ParameterizedType) returnType;
                Type[] typeArguments = ptype.getActualTypeArguments();
                assertEquals(msg, 1, typeArguments.length);
                assertEquals(msg, typePar, typeArguments[0]);
            }
            else 
            {
                fail(msg);
            }
        }
        catch (NoSuchMethodException e)
        {
            fail("Method "+fullname+" does not exist or has incorrect signature");
        }

    }

    /**
     * helper method for testing whether a (static) function is declared
     */
    private void testFunctionDeclared(Class<?> targetClass, String methodName, Class returnType, Class... paramType)
    {
        assertDefaultPackage(targetClass);

        String fullname = targetClass.getName()+"."+methodName; 
        try
        {
            Method method = targetClass.getDeclaredMethod(methodName, paramType);

            assertEquals("Return type of "+fullname, returnType, method.getReturnType());

            assertTrue(fullname+" not declared static", Modifier.isStatic(method.getModifiers()) );
            assertFalse(fullname+" must not be declared private", Modifier.isPrivate(method.getModifiers()) );
        }
        catch (NoSuchMethodException e)
        {
            fail("Function "+fullname+" does not exist or has incorrect signature");
        }                

    }

    /**
     * helper method for testing whether a method is declared
     */
    private void testMethodDeclared(Class<?> targetClass, String methodName, Class returnType, Class... paramType)
    {
        assertDefaultPackage(targetClass);

        String fullname = targetClass.getName()+"."+methodName; 
        try
        {
            Method method = targetClass.getDeclaredMethod(methodName, paramType);

            assertEquals("Return type of "+fullname, returnType, method.getReturnType());

            assertFalse(fullname+" must not be declared static", Modifier.isStatic(method.getModifiers()) );
            assertFalse(fullname+" must not be declared private", Modifier.isPrivate(method.getModifiers()) );
        }
        catch (NoSuchMethodException e)
        {
            fail("Method "+fullname+" does not exist or has incorrect signature");
        }                

    }

    /**
     * helper method for testing whether a constructor is declared
     */
    private void testConstructorDeclared(Class<?> targetClass, Class... paramType)
    {
        assertDefaultPackage(targetClass);

        String fullname = "Constructor in "+targetClass.getName(); 
        try
        {
            Constructor constructor = targetClass.getConstructor(paramType);

            assertFalse(fullname+" must not be declared private", Modifier.isPrivate(constructor.getModifiers()) );
        }
        catch (NoSuchMethodException e)
        {
            fail(fullname+" does not exist or has incorrect signature");
        }    
    }

    /**
     * helper method for testing whether *no* constructor is declared 
     * (for composite data types)
     */
    private void testNoConstructorDeclared(Class<?> targetClass)
    {
        assertDefaultPackage(targetClass);

        boolean hasConstructors = targetClass.getConstructors().length > 0;

        String msg = targetClass.getName() + " must not have constructors";

        assertFalse(msg, hasConstructors);
    }

    /**
     * helper method for testing whether a field is declared
     */
    private void testFieldDeclared(Class<?> targetClass, String fieldName, Class type)
    {
        assertDefaultPackage(targetClass);

        String fullname = targetClass.getName()+"."+fieldName; 
        try
        {
            Field field = targetClass.getDeclaredField(fieldName);

            assertEquals("Type of field "+fullname, type, field.getType());

            assertFalse(fullname+" must not be declared static", Modifier.isStatic(field.getModifiers()) );
            assertFalse(fullname+" must not be declared private", Modifier.isPrivate(field.getModifiers()) );
        }
        catch (NoSuchFieldException e)
        {
            fail("Field "+fullname+" does not exist");
        }      
    }

    /**
     * helper method for testing whether a class is a subclass of another one
     */
    private void testSubclassOf(Class<?> targetClass, Class<?> superClass)
    {
        assertDefaultPackage(targetClass);

        String msg = "Checking whether "+targetClass.getName()+ " is a subclass of "+superClass.getName();
        assertTrue(msg, superClass.isAssignableFrom(targetClass));
    }

    /**
     * Tests whether the constructor of FourierPolynomial are declared correctly.
     */
    @Test
    public void fourierPolyConstructorDeclaredTest()
    {
        testConstructorDeclared(FourierPolynomial.class, double.class, doubleArray, doubleArray);
    }

    /**
     * Tests whether the getCoefficient method of a FourierPolynomial are declared.
     */
    @Test
    public void getCoefficientDeclaredTest()
    {
        testMethodDeclared(FourierPolynomial.class, "getCoefficient", double.class, int.class, boolean.class);
    }

    /**
     * Tests whether the valueAt and derviativeValueAt methods are declared.
     */
    @Test
    public void valueAtDeclaredTest()
    {
        testMethodDeclared(FourierPolynomial.class, "valueAt", double.class, double.class);
        testMethodDeclared(FourierPolynomial.class, "derivativeValueAt", double.class, double.class);
    }

    /**
     * Tests whether the method add() is declared.
     */
    @Test
    public void addDeclaredTest()
    {
        testMethodDeclared(FourierPolynomial.class, "add", FourierPolynomial.class, FourierPolynomial.class);
    }

    /**
     * Tests whether the method multiply() is declared.
     */
    @Test
    public void multiplyDeclaredTest()
    {
        testMethodDeclared(FourierPolynomial.class, "multiply", FourierPolynomial.class, FourierPolynomial.class);
    }

    /**
     * Tests whether the methods derivative() and antiderivative() are declared.
     */
    @Test
    public void derivativeDeclaredTest()
    {
        testMethodDeclared(FourierPolynomial.class, "derivative", FourierPolynomial.class);
        testMethodDeclared(FourierPolynomial.class, "antiderivative", FourierPolynomial.class);
    }

    /**
     * Tests whether the class FourierTransformer and its methods are declared.
     */
    @Test
    public void transformerDeclaredTest()
    {
        testFunctionDeclared(FourierTransformer.class, "approximate", FourierPolynomial.class, RealFunction.class, int.class);
    }

    /**
     * Tests whether the class HeatEquation and its methods are declared.
     */
    @Test
    public void heatEquationDeclaredTest()
    {
        testConstructorDeclared(HeatEquation.class, double.class, RealFunction.class, int.class);        
        testMethodDeclared(HeatEquation.class, "evaluateSolution", double.class, double.class, double.class);
        testMethodDeclared(HeatEquation.class, "getSolution", FourierPolynomial.class, double.class);
    }

}
