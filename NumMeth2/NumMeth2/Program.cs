using System;

namespace NumMeth2
{
    internal static class Program
    {
        private static readonly double Eps = 1e-5;

        private static double FirstLineFunction(double x, double y)
        {
            return Math.Pow(x, 2) * Math.Pow(y, 2) - 2 * Math.Pow(x, 3) - 5 * Math.Pow(y, 3) + 10;
        }

        private static double SecondLineFunction(double x, double y)
        {
            return Math.Pow(x, 4) - 8 * y + 1;
        }

        private static double XDerivativeFunctionFirstLine(double x, double y)
        {
            return 2 * x * Math.Pow(y, 2) - 6 * Math.Pow(x, 2);
        }

        private static double YDerivativeFunctionFirstLine(double x, double y)
        {
            return 2 * Math.Pow(x, 2) * y - 15 * Math.Pow(y, 2);
        }

        private static double XDerivativeFunctionSecondLine(double x)
        {
            return 4 * Math.Pow(x, 3);
        }
        
        private static double YDerivativeFunctionSecondLine()
        {
            return -8;
        }

        private static double XDerivativeNumFormFunctionFirstLine(double x, double y)
        {
            return (1 / Eps) * (FirstLineFunction(x + Eps, y) - FirstLineFunction(x, y));
        }
        private static double YDerivativeNumFormFunctionFirstLine(double x, double y)
        {
            return (1 / Eps) * (FirstLineFunction(x, y + Eps) - FirstLineFunction(x, y));
        }
        private static double XDerivativeNumFormFunctionSecondLine(double x, double y)
        {
            return (1 / Eps) * (SecondLineFunction(x + Eps, y) - SecondLineFunction(x, y)) ;
        }
        private static double YDerivativeNumFormFunctionSecondLine(double x, double y)
        {
            return (1 / Eps) * (SecondLineFunction(x, y + Eps) - SecondLineFunction(x, y)) ;
        }

        private static void InvertMatrix(double[ , ] a)
        {
            double det = a[0, 0] * a[1, 1] - a [0, 1] * a [1, 0];
            double temp = a[0, 0];
            a[0, 0] = a [1, 1] / det;
            a[1, 1] = temp / det;
            a[1, 0] = -a[ 1, 0 ] / det;
            a[ 0, 1 ] = -a [ 0, 1 ] / det;
        }

        private static void FuncAnalyticNumeric(double x, double y, bool choice, ref double xRes, ref double yRes,
            ref int iteration)
        {
            double[,] a = new double[2, 2];
            double norm;
            double[] b = new double[2];
            xRes = 0; yRes = 0; iteration = 0;
            do
            {
                if (choice)
                {
                    a[0, 0] = XDerivativeFunctionFirstLine(x, y);
                    a[0, 1] = YDerivativeFunctionFirstLine(x, y);
                    a[1, 0] = XDerivativeFunctionSecondLine(x);
                    a[1, 1] = YDerivativeFunctionSecondLine();
                }
                else
                {
                    a[0, 0] = XDerivativeNumFormFunctionFirstLine(x, y);
                    a[0, 1] = YDerivativeNumFormFunctionFirstLine(x, y);
                    a[1, 0] = XDerivativeNumFormFunctionSecondLine(x, y);
                    a[1, 1] = YDerivativeNumFormFunctionSecondLine(x, y);
                }
                InvertMatrix(a);
                double dx = -a[0, 0] * FirstLineFunction(x, y) + -a[0, 1] * SecondLineFunction(x, y);
                double dy = -a[1, 0] * FirstLineFunction(x, y) + -a[1, 1] * SecondLineFunction(x, y);
                x += dx;
                y += dy;
                b[0] = FirstLineFunction(x, y);
                b[1] = SecondLineFunction(x, y);
                norm = Math.Sqrt(b[0] * b[0] + b[1] * b[1]);
                iteration++;
            } while (norm >= Eps);
            xRes = x;
            yRes = y;
        }

        private static void Main()
        {
            double xRes = 0;
            double yRes = 0;
            int iteration = 0;
            
            FuncAnalyticNumeric(-3, 3, true, ref xRes, ref yRes, ref iteration);
            Console.WriteLine("Analytical Derivative");
            Console.WriteLine("Number of iterations: {0}\n x = {1:F6}, y = {2:F6} ", iteration,
                xRes, yRes);
            FuncAnalyticNumeric(1, 1, false, ref xRes, ref yRes, ref iteration);
            Console.WriteLine("Number of iterations: {0}\n x = {1:F6}, y = {2:F6} ", iteration,
                xRes, yRes);
            
            FuncAnalyticNumeric(-3, 3, true, ref xRes, ref yRes, ref iteration);
            Console.WriteLine("\nNumerical Derivative");
            Console.WriteLine("Number of iterations: {0}\n x = {1:F6}, y = {2:F6} ", iteration,
                xRes, yRes);
            FuncAnalyticNumeric(1, 1, false, ref xRes, ref yRes, ref iteration);
            Console.WriteLine("Number of iterations: {0}\n x = {1:F6}, y = {2:F6} ", iteration,
                xRes, yRes);
        }
    }
}