using System;

namespace NumMeth1
{
    internal static class Program
    {
        private static double Function(double x)
        {
            return Math.Exp(x) - 6 * x - 3;
        }

        private static double DerivativeOfAFunction(double x)
        {
            return Math.Exp(x) - 6;
        }

        private static Tuple<double, int> Dichotomy(double lowerBound, double upperBound)
        {
            const double eps = 1e-5;
            int count = 0;
            while (upperBound - lowerBound > eps)
            {
                double tmpValue = (lowerBound + upperBound) / 2;
                if (Function(upperBound) * Function(tmpValue) < 0)
                    lowerBound = tmpValue;
                else
                    upperBound = tmpValue;
                count++;
            }
            
            return new Tuple<double, int>((lowerBound + upperBound) / 2, count);
        }

        private static Tuple<double, int> SimpleIteration(double lowerBound, double upperBound)
        {
            int count = 0;
            const double eps = 1e-5;
            int sign = Function(lowerBound) < 0 ? -1 : 1;
            double tau = - 2 / Math.Abs(DerivativeOfAFunction(lowerBound) + DerivativeOfAFunction(upperBound));
            double x = upperBound;
            
            while (Math.Abs(Function(x)) > eps)
            {
                count++;
                x += sign * tau * Function(x);
            }

            return new Tuple<double, int>(x, count);
        }

        static void Main()
        {
            Tuple<double, int> firstDichotomy = Dichotomy(-3, 0);
            Tuple<double, int> secondDichotomy = Dichotomy(0, 4);
            Tuple<double, int> firstIteration = SimpleIteration(-3, 0);
            Tuple<double, int> secondIteration = SimpleIteration(0,3);
            
            Console.WriteLine("Dichotomy method");
            Console.WriteLine("Value = {0:F6}, Number of iterations {1} (-3 - lower bound, 0 - upper bound)",
                firstDichotomy.Item1, firstDichotomy.Item2);
            Console.WriteLine("Value = {0:F6}, Number of iterations {1} (0 - lower bound, 4 - upper bound)",
                secondDichotomy.Item1, secondDichotomy.Item2);
            
            Console.WriteLine("Simple iteration method");
            Console.WriteLine("Value = {0:F6}, Number of iterations {1} (-3 - lower bound, 0 - upper bound)",
                firstIteration.Item1, firstIteration.Item2);
            Console.WriteLine("Value = {0:F6}, Number of iterations {1} (0 - lower bound, 3 - upper bound)", 
                secondIteration.Item1, secondIteration.Item2);
        }
    }
}
