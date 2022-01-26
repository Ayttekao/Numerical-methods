using System;

namespace NumMeth4
{
    static class Program
    {
        static void Main(string[] args)
        {
            double eps = 1e-5;
            double res;
            double prev;
            long segCount = 2;
            double a = 0.0;
            double b = 1.0;
            
            do {
                prev = TrapezoidIntegral(a, b, segCount);
                segCount *= 2;
                res = TrapezoidIntegral(a, b, segCount);
		
            } while(Math.Abs(res - prev) > eps / 2);
            
            Console.WriteLine("Integral by the trapezoidal method: {0:f5}\n" +
                              "It took partitions: {1}\n" +
                              "Richardson extrapolation: {2:f5}", res, segCount, (4 * TrapezoidIntegral(a, b, segCount * 2) - res) / 3);
            
            segCount = 2;
            do {
                prev = SimpsonIntegral(a, b, segCount);
                segCount *= 2;
                res = SimpsonIntegral(a, b, segCount);
			
            } while(Math.Abs(res - prev) > eps);
            
            Console.WriteLine("Integral by Simpson's method: {0:f5}\n" +
                              "It took partitions: {1}\n" +
                              "Richardson Extrapolation: {2:f5}", res, segCount, (16 * SimpsonIntegral(a, b, segCount * 2) - res) / 15);
        }

        private static double Function(double x)
        {
            return Math.Sin(x) / (1 + Math.Pow(x, 2));
        }

        private static double TrapezoidIntegral(double left, double right, long segment)
        {
            double sum = 0;
            var step = (right - left) / (1.0 * segment);
            if (segment == 0) 
                return sum;
            
            for (uint count = 1; count < segment; ++count) 
            {
                sum += Function(left + count * step);
            }
	
            sum += (Function(left) + Function(right)) / 2;
            sum *= step;
            
            return sum;
        }
        
        static double SimpsonIntegral(double left, double right, long segments)

        {
            double integral = 0.0;
            double lenghtOfNewSegment = (left + right) / (2 * segments);
            for (long count = 0; count < segments; count++)
            {
                var tmp = left + count * 2 * lenghtOfNewSegment;
                integral += lenghtOfNewSegment / 3 * (Function(tmp) + 4 * Function(tmp + lenghtOfNewSegment) +
                                                      Function(tmp + 2 * lenghtOfNewSegment));
            }
            return integral;
        }
    }
}