using System;

namespace Numerics5
{
    internal abstract class BorderExpression
    {
        public BorderExpression(double alpha, double beta, Func<double, double> fi)
        {
            Alpha = alpha;
            Beta = beta;
            Phi = fi;
        }

        public double Alpha { get; }
        public double Beta { get; }
        public Func<double, double> Phi { get; }
    }
}
