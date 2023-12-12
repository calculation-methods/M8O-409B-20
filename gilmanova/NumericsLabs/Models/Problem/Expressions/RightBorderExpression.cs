using System;

namespace Numerics5
{
    internal class RightBorderExpression : BorderExpression
    {
        public RightBorderExpression(double alpha, double beta, Func<double, double> phi) : base(alpha, beta, phi) { }

        public Equation Approximate(double t, int n)
        {
            var coefs = ListUtils.New<double>(n + 1);
            coefs[n] = Beta;
            var result = Phi(t);
            return new Equation(coefs, result);
        }
    }
}
