using System;

namespace Numerics5
{
    internal class LeftBorderExpression : BorderExpression
    {
        public LeftBorderExpression(double alpha, double beta, Func<double, double> phi) : base(alpha, beta, phi) { }

        public Equation Approximate(double t, int n)
        {
            var coefs = ListUtils.New<double>(n + 1);
            coefs[0] = Beta;
            var result = Phi(t);
            return new Equation(coefs, result);
        }
    }
}
