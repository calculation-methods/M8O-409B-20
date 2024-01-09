using System;
using System.Reflection;

namespace Numerics7
{
    internal class TopBorderExpression : BorderExpression
    {
        public TopBorderExpression(double alpha, double beta, Func<double, double> fi) : base(alpha, beta, fi) { }

        public override Equation Approximate1D(IssueParameters parameters, int xIndex)
        {
            var coefs = ListUtils.New<double>(parameters.Ny + 1);
            var x = xIndex * parameters.Hx;

            coefs[0] = -Alpha / parameters.Hy + Beta;
            coefs[1] = Alpha / parameters.Hy;

            return new Equation(coefs, Fi(x));
        }
    }
}
