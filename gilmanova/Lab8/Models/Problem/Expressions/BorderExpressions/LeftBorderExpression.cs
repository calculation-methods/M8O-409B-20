using System;

namespace Numerics8
{
    internal class LeftBorderExpression : BorderExpression
    {
        public LeftBorderExpression(double alpha, double beta, Func<double, double, double> fi) : base(alpha, beta, fi) { }


        public override Equation Approximate1D(IssueParameters parameters, int yIndex, double t)
        {
            var coefs = ListUtils.New<double>(parameters.Nx + 1);
            var y = yIndex * parameters.Hy;

            coefs[0] = -Alpha / parameters.Hx + Beta;
            coefs[1] =  Alpha / parameters.Hx;

            return new Equation(coefs, Fi(y, t));
        }
    }
}
