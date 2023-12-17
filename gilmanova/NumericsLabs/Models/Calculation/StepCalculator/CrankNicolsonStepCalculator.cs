using System.Collections.Generic;
using System.Linq;

using Solution = System.Collections.Generic.List<double>;
using ReadOnlySolution = System.Collections.Generic.List<double>;

namespace Numerics5
{
    internal class CrankNicolsonStepCalculator : IStepCalculator
    {
        private readonly Issue issue;

        public CrankNicolsonStepCalculator(Issue issue)
        {
            this.issue = issue;
        }

        public List<Solution> Calculate()
        {
            var result = ListUtils.New<Solution>(issue.Parameters.K + 1);

            result[0] = ListUtils.New<double>(issue.Parameters.N + 1);
            for (int j = 0; j <= issue.Parameters.N; j++)
            {
                var x = issue.Parameters.LeftX + issue.Parameters.H * j;
                result[0][j] = issue.StartExpression.f(x);
            }

            for (int k = 1; k <= issue.Parameters.K; k++)
            {
                var t = issue.Parameters.Tau * k;
                result[k] = CalculateNextStep(result[k - 1], t);
            }

            return result;
        }

        private Solution CalculateNextStep(ReadOnlySolution previousStep, double t)
        {
            var equations = new List<Equation>();
            equations.Add(issue.LeftBorderExpression.Approximate(t, issue.Parameters.N));
            equations.AddRange(Enumerable.Range(1, issue.Parameters.N - 1).Select(j => GetJthEquation(previousStep, t, j)));
            equations.Add(issue.RightBorderExpression.Approximate(t, issue.Parameters.N));

            var system = new LinearSystem(equations);
            return system.Solve3Dim();
        }

        private Equation GetJthEquation(ReadOnlySolution u, double t, int j)
        {
            var coefs = ListUtils.New<double>(u.Count);
            var x = issue.Parameters.LeftX + issue.Parameters.H * j;

            var coefClearPrev = 0.0;
            var coefClearCur = 0.0;
            var coefClearNext = 0.0;
            var dClear = issue.BaseExpression.a * (u[j + 1] - 2 * u[j] + u[j - 1]) / (issue.Parameters.H * issue.Parameters.H)
                + issue.BaseExpression.b * (u[j + 1] - u[j - 1]) / (2.0 * issue.Parameters.H)
                + issue.BaseExpression.c * u[j]
                + issue.BaseExpression.f(x, t);

            var coefUnclearPrev = -issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) + issue.BaseExpression.b / (2.0 * issue.Parameters.H); ;
            var coefUnclearCur = issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) - issue.BaseExpression.c;
            var coefUnclearNext = -issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) - issue.BaseExpression.b / (2.0 * issue.Parameters.H);
            var dUnclear = issue.BaseExpression.f(x, t);


            var d = u[j] / issue.Parameters.Tau + (dClear + dUnclear) / 2.0;
            coefs[j - 1] = (coefClearPrev + coefUnclearPrev) / 2.0;
            coefs[j] = 1 / issue.Parameters.Tau + (coefClearCur + coefUnclearCur);
            coefs[j + 1] = (coefClearNext + coefUnclearNext) / 2.0;

            return new Equation(coefs, d);
        }
    }
}
