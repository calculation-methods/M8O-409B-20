using System.Collections.Generic;
using System.Linq;

using Solution = System.Collections.Generic.List<double>;
using ReadOnlySolution = System.Collections.Generic.List<double>;
using NumericsLabs.Models;

namespace Numerics6
{

    internal class UnclearStepCalculator : IStepCalculator
    {
        private readonly Issue issue;
        private readonly IApproximation approximation;
        private readonly IApproximationByT approximationByT;

        public UnclearStepCalculator(Issue issue, IApproximation approximation, IApproximationByT approximationByT)
        {
            this.issue = issue;
            this.approximation = approximation;
            this.approximationByT = approximationByT;
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

            result[1] = approximationByT.Approximate(issue);

            for (int k = 2; k <= issue.Parameters.K; k++)
            {
                var t = issue.Parameters.Tau * (k + 1);
                result[k] = CalculateNextStep(result[k - 1], result[k - 2], t);
            }

            return result;
        }

        private Solution CalculateNextStep(ReadOnlySolution previousStep, ReadOnlySolution prepreviousStep, double t)
        {
            var equations = new List<Equation>();
            equations.Add(approximation.Approximate(issue.LeftBorderExpression, issue.Parameters, t, new LeftBorderInfo(issue.BaseExpression, previousStep[0])));
            equations.AddRange(Enumerable.Range(1, issue.Parameters.N - 1).Select(j => GetJthEquation(previousStep, prepreviousStep, t, j)));
            equations.Add(approximation.Approximate(issue.RightBorderExpression, issue.Parameters, t, new RightBorderInfo(issue.BaseExpression, previousStep[issue.Parameters.N])));

            var system = new LinearSystem(equations);
            return system.Solve3Dim();
        }

        private Equation GetJthEquation(ReadOnlySolution up, ReadOnlySolution upp, double t, int j)
        {
            var coefs = ListUtils.New<double>(up.Count);
            coefs[j - 1] = -issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) + issue.BaseExpression.b / (2.0 * issue.Parameters.H);
            coefs[j] = 1 / (issue.Parameters.Tau * issue.Parameters.Tau) + issue.BaseExpression.d / (2.0 * issue.Parameters.Tau) + 2.0 * issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) - issue.BaseExpression.c;
            coefs[j + 1] = -issue.BaseExpression.a / (issue.Parameters.H * issue.Parameters.H) - issue.BaseExpression.b / (2.0 * issue.Parameters.H);
            var x = issue.Parameters.LeftX + issue.Parameters.H * j;

            var d = 2.0 * up[j] / (issue.Parameters.Tau * issue.Parameters.Tau)
                + upp[j] * (issue.BaseExpression.d / (2.0 * issue.Parameters.Tau) - 1.0 / (issue.Parameters.Tau * issue.Parameters.Tau))
                + issue.BaseExpression.f(x, t);

            return new Equation(coefs, d);
        }
    }
}
