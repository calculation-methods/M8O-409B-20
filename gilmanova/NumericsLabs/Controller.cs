using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media.Animation;
using Numerics5;
using Matrix = System.Collections.Generic.List<System.Collections.Generic.List<double>>;

namespace NumericsLabs
{
    internal class Controller
    {
        private const double a = 2;
        private const double sigma = 0.4;
        private const int n = 50;

        private readonly Issue issue = new Issue()
        {
            AnalyticalExpression = new AnalyticalExpression((x, t) => x + Math.Exp(-Math.PI * Math.PI * a * t) * Math.Sin(Math.PI * x)),
            BaseExpression = new BaseExpression(a, 0, 0, (x, t) => 0),
            LeftBorderExpression = new LeftBorderExpression(0, 1, t => 0),
            RightBorderExpression = new RightBorderExpression(0, 1, t => 1),
            StartExpression = new StartExpression(x => x + Math.Sin(x)),
            Parameters = new IssueParameters()
            {
                LeftX = 0,
                RightX = 1,
                Tau = sigma * (1.0 / n) * (1.0 / n) / a,
                H = 1.0 / n,
                K = 5000,
                N = n
            }
        };

        private readonly Dictionary<string, CalculationMethod> calculationMethods = new Dictionary<string, CalculationMethod>()
        {
            { "Явный конечно-разностный метод", CalculationMethod.CLEAR },
            { "Невный конечно-разностный метод", CalculationMethod.UNCLEAR },
            { "Метод Кранка-Никольсона", CalculationMethod.CRANK_NIKOLSON }
        };

        private Solver solver;
        private CalculationMethod calculationMethod;


        public IEnumerable<string> CalculationNames => calculationMethods.Keys;

        public double[] RangeOfX => Enumerable.Range(0, Issue.Parameters.N + 1).Select(x => Issue.Parameters.LeftX + Issue.Parameters.H * x).ToArray();
        public double[] RangeOfT => Enumerable.Range(0, Issue.Parameters.K + 1).Select(x => Issue.Parameters.Tau * x).ToArray();

        public Issue Issue => issue;
        public Matrix Result { get; private set; }
        public Matrix AnalyticalResult { get; private set; }

        public double[] GetResultByTIndex(int t) => Result[t].ToArray();
        public double[] GetAnalyticalResultByTIndex(int t) => AnalyticalResult[t].ToArray();
        
        public double[] GetError()
        {
            var error = new double[Issue.Parameters.K + 1];

            for (int k = 0; k <= Issue.Parameters.K; k++)
            {
                error[k] = 0;

                for (int j = 0; j < Issue.Parameters.N; j++)
                {
                    var tmp = Math.Abs(AnalyticalResult[k][j] - Result[k][j]);
                    
                    if (tmp > error[k])
                    {
                        error[k] = tmp;
                    }
                }
            }

            return error;
        }

        public void SetCalculator(string calculationMethod)
        {
            this.calculationMethod = calculationMethods[calculationMethod];
        }

        public void Refresh()
        {
            solver = new Solver(Issue);
            Calculate();
        }

        private void Calculate()
        {
            Result = solver.SolveWith(calculationMethod);
            AnalyticalResult = solver.GetAnalyticalResult();
        }

    }
}
