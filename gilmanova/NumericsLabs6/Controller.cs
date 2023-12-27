using System;
using System.Collections.Generic;
using System.Linq;
using Numerics6;
using Matrix = System.Collections.Generic.List<System.Collections.Generic.List<double>>;

namespace NumericsLabs
{
    internal class Controller
    {
        private const double a = 1;
        private const double sigma = 0.051;
        private const int n = 500;

        private readonly Issue issue = new Issue()
        {
            AnalyticalExpression = new AnalyticalExpression((x, t) => Math.Sin(x - a * t) + Math.Cos(x + a * t)),
            BaseExpression = new BaseExpression(0, a * a, 0, 0, (x, t) => 0),
            LeftBorderExpression = new LeftBorderExpression(1, -1, _ => 0),
            RightBorderExpression = new RightBorderExpression(1, -1, _ => 0),
            StartExpression = new StartExpression(x => Math.Sin(x) + Math.Cos(x)),
            StartDiffExpression = new StartDiffExpression(x => -a * (Math.Sin(x) + Math.Cos(x)),
                                                              x => -a * (Math.Cos(x) - Math.Sin(x)),
                                                              x => -a * (-Math.Sin(x) - Math.Cos(x))),
            Parameters = new IssueParameters()
            {
                LeftX = 0,
                RightX = Math.PI,
                Tau = sigma * (Math.PI / n) / a,
                H = Math.PI / n,
                K = 300,
                N = n
            }
        };


        private readonly Dictionary<string, ApproximationMethod> approximationMethods = new Dictionary<string, ApproximationMethod>()
        {
            { "Двухточечная аппроксимация первого порядка", ApproximationMethod.POINTS2_ORDER1 },
            { "Двухточечная аппроксимация второго порядка", ApproximationMethod.POINTS2_ORDER2 },
            { "Трехточечная аппроксимация второго порядка", ApproximationMethod.POINTS3_ORDER2 }
        };

        private readonly Dictionary<string, ApproximationByTMethod> approximationByTMethods = new Dictionary<string, ApproximationByTMethod>()
        {
            { "Аппроксимация по времени первого порядка", ApproximationByTMethod.ORDER1 },
            { "Аппроксимация по времени второго порядка", ApproximationByTMethod.ORDER2 },
        };

        private readonly Dictionary<string, CalculationMethod> calculationMethods = new Dictionary<string, CalculationMethod>()
        {
            { "Явный конечно-разностный метод", CalculationMethod.CLEAR },
            { "Невный конечно-разностный метод", CalculationMethod.UNCLEAR }
        };

        private Solver solver;
        private ApproximationMethod approximationMethod;
        private ApproximationByTMethod approximationByTMethod;
        private CalculationMethod calculationMethod;


        public IEnumerable<string> ApproximationNames => approximationMethods.Keys;
        public IEnumerable<string> ApproximationByTMethods => approximationByTMethods.Keys;
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

        public void SetApproximation(string approximationMethod)
        {
            this.approximationMethod = approximationMethods[approximationMethod];
        }

        public void SetApproximationByT(string approximationByTMethod)
        {
            this.approximationByTMethod = approximationByTMethods[approximationByTMethod];
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
            Result = solver.SolveWith(calculationMethod, approximationMethod, approximationByTMethod);
            AnalyticalResult = solver.GetAnalyticalResult();
        }

    }
}
