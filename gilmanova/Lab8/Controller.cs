using Numerics8;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsLabs3
{
    internal class Controller
    {
        private const int nx = 50;
        private const int ny = 50;
        private const double a = 1;
        private const double mu1 = 1.0;
        private const double mu2 = 1.0;

        private readonly Dictionary<string, CalculationMethod> calculationMethods = new Dictionary<string, CalculationMethod>()
        {
            { "Метод переменных направлений", CalculationMethod.MPN },
            { "Метод дробных шагов", CalculationMethod.MDSh },
        };

        private readonly Issue issue = new()
        {
            AnalyticalExpression = new AnalyticalExpression((x, y, t) => Math.Cos(mu1 * x) * Math.Cos(mu2 * y) * Math.Exp(-(mu1 * mu1 + mu2 * mu2) * a * t)),
            BaseExpression = new BaseExpression(a, a, 0, 0, 0, (x, y, t) => 0),
            LeftBorderExpression = new LeftBorderExpression(0, 1, (y, t) => Math.Cos(mu2 * y) * Math.Exp(-(mu1 * mu1 + mu2 * mu2) * a * t)),
            RightBorderExpression = new RightBorderExpression(0, 1, (y, t) => 0),
            TopBorderExpression = new TopBorderExpression(0, 1, (x, t) => Math.Cos(mu1 * x) * Math.Exp(-(mu1 * mu1 + mu2 * mu2) * a * t)),
            BottomBorderExpression = new BottomBorderExpression(0, 1, (x, t) => 0),
            StartExpression = new StartExpression((x, y) => Math.Cos(mu1 * x) * Math.Cos(mu2 * y)),
            Parameters = new IssueParameters()
            {
                RightX = Math.PI / 2.0 * mu1,
                TopY = Math.PI / 2.0 * mu2,
                Ny = ny,
                Nx = nx,
                Hx = Math.PI / 2.0 * mu1 / nx,
                Hy = Math.PI / 2.0 * mu2 / ny,
                Tau = 0.01,
                K = 50
            }
        };
        
        private Solver solver;
        private CalculationMethod calculationMethod;

        public double[] RangeOfX => Enumerable.Range(0, Issue.Parameters.Nx + 1).Select(x => Issue.Parameters.Hx * x).ToArray();
        public double[] RangeOfY => Enumerable.Range(0, Issue.Parameters.Ny + 1).Select(x => Issue.Parameters.Hy * x).ToArray();
        public double[] RangeOfT => Enumerable.Range(0, Issue.Parameters.K + 1).Select(x => Issue.Parameters.Tau * x).ToArray();

        public IEnumerable<string> CalculationNames => calculationMethods.Keys;

        public Issue Issue => issue;
        public List<double[,]> Result { get; private set; }
        public List<double[,]> AnalyticalResult { get; private set; }

        public double[,] GetResultByTIndex(int t) => Result[t];
        public double[,] GetAnalyticalResultByTIndex(int t) => AnalyticalResult[t];

        public double[] GetResultByTandYIndex(int t, int y)
        {
            var tresult = GetResultByTIndex(t);
            var res = new double[tresult.GetLength(1)];
            for (int i = 0; i < res.GetLength(0); i++)
            {
                res[i] = tresult[y, i];
            }

            return res;
        }

        public double[] GetAnalyticalResultByTandYIndex(int t, int y)
        {
            var tresult = GetAnalyticalResultByTIndex(t);
            var res = new double[tresult.GetLength(1)];
            for (int i = 0; i < res.GetLength(0); i++)
            {
                res[i] = tresult[y, i];
            }

            return res;
        }


        public double[] GetError()
        {
            var error = new double[Result.Count];

            for (int k = 0; k < Result.Count; k++)
            {
                error[k] = 0;
                for (int j = 0; j <= Issue.Parameters.Ny; j++)
                {
                    for (int i = 0; i <= Issue.Parameters.Nx; i++)
                    {
                        var tmp = Math.Abs(AnalyticalResult[k][j, i] - Result[k][j, i]);

                        if (tmp > error[k])
                        {
                            error[k] = tmp;
                        }
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
            Result = solver.SolveWith(calculationMethod).Select(el => ListUtils.To2DArray(el)).ToList();
            AnalyticalResult = solver.GetAnalyticalResult().Select(el => ListUtils.To2DArray(el)).ToList();
        }
    }
}
