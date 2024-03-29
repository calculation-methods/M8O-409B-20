﻿using System;
using System.Collections.Generic;
using Matrix = System.Collections.Generic.List<System.Collections.Generic.List<double>>;
using Solution = System.Collections.Generic.List<double>;

namespace Numerics8
{
    internal class Solver
    {
        public Solver(Issue issue)
        {
            Issue = issue;
        }
        public Issue Issue { get; set; }

        public List<Matrix> GetAnalyticalResult()
        {
            var result = ListUtils.New<Matrix>(Issue.Parameters.K + 1);
            
            for (int k = 0; k < result.Count; k++)
            {
                result[k] = ListUtils.New<Solution>(Issue.Parameters.Ny + 1);
                for (int j = 0; j < result[k].Count; j++)
                {
                    result[k][j] = ListUtils.New<double>(Issue.Parameters.Nx + 1);

                    for (int i = 0; i < result[k][j].Count; i++)
                    {
                        double x = i * Issue.Parameters.Hx;
                        double y = j * Issue.Parameters.Hy;
                        double t = k * Issue.Parameters.Tau;
                        result[k][j][i] = Issue.AnalyticalExpression.U(x, y, t);
                    }
                }
            }

            return result;
        }

        public List<Matrix> SolveWith(CalculationMethod calculationMethod)
        {
            IStepCalculator stepCalculator;
            switch (calculationMethod)
            {
                case CalculationMethod.MPN:
                    stepCalculator = new MPNStepCalculator(Issue);
                    break;
                case CalculationMethod.MDSh:
                    stepCalculator = new MDShStepCalculator(Issue);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(calculationMethod));
            }

            return stepCalculator.Calculate();
        }
    }
}
