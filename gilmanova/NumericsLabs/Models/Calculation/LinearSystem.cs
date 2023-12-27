using System;
using System.Collections.Generic;
using System.Linq;
using Matrix = System.Collections.Generic.List<System.Collections.Generic.List<double>>;
using Column = System.Collections.Generic.List<double>;
using Solution = System.Collections.Generic.List<double>;

namespace Numerics5
{
    internal class LinearSystem
    {
        private readonly Matrix A;
        private readonly Column D;

        public LinearSystem(List<Equation> equations)
        {
            A = new Matrix(equations.Select(eq => eq.Coefs.ToList()));
            D = new Column(equations.Select(eq => eq.Result));
        }

        public Solution SolveOneCoef()
        {
            int n = D.Count;
            var result = ListUtils.New<double>(n);

            for (int i = 0; i < n; i++)
            {
                result[i] = D[i] / A[i][i];
            }

            return result;
        }

        public Solution Solve3Dim()
        {
            int n = D.Count;

            Make3Diag();
            var abc = Decomposite(ListToArray2D(A));
            var pq = FindPQ(abc, D.ToArray());
            var xn = (D[n - 1] - abc.a[n - 1] * pq.Q[n - 2]) / (abc.a[n - 1] * pq.P[n - 2] + abc.b[n - 1]);
            var x = FindXWithPQ(pq, xn);

            return x.ToList();
        }

        private void Make3Diag()
        {
            int n = D.Count;

            if (A[0][2] != 0)
            {
                var koef = -A[0][2] / A[1][2];

                for (var i = 0; i < 3; i++)
                {
                    A[0][i] += koef * A[1][i];
                }

                D[0] += koef * D[1];
            }

            if (A[n - 1][n - 3] != 0)
            {
                var koef = -A[n - 1][n - 3] / A[n - 2][n - 3];

                for (var i = 1; i <= 3; i++)
                {
                    A[n - 1][n - i] += koef * A[n - 2][n - i];
                }

                D[n - 1] += koef * D[n - 2];
            }
        }

        private T[,] ListToArray2D<T>(List<List<T>> lists)
        {
            T[,] arrays = new T[lists.Count, lists[0].Count];
            for (int i = 0; i < lists.Count; i++)
            {
                for (int j = 0; j < lists[i].Count; j++)
                {
                    arrays[i, j] = lists[i][j];
                }

            }
            return arrays;
        }

        private static (double[] P, double[] Q) FindPQ((double[] a, double[] b, double[] c) abc, double[] d)
        {
            int n = abc.a.Length;
            var P = new double[n];
            var Q = new double[n];
            P[0] = -abc.c[0] / abc.b[0];
            Q[0] = d[0] / abc.b[0];

            for (int i = 1; i < n - 1; i++)
            {
                P[i] = -abc.c[i] / (abc.a[i] * P[i - 1] + abc.b[i]);
                Q[i] = (d[i] - abc.a[i] * Q[i - 1]) / (abc.a[i] * P[i - 1] + abc.b[i]);
            }

            return (P, Q);
        }

        private static double[] FindXWithPQ((double[] P, double[] Q) pq, double xn)
        {
            var n = pq.P.Length;
            var x = new double[n];

            x[n - 1] = xn;
            for (int i = n - 2; i >= 0; i--)
            {
                x[i] = pq.P[i] * x[i + 1] + pq.Q[i];
            }
            return x;
        }

        private static (double[] a, double[] b, double[] c) Decomposite(double[,] A)
        {
            var n = A.GetLength(0);
            var a = new double[n];
            var b = new double[n];
            var c = new double[n];
            b[0] = A[0, 0];
            c[0] = A[0, 1];
            a[n - 1] = A[n - 1, n - 2];
            b[n - 1] = A[n - 1, n - 1];

            for (int i = 1; i < n - 1; i++)
            {
                a[i] = A[i, i - 1];
                b[i] = A[i, i];
                c[i] = A[i, i + 1];
            }

            return (a, b, c);
        }

    }
}
