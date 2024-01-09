using Numerics6;

namespace NumericsLabs.Models
{
    internal class LeftBorderInfo : ApproximationInfo
    {
        public LeftBorderInfo(BaseExpression expression, double previousU0Value)
        {
            Expression = expression;
            PreviousU0Value = previousU0Value;
        }

        public BaseExpression Expression { get; }

        public double PreviousU0Value { get; }
    }
}
