using System.Collections.Generic;

namespace Numerics5
{
    internal static class ListUtils
    {
        public static List<T> New<T>(int capacity) => new List<T>(new T[capacity]);
    }
}
