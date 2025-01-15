using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace simplecsMethod
{
    internal class SimplecsTable
    {
        Matrix<double> table;
       
        int zRowIdx;
        int mRowIdx = -1;
        int[] basis;
        int[] others;
        HashSet<int> artifitials;

        public SimplecsTable(Matrix<double> m, Vector<double> b, Vector<double> z)
        {
            table = b.ToColumnMatrix().Append(m).Stack(z.ToRowMatrix());
            zRowIdx = table.RowCount - 1;
            basis = Enumerable.Repeat(-1, m.RowCount).ToArray();
            others = Enumerable.Range(0, m.ColumnCount + 1).ToArray();
            artifitials = new HashSet<int>();
            changeRightPart();
            int numberOfX = table.ColumnCount - 1;
            
            for (int i = 0; i < zRowIdx; i++)
            {
                int t = findBasis(i);
                if (t < 0)
                {
                    numberOfX++;
                    artifitials.Add(numberOfX);
                    basis[i] = numberOfX;
                    
                }
                else
                {
                    basis[i] = others[t];
                    others = others.Where(x => x != basis[i]).ToArray();
                    table.SetRow(i, table.Row(i) / table[i, t]);
                    double c = table[zRowIdx, t];
                    table = table.RemoveColumn(t);
                    Vector<double> Z = -table.Row(zRowIdx);
                    Z[0] = -Z[0];
                    Z += table.Row(i) * c;
                    
                    table.SetRow(zRowIdx, Z);
                }
            }
            if (artifitials.Count() !=0) 
            {
                Vector<double> M = Vector<double>.Build.Dense(table.ColumnCount);
                int count=0;
                for (int i = 0; i < basis.Length; i++) 
                {
                    if (artifitials.Contains(basis[i]))
                        M += table.Row(i);
                }
                mRowIdx = table.RowCount;
                M = -M;
                table=table.Stack(M.ToRowMatrix());
            }
            Console.WriteLine(ToString());

        }



        void changeRightPart()
        {
            for (int i = 0; i < table.RowCount - 1; i++)
            {
                if (table[i, 0] < 0)
                    table.SetRow(i, table.Row(i) * -1);
            }
        }

        int findBasis(int row)
        {
            int pos = -1;
            for (int i = 1; i < table.ColumnCount; i++)
            {
                if (table[row, i] != 0)
                {
                    int count = 0;
                    for (int j = 0; j < zRowIdx; j++)
                    {
                        if (table[j, i] == 0)
                            count++;
                    }
                    if (count == zRowIdx - 1 && table[row,i] > 0)
                        return i;
                }
            }
            return pos;
        }
        Tuple<int, int> findResolvingElement(int funcRow) 
        {
            double min=0;
            int col=-1;

            for (int i = 1; i < table.ColumnCount; i++) 
            {
                if (table[funcRow, i] < min)
                {
                    min = table[funcRow, i];
                    col = i;
                }
                   
            }

            if(col==-1)
                return new Tuple<int, int>(-1, -1);

            min = double.MaxValue;
            int row = -1;
            for (int i = 0; i < zRowIdx; i++) 
            {
                if (table[i, col] <= 0)
                    continue;
                double t = table[i, 0] / table[i, col];
                if (t < min) 
                {
                    row= i;
                    min = t;
                }
            }
            return new Tuple<int, int>(row,col);
        }

        bool isOptimal(int funcRow) 
        {
            for (int i = 1; i < table.ColumnCount; i++) 
            {
                if (table[funcRow, i] < -0.00005 ) 
                {
                    return false;
                }
            }
            return true;
        }
        void transform(int row, int col) 
        {
            Matrix<double> newTable = Matrix<double>.Build.Dense(table.RowCount, table.ColumnCount);
            int t = basis[row];
            basis[row] = others[col];
            others[col] = t;
            
            for (int i = 0; i < table.RowCount; i++) 
            {
                for (int j = 0; j < table.ColumnCount; j++) 
                {
                    if (i == row && j == col) 
                    {
                        newTable[row, col] = 1 / table[row, col];
                        continue;
                    }
                    if (i == row) 
                    {
                        newTable[i, j] = table[i, j] / table[row, col];  
                        continue;
                    }
                    if (j == col) 
                    {
                        newTable[i, j] = -table[i, j] / table[row, col]; 
                        continue;
                    }
                    newTable[i, j] =table[i, j] - (table[row, j]* table[i, col] / table[row, col]); 
                }
            }
            table= newTable;
        }

        public void solve(int funcRow) 
        {
            while (!isOptimal(funcRow)) 
            {
                (int row,int col)=findResolvingElement(funcRow);
                transform(row,col);
                if (mRowIdx > 0)
                {
                    
                    for (int i = 1; i < others.Length; i++)
                    {
                        if (artifitials.Contains(others[i]))
                        {
                            table = table.RemoveColumn(i);
                            artifitials.Remove(others[i]);
                            others = others.Where(x => x != others[i]).ToArray() ;
                            
                            
                        }

                    }
                    
                }
                Console.WriteLine(ToString());
            }
        }
        public void solve() 
        {
            if (mRowIdx > 0) 
            {
                solve(mRowIdx);
                if (artifitials.Count == 0) 
                {
                    table = table.RemoveRow(mRowIdx);
                    mRowIdx = -1;
                    
                }
                    
                
            }
            solve(zRowIdx);
        }

        public override string ToString() 
        {
            int maxWidth = table.Enumerate().Select(x=> string.Format("{0:0.0000}",x).Length).Max()+1;
            string fieldFMT = string.Format("{{0,{0}}}", maxWidth);
            string hLine = new string('-', maxWidth * (table.ColumnCount + 1)+6);
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat(fieldFMT, "");
            sb.Append(" | ");
            sb.AppendFormat(fieldFMT, "1");
            sb.Append(" | ");
            for (int i = 1; i < others.Length; i++) 
            {
                string t = string.Format("-x{0}", others[i]);
                sb.AppendFormat(fieldFMT,t);
            }
            sb.AppendLine();
            sb.AppendLine(hLine);
            for (int i = 0; i < zRowIdx; i++) 
            {
                string t = string.Format("x{0} = ", basis[i]);
                sb.AppendFormat(fieldFMT, t);
                sb.Append(" | ");
                sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[i,0]));
                sb.Append(" | ");
                for (int j = 1; j < table.ColumnCount; j++) 
                {
                    sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[i, j]));
                }
                sb.AppendLine();
            }
            sb.AppendLine(hLine);
            sb.AppendFormat(fieldFMT, "Z = ");
            sb.Append(" | ");
            sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[zRowIdx, 0]));
            sb.Append(" | ");
            for (int i = 1; i < table.ColumnCount; i++) 
            {
                sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[zRowIdx, i]));
            }
            sb.AppendLine();

            if (mRowIdx > 0) 
            {
                sb.AppendLine(hLine);
                sb.AppendFormat(fieldFMT, "M = ");
                sb.Append(" | ");
                sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[mRowIdx, 0]));
                sb.Append(" | ");
                for (int i = 1; i < table.ColumnCount; i++)
                {
                    sb.AppendFormat(fieldFMT, string.Format("{0:0.0000}", table[mRowIdx, i]));
                }
                sb.AppendLine();
            }

            return sb.ToString();
        }


        //    bool check() 
        //    {
        //        double t = 0;
        //        for (int i = 0; i < b.Length; i++) {
        //            t = t * b[i];
        //        }
        //        if(t<0)
        //            return false;
        //        if(basis())
        //        return true;
        //    }

        //    int[] bazis() 
        //    {
        //        List<int> pos = new List<int>();
        //        List<int> baz = new List<int>();
        //        for (int i = 0; i < borders.GetLength(1); i++)
        //        {
        //            int count = 0;
        //            int p = -1;
        //            for (int j = 0; j < borders.GetLength(0); j++) 
        //            {
        //                if (borders[j, i] != 0) 
        //                {
        //                    count ++;
        //                    p = j;
        //                    if (count > 1)
        //                        break;

        //                }
        //            }
        //            if (count == 1 && !pos.Contains(p)) 
        //            {
        //                baz.Add(i);
        //                pos.Add(p);
        //            }
        //        }
        //        return baz.ToArray();
        //    }

    }
}
