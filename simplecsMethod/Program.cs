using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra;


namespace simplecsMethod
{
    internal static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Form1 form1 = new Form1();

            programm();

            //Application.Run(form1);
        }

        static void programm()
        {
            //Matrix<double> m = Matrix<double>.Build.DenseOfArray(new double[,]{
            //    {2, 3, 1, 0, 0},
            //    {1, 0, 0, 1, 0},
            //    {1,-1, 0, 0,-1}
            //});

            //Vector<double> b = Vector<double>.Build.DenseOfArray(new double[]{6, 1, -1});

            //Vector<double> z = Vector<double>.Build.Dense(new double[] {0, -1, -2,0,0,0});

            Matrix<double> m = Matrix<double>.Build.DenseOfArray(new double[,]{
                {5,  3, -1,  0,  0},
                {2,  4,  0, -1,  0},
                {3, 11,  0,  0, -1}
            });

            Vector<double> b = Vector<double>.Build.DenseOfArray(new double[] { 30, 26, 54 });

            Vector<double> z = Vector<double>.Build.Dense(new double[] { 0, 5, 2, 0, 0, 0 });

            SimplecsTable st = new SimplecsTable(m, b, z);

            //Matrix<double> m = Matrix<double>.Build.DenseOfArray(new double[,]{
            //    {1, -4,  2, -5,  9},
            //    {0,  1, -3,  4, -5},
            //    {0,  1, -1,  1, -1}
            //});

            //Vector<double> b = Vector<double>.Build.DenseOfArray(new double[] { 3, 6, 1 });

            //Vector<double> z = Vector<double>.Build.Dense(new double[] { 0, -2, -6, 5, -1, -4 });

            //SimplecsTable st = new SimplecsTable(m, b, z);

            st.solve();
        }
    }
}
