using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Collections.ObjectModel;
using ZedGraph;

namespace Prognostication2
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double DELTA = 0.0001;
        Int32 K, N0;
        Double[] P;
        Double dt, Pg;
        ObservableCollection<Results> ResCol = new ObservableCollection<Results>();
        public MainWindow()
        {
            InitializeComponent();
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            ExpZedGraph.GraphPane.Title.Text = "Экспоненциальный закон";
            ExpZedGraph.GraphPane.XAxis.Title.Text = "t";
            ExpZedGraph.GraphPane.YAxis.Title.Text = "P";
            ErlZedGraph.GraphPane.Title.Text = "Закон Эрланга";
            ErlZedGraph.GraphPane.XAxis.Title.Text = "t";
            ErlZedGraph.GraphPane.YAxis.Title.Text = "P";
            RelZedGraph.GraphPane.Title.Text = "Закон Рэлея";
            RelZedGraph.GraphPane.XAxis.Title.Text = "t";
            RelZedGraph.GraphPane.YAxis.Title.Text = "P";
            VejbZedGraph.GraphPane.Title.Text = "Закон Вейбулла";
            VejbZedGraph.GraphPane.XAxis.Title.Text = "t";
            VejbZedGraph.GraphPane.YAxis.Title.Text = "P";
            NormZedGraph.GraphPane.Title.Text = "Нормальный закон";
            NormZedGraph.GraphPane.XAxis.Title.Text = "t";
            NormZedGraph.GraphPane.YAxis.Title.Text = "P";
            ShortNormZedGraph.GraphPane.Title.Text = "Усеченный нормальный закон";
            ShortNormZedGraph.GraphPane.XAxis.Title.Text = "t";
            ShortNormZedGraph.GraphPane.YAxis.Title.Text = "P";
            dataGrid1.ItemsSource = ResCol;
        }

        private void KTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Int32.TryParse(KTextBox.Text, out K);
            if ( K < 1 )
                return;
            ResCol.Clear();
            for (int i = 0; i < K; i++)
            {
                ResCol.Add(new Results(i + 1));
            }
            dataGrid1.ItemsSource = ResCol;
        }

        private void StartButton_Click(object sender, RoutedEventArgs e)
        {
            Clear();
            Double.TryParse(dtTextBox.Text, out dt);
            if (dt <= 0)
            {
                MessageBox.Show("Неверно введено dt.");
                return;
            }
            Double.TryParse(PgTextBox.Text, out Pg);
            if (!(Pg >= 0 && Pg <= 1))
            {
                MessageBox.Show("Неверно введено Pg.");
                return;
            }
            Int32.TryParse(N0TextBox.Text, out N0);
            if (N0 <= 0)
            {
                MessageBox.Show("Неверно введено N0.");
                return;
            }
            for (int i = 1; ; i *= 10 )
            {
                if ((int)dt / i == 0)
                {
                    DELTA *= i;
                    break;
                }
            }
            Int32.TryParse(KTextBox.Text, out K);
            if (K < 1)
                return;
            P = new Double[K];
            for(int i = 0; i < K; i++)
            {
                P[i] = CountProbability(i);
            }
            if (ExpRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Экпоненциальное";
                DrawExpGraphics();
            }
            else if (ErlRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Эрланга";
                DrawErlGraphics();
            }
            else if (RelRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Рэлея";
                DrawRelGraphics();
            }
            else if (VejbRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Вейбулла";
                DrawVejbGraphics();
            }
            else if (NormRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Нормальный";
                DrawNormGraphics();
            }
            else if (ShortNormRadioButton.IsChecked == true)
            {
                LowLabel.Content = "Усеченный нормальный";
                DrawShortNormGraphics();
            }
        }

        int CountNi(int idx)
        {
            idx--;
            if (idx < 0)
                return N0;
            int Ni = N0;
            for (int i = 0; i < idx + 1; i++)
            {
                Ni -= ResCol[i].ni;
            }
            return Ni;
        }

        double CountProbability(int i)
        {
            int i1 = i;
            int i2 = i + 1;
            return (double)(CountNi(i1) + CountNi(i2)) / (2.0 * N0);
        }

        double CountLambda(int i)
        {
            if (i < 0)
                return 0;
            int i1 = i;
            int i2 = i + 1;
            return (double)(2.0 * ResCol[i].ni / (dt * (CountNi(i1) + CountNi(i2))));
        }

        double CountF(int i)
        {
            i--;
            return ((double)ResCol[i].ni / (N0 * dt));
        }

        void DrawExpGraphics()
        {
            double t, a1 = 0, a2 = 0;
            for (int i = 0; i < K; i++)
            {
                a1 += (double)ResCol[i].ni / (CountNi(i) + CountNi(i + 1));
                a2 += (1.0 - P[i]) / (i + 1);
            }
            a1 *= 2.0 / (K * dt);
            a2 *= 1.0 / (K * dt);
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            GraphPane pane = ExpZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowExp(a1, t));
                list2.Add(t, FuncLowExp(a2, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("aэ*", list2, System.Drawing.Color.Green, SymbolType.None);
            ExpZedGraph.AxisChange();
            // Обновляем график
            ExpZedGraph.Invalidate();

            Inf1Label.Content = "a э = " + a1.ToString();
            Inf2Label.Content = "a э* = " + a2.ToString();
            Inf3Label.Content = "T гэ = " + ((1.0 - Pg) / a1).ToString();
            Inf4Label.Content = "T оэ = " + (1.0 / a1).ToString();
            Inf5Label.Content = "σ тэ = " + (1.0 / a1).ToString();
            //D[0] = CountDExp(a1);
            //D[1] = CountDExp(a2);
            //D[2] = CountDExp(a3);
        }

        void DrawErlGraphics()
        {
            double t, a1 = 0, a2 = 0, sqr;
            for (int i = 0; i < K; i++)
            {
                sqr = CountLambda(i) - CountLambda(i - 1);/*
sqr = (CountLambda(i) - CountLambda(i - 1)) / dt;
if (sqr < 0)
    sqr *= -1;
sqr = Math.Sqrt(sqr)*/
                a1 += Math.Sqrt(sqr) / (1 - (i + 1 - 0.5) * Math.Sqrt(sqr * dt));
                a2 += Math.Sqrt((1 - P[i]) / ((i + 1) * (i + 1)));
            }
            a1 /= (K * Math.Sqrt(dt));
            a2 /= (K * dt);
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            GraphPane pane = ErlZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowErl(a1, t));
                list2.Add(t, FuncLowErl(a2, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("aэ*", list2, System.Drawing.Color.Green, SymbolType.None);
            ErlZedGraph.AxisChange();
            // Обновляем график
            ErlZedGraph.Invalidate();

            Inf1Label.Content = "a э = " + a1.ToString();
            Inf2Label.Content = "a э* = " + a2.ToString();
            Inf3Label.Content = "T гэ = " + (Math.Sqrt(1 - Pg) / a1).ToString();
            Inf4Label.Content = "T оэ = " + (2.0 / a1).ToString();
            Inf5Label.Content = "σ тэ = " + (Math.Sqrt(2.0) / a1).ToString();

            //D[3] = CountDErl(a1);
            //D[4] = CountDErl(a2);
        }

        void DrawRelGraphics()
        {
            double t, a1 = 0, a2 = 0, a3 = 0;
            for (int i = 0; i < K; i++)
            {
                a1 += ResCol[i].ni * (i + 1 - 0.5) / (CountNi(i) + CountNi(i + 1));
                a3 += (1.0 - P[i]) / ((i + 1) * (i + 1));
            }
            a1 *= 12 / (K * (4 * K * K - 1) * dt * dt);
            a2 = (CountLambda(K - 1) - CountLambda(0)) / (2 * K * dt);
            a3 /= (K * dt * dt);
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            GraphPane pane = RelZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            PointPairList list2 = new PointPairList();
            PointPairList list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowRel(a1, t));
                list2.Add(t, FuncLowRel(a2, t));
                list3.Add(t, FuncLowRel(a3, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            LineItem Curve2 = pane.AddCurve("aэ2", list2, System.Drawing.Color.Green, SymbolType.None);
            LineItem Curve3 = pane.AddCurve("aэ*", list3, System.Drawing.Color.Blue, SymbolType.None);
            RelZedGraph.AxisChange();
            // Обновляем график
            RelZedGraph.Invalidate();

            Inf1Label.Content = "a э1 = " + a1.ToString();
            Inf2Label.Content = "T гэ1 = " + (Math.Sqrt((1.0 - Pg) / a1)).ToString();
            Inf3Label.Content = "T оэ1 = " + (Math.Sqrt(Math.PI / (4 * a1))).ToString();
            Inf4Label.Content = "σ тэ1 = " + (Math.Sqrt((4.0 - Math.PI) / (4 * a1))).ToString();
            Inf5Label.Content = "a э2 = " + a2.ToString();
            Inf6Label.Content = "T гэ2 = " + (Math.Sqrt((1.0 - Pg) / a2)).ToString();
            Inf7Label.Content = "T оэ2 = " + (Math.Sqrt(Math.PI / (4 * a2))).ToString();
            Inf8Label.Content = "σ тэ2 = " + (Math.Sqrt((4.0 - Math.PI) / (4 * a2))).ToString();
            Inf9Label.Content = "a э* = " + a3.ToString();
            Inf10Label.Content = "D1 = " + CountDRel(a1).ToString();
            Inf11Label.Content = "D2 = " + CountDRel(a2).ToString();
            if (CountDRel(a1) < CountDRel(a2))
                Inf10Label.Foreground = Brushes.Red;
            else if (CountDRel(a1) > CountDRel(a2))
                Inf11Label.Foreground = Brushes.Red;
            //D[6] = CountDRel(a1);
            //D[7] = CountDRel(a2);
            //D[8] = CountDRel(a3);
        }

        void DrawVejbGraphics()
        {
            double t, x, y, a = 0, b = 0, c = 0, d = 0, e = 0;
            for (int i = 0; i < K; i++)
            {
                a += Math.Log(CountLambda(i)); // i+1
                b += Math.Log(i + 1 - 0.5) + Math.Log(dt);
                c += Math.Pow(Math.Log(i + 1 - 0.5) + Math.Log(dt), 2);
                e += Math.Log(CountLambda(i)) * (Math.Log(i + 1 - 0.5) + Math.Log(dt));
            }
            d = b * b;
            x = (a  * c - b * e) / (K * c - b * b);
            y = (a * b - K * e) / (b * b - K * c);
            double B = 1 + y;
            double A = Math.Exp(x) / B;
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            GraphPane pane = VejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowVejb(A, B, t));
            }//
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            VejbZedGraph.AxisChange();
            // Обновляем график
            VejbZedGraph.Invalidate();

            double Tg = Math.Pow((1.0 - Pg) / A, (1.0 / B));

            if (A.ToString() != "NaN")
                Inf1Label.Content = "a э = " + A.ToString();
            if (B.ToString() != "NaN")
                Inf2Label.Content = "b э = " + B.ToString();
            if (Tg.ToString() != "NaN")
            {
                Inf3Label.Content = "Т гэ = " + Tg.ToString();
                Inf4Label.Content = "σ находится по фомулам 18, д-е.";
            }
        }

        void DrawNormGraphics()
        {
            double t, T4 = 0, q4 = 0, a, b = 0, g = 1, Y, X;
            a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i <= K - 1; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            Y = (double)(K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            T4 = X * dt;
            q4 = Math.Sqrt(Y ) * dt;
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            GraphPane pane = NormZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list4.Add(t, FuncLowNorm(T4, q4, t));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve4 = pane.AddCurve("a54", list4, System.Drawing.Color.Green, SymbolType.None);
            NormZedGraph.AxisChange();
            // Обновляем график
            NormZedGraph.Invalidate();
            if (T4.ToString() != "NaN")
                Inf1Label.Content = "T нэ = " + T4.ToString();
            if (q4.ToString() != "NaN")
            {
                Inf2Label.Content = "σ нэ = " + q4.ToString();
                Inf3Label.Content = "T гэ находится по формуле 19, ж.";
            }
            if (q4.ToString() != "NaN")
                Inf4Label.Content = "σ = " + q4.ToString();
        }

        void DrawShortNormGraphics()
        {
            double t, C, T = 0, q = 0, a, b = 0, g = 1, Y, X;
            a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i < K; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            Y = (double)(K - 1) * (0.5 * K * a - g) / (a * a - b * (K - 1));
            X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            T = X * dt;
            q = Math.Sqrt(Y) * dt;
            C = 1.0 / (0.5 + FLaplas(T / q));
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            GraphPane pane = ShortNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            PointPairList list0 = new PointPairList();
            PointPairList list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                t = (ResCol[i].i - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowShortNorm(T, q, t, C));
            }
            LineItem Curve0 = pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            LineItem Curve1 = pane.AddCurve("a6", list1, System.Drawing.Color.Black, SymbolType.None);
            ShortNormZedGraph.AxisChange();
            // Обновляем график
            ShortNormZedGraph.Invalidate();

            double k = C / Math.Sqrt(2 * Math.PI) * Math.Exp(T * T / (-2 * q * q));
            double q1 = q * Math.Sqrt(1 - k * k + k * T / q);

            if (T.ToString() != "NaN")
                Inf1Label.Content = "T ун.э = " + T.ToString();
            if (q.ToString() != "NaN")
                Inf2Label.Content = "σ ун.э = " + q.ToString();
            if (C.ToString() != "NaN")
                Inf3Label.Content = "Cy = " + C.ToString();
            if (q1.ToString() != "NaN")
            {
                Inf5Label.Content = "σ = " + q1.ToString();
                Inf4Label.Content = "Tгэ находится по формуле 20, е.";
            }
        }

        double FuncLowExp(double a, double t)
        {
            return Math.Exp(-a * t);
        }

        double FuncLowErl(double a, double t)
        {
            return (1 + a * t) * Math.Exp(-a * t);
        }

        double FuncLowRel(double a, double t)
        {
            return Math.Exp(-a * t * t);
        }

        double FuncLowVejb(double a, double b, double t)
        {
            return Math.Exp(-a * Math.Pow(t, b));
        }

        double FuncNorm(double T, double q, double t)
        {
            return (1 / (q * Math.Sqrt(2 * Math.PI))) * Math.Exp(-Math.Pow(t - T, 2) / (2 * Math.Pow(q, 2)));
        }

        double FuncLowNorm(double T, double q, double t)
        {
            return 1 - IntegralForNorm(0, t, T, q);
        }

        double FuncLowShortNorm(double T, double q, double t, double C)
        {
            return C * (1 - IntegralForNorm(0, t, T, q));
        }

 /*       double FuncLowShortNorm(double T, double q, double t, double C)
        {
            double res = C * (0.5 - FLaplas((t - T) / q));
            return res;
        }
*/
        double IntegralForNorm(double lim1, double lim2, double T, double q)
        {
            double result = 0;
            double s = DELTA;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += FuncNorm(T, q, i) * 2 * s;
            }
            return result;
        }
/*
        double CountDExp(double a)
        {
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = a * Math.Exp(-a * (i - 0.5) * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDErl(double a)
        {
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = a * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }*/

        double CountDRel(double a)
        {/*
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = 2 * a * (i - 0.5) * dt * Math.Exp(-a * (i - 0.5) * (i - 0.5) * dt * dt) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;*/
            Double D = 0;
            Double temp;
            for (int i = 1; i <= K; i++)
            {
                temp = Math.Exp(-a * (i - 0.5) * dt * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }
        /*

        double CountDVejb(double a, double b)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = a * b * Math.Pow(t, b - 1) * Math.Exp(-a * Math.Pow(t, b)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDNorm(double q, double T)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = 1 / (q * Math.Sqrt(2 * Math.PI)) * Math.Exp(-(t - T) * (t - T) / (2 * q * q)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }

        double CountDShortNorm(double q, double T, double C)
        {
            Double D = 0;
            Double temp;
            Double t;
            for (int i = 1; i <= K; i++)
            {
                t = (i - 0.5) * dt;
                temp = C / (q * Math.Sqrt(2 * Math.PI)) * Math.Exp(-(t - T) * (t - T) / (2 * q * q)) - ResCol[i - 1].ni / (N0 * dt);
                D += Math.Pow(temp, 2);
            }
            D /= K;
            return D;
        }*/

        double Factorial(int x)
        {
            double y = 1;
            if (x == 0)
                return 0;
            for (int i = x; i != 0; i--)
            {
                y *= i;
            }
                return y;
        }

        double FLaplas(double x)
        {
            double res = 1.0 / Math.Sqrt(2 * Math.PI) * IntegralForLaplas(0, x, x);
            return res;
        }

        double IntegralForLaplas(double lim1, double lim2, double x)
        {/*
            if (lim2 < 0)
            {
                double temp = lim2;
                lim2 = lim1;
                lim1 = temp;
            }*/
            double result = 0;
            double s = DELTA;
            for (double i = lim1 + s; i <= lim2 - s; i += 2 * s)
            {
                result += Math.Exp(-(x * x) / 2) * s * 2;
            }
            return result;
        }

        void Clear()
        {
            ZedGraphControl ExpZedGraph = ExpWFH.Child as ZedGraphControl;
            ZedGraphControl ErlZedGraph = ErlWFH.Child as ZedGraphControl;
            ZedGraphControl RelZedGraph = RelWFH.Child as ZedGraphControl;
            ZedGraphControl VejbZedGraph = VejbWFH.Child as ZedGraphControl;
            ZedGraphControl NormZedGraph = NormWFH.Child as ZedGraphControl;
            ZedGraphControl ShortNormZedGraph = ShortNormWFH.Child as ZedGraphControl;
            GraphPane pane = ExpZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = ErlZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = RelZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = VejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = NormZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = ShortNormZedGraph.GraphPane;
            pane.CurveList.Clear();

            LowLabel.Content = "";
            Inf1Label.Content = "";
            Inf2Label.Content = "";
            Inf3Label.Content = "";
            Inf4Label.Content = "";
            Inf5Label.Content = "";
            Inf6Label.Content = "";
            Inf7Label.Content = "";
            Inf8Label.Content = "";
            Inf9Label.Content = "";
            Inf10Label.Content = "";
            Inf11Label.Content = "";
            Inf10Label.Foreground = Brushes.Black;
            Inf11Label.Foreground = Brushes.Black;
        }
    }
}
