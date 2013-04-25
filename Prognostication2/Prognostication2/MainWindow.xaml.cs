using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Collections.ObjectModel;
using ZedGraph;

namespace Prognostication2
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private double DELTA;
        private int K;
        private int N0;
        private double[] P;
        private double dt;
        private double Pg;
        private readonly ObservableCollection<Results> _resCol = new ObservableCollection<Results>();
        
        public MainWindow()
        {
            InitializeComponent();

            var expZedGraph = (ZedGraphControl)ExpWFH.Child;
            expZedGraph.GraphPane.Title.Text = "Экспоненциальный закон";
            expZedGraph.GraphPane.XAxis.Title.Text = "t";
            expZedGraph.GraphPane.YAxis.Title.Text = "P";

            var erlZedGraph = (ZedGraphControl)ErlWFH.Child;
            erlZedGraph.GraphPane.Title.Text = "Закон Эрланга";
            erlZedGraph.GraphPane.XAxis.Title.Text = "t";
            erlZedGraph.GraphPane.YAxis.Title.Text = "P";

            var relZedGraph = (ZedGraphControl)RelWFH.Child;
            relZedGraph.GraphPane.Title.Text = "Закон Рэлея";
            relZedGraph.GraphPane.XAxis.Title.Text = "t";
            relZedGraph.GraphPane.YAxis.Title.Text = "P";

            var vejbZedGraph = (ZedGraphControl)VejbWFH.Child;
            vejbZedGraph.GraphPane.Title.Text = "Закон Вейбулла";
            vejbZedGraph.GraphPane.XAxis.Title.Text = "t";
            vejbZedGraph.GraphPane.YAxis.Title.Text = "P";

            var normZedGraph = (ZedGraphControl)NormWFH.Child;
            normZedGraph.GraphPane.Title.Text = "Нормальный закон";
            normZedGraph.GraphPane.XAxis.Title.Text = "t";
            normZedGraph.GraphPane.YAxis.Title.Text = "P";

            var shortNormZedGraph = (ZedGraphControl)ShortNormWFH.Child;
            shortNormZedGraph.GraphPane.Title.Text = "Усеченный нормальный закон";
            shortNormZedGraph.GraphPane.XAxis.Title.Text = "t";
            shortNormZedGraph.GraphPane.YAxis.Title.Text = "P";

            dataGrid1.ItemsSource = _resCol;
        }

        private void KTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            Int32.TryParse(KTextBox.Text, out K);
            if ( K < 1 )
                return;
            _resCol.Clear();
            for (int i = 0; i < K; i++)
            {
                _resCol.Add(new Results(i + 1));
            }
            dataGrid1.ItemsSource = _resCol;
        }

        private void StartButton_Click(object sender, RoutedEventArgs e)
        {
            try
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

                double ni = 0;
                for (int i = 0; i < _resCol.Count; i++)
                {
                    ni += _resCol[i].Ni;
                }
                if (ni != N0)
                {
                    MessageBox.Show("Сумма ni не равна N0.");
                    return;
                }

                DELTA = 0.0001;
                for (int i = 1; ; i *= 10)
                {
                    if ((int)dt / i == 0)
                    {
                        DELTA *= i;
                        break;
                    }
                }

                Int32.TryParse(KTextBox.Text, out K);
                if (K < 1)
                {
                    MessageBox.Show("Неправильно указано число интервалов m.");
                    return;
                }

                P = new Double[K];
                for (int i = 0; i < K; i++)
                {
                    P[i] = CountProbability(i);
                }

                if (ExpRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Экпоненциальное";

                    ExpTabItem.IsEnabled = true;
                    ErlTabItem.IsEnabled = false;
                    RelTabItem.IsEnabled = false;
                    VejbTabItem.IsEnabled = false;
                    NormTabItem.IsEnabled = false;
                    ShortNormTabItem.IsEnabled = false;

                    DrawExpGraphics();
                }
                else if (ErlRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Эрланга";

                    ExpTabItem.IsEnabled = false;
                    ErlTabItem.IsEnabled = true;
                    RelTabItem.IsEnabled = false;
                    VejbTabItem.IsEnabled = false;
                    NormTabItem.IsEnabled = false;
                    ShortNormTabItem.IsEnabled = false;

                    DrawErlGraphics();
                }
                else if (RelRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Рэлея";

                    ExpTabItem.IsEnabled = false;
                    ErlTabItem.IsEnabled = false;
                    RelTabItem.IsEnabled = true;
                    VejbTabItem.IsEnabled = false;
                    NormTabItem.IsEnabled = false;
                    ShortNormTabItem.IsEnabled = false;

                    DrawRelGraphics();
                }
                else if (VejbRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Вейбулла";

                    ExpTabItem.IsEnabled = false;
                    ErlTabItem.IsEnabled = false;
                    RelTabItem.IsEnabled = false;
                    VejbTabItem.IsEnabled = true;
                    NormTabItem.IsEnabled = false;
                    ShortNormTabItem.IsEnabled = false;

                    DrawVejbGraphics();
                }
                else if (NormRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Нормальный";

                    ExpTabItem.IsEnabled = false;
                    ErlTabItem.IsEnabled = false;
                    RelTabItem.IsEnabled = false;
                    VejbTabItem.IsEnabled = false;
                    NormTabItem.IsEnabled = true;
                    ShortNormTabItem.IsEnabled = false;

                    DrawNormGraphics();
                }
                else if (ShortNormRadioButton.IsChecked == true)
                {
                    LowLabel.Content = "Усеченный нормальный";

                    ExpTabItem.IsEnabled = false;
                    ErlTabItem.IsEnabled = false;
                    RelTabItem.IsEnabled = false;
                    VejbTabItem.IsEnabled = false;
                    NormTabItem.IsEnabled = false;
                    ShortNormTabItem.IsEnabled = true;

                    DrawShortNormGraphics();
                }
                tabControl1.SelectedIndex = 7;
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
                Clear();
            }
        }

        private double CountNi(int idx)
        {
            idx--;
            if (idx < 0)
                return N0;
            double Ni = N0;
            for (int i = 0; i < idx + 1; i++)
            {
                Ni -= _resCol[i].Ni;
            }
            return Ni;
        }

        private double CountProbability(int i)
        {
            int i1 = i;
            int i2 = i + 1;
            return (CountNi(i1) + CountNi(i2)) / (2.0 * N0);
        }

        private double CountLambda(int i)
        {
            if (i < 0)
                return 0;
            int i1 = i;
            int i2 = i + 1;
            return 2.0 * _resCol[i].Ni / (dt * (CountNi(i1) + CountNi(i2)));
        }

        private double CountF(int i)
        {
            i--;
            return _resCol[i].Ni / (N0 * dt);
        }

        private void DrawExpGraphics()
        {
            double a1 = 0, a2 = 0;
            for (int i = 0; i < K; i++)
            {
                a1 += _resCol[i].Ni / (CountNi(i) + CountNi(i + 1));
                a2 += (1.0 - P[i]) / (i + 1);
            }
            a1 *= 2.0 / (K * dt);
            a2 *= 1.0 / (K * dt);
            var expZedGraph = (ZedGraphControl) ExpWFH.Child;
            GraphPane pane = expZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowExp(a1, t));
                list2.Add(t, FuncLowExp(a2, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("aэ*", list2, System.Drawing.Color.Green, SymbolType.None);
            expZedGraph.AxisChange();
            // Обновляем график
            expZedGraph.Invalidate();

            Inf1Label.Content = "a э = " + a1.ToString();
            Inf2Label.Content = "a э* = " + a2.ToString();
            Inf3Label.Content = "T гэ = " + ((1.0 - Pg) / a1).ToString();
            Inf4Label.Content = "T оэ = " + (1.0 / a1).ToString();
            Inf5Label.Content = "σ тэ = " + (1.0 / a1).ToString();
            //D[0] = CountDExp(a1);
            //D[1] = CountDExp(a2);
            //D[2] = CountDExp(a3);
        }

        private void DrawErlGraphics()
        {
            double a1 = 0, a2 = 0, a3 = 0;
            for (int i = 0; i < K; i++)
            {
                double sqr = CountLambda(i) - CountLambda(i - 1);
                a1 += Math.Sqrt(sqr) / (1 - (i + 1 - 0.5) * Math.Sqrt(sqr * dt));
                a2 += Math.Sqrt((1 - P[i]) / ((i + 1) * (i + 1)));
            }
            a1 /= (K * Math.Sqrt(dt));
            a2 /= (K * dt);
            for (int i = 0; i < K - 1; i++)
            {
                a3 += (CountLambda(i) / (i + 1 - 0.5) - (CountLambda(i + 1) - CountLambda(i))) / ((i + 1 + 0.5) * (CountLambda(i + 1) - CountLambda(i)));
            }
            a3 /= ((K - 1) * dt);
            var erlZedGraph = (ZedGraphControl) ErlWFH.Child;
            GraphPane pane = erlZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowErl(a1, t));
                list2.Add(t, FuncLowErl(a2, t));
                list3.Add(t, FuncLowErl(a3, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("aэ*", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("a2", list2, System.Drawing.Color.Blue, SymbolType.None);
            erlZedGraph.AxisChange();
            // Обновляем график
            erlZedGraph.Invalidate();

            Inf1Label.Content = "a э = " + a1.ToString();
            Inf2Label.Content = "a э* = " + a2.ToString();
            Inf3Label.Content = "T гэ = " + (Math.Sqrt(1 - Pg) / a1).ToString();
            Inf4Label.Content = "T оэ = " + (2.0 / a1).ToString();
            Inf5Label.Content = "σ тэ = " + (Math.Sqrt(2.0) / a1).ToString();

            //D[3] = CountDErl(a1);
            //D[4] = CountDErl(a2);
        }

        private void DrawRelGraphics()
        {
            double a1 = 0, a3 = 0;
            for (int i = 0; i < K; i++)
            {
                a1 += _resCol[i].Ni * (i + 1 - 0.5) / (CountNi(i) + CountNi(i + 1));
                a3 += (1.0 - P[i]) / ((i + 1) * (i + 1));
            }
            a1 *= 12 / (K * (4 * K * K - 1) * dt * dt);
            double a2 = (CountLambda(K - 1) - CountLambda(0)) / (2 * K * dt);
            a3 /= (K * dt * dt);
            var relZedGraph = (ZedGraphControl) RelWFH.Child;
            GraphPane pane = relZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            var list2 = new PointPairList();
            var list3 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowRel(a1, t));
                list2.Add(t, FuncLowRel(a2, t));
                list3.Add(t, FuncLowRel(a3, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("aэ1", list1, System.Drawing.Color.Black, SymbolType.None);
            pane.AddCurve("aэ2", list2, System.Drawing.Color.Green, SymbolType.None);
            pane.AddCurve("aэ*", list3, System.Drawing.Color.Blue, SymbolType.None);
            relZedGraph.AxisChange();
            // Обновляем график
            relZedGraph.Invalidate();

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

        private void DrawVejbGraphics()
        {
            double a = 0, b = 0, c = 0, e = 0;
            for (int i = 0; i < K; i++)
            {
                a += Math.Log(CountLambda(i));
                b += Math.Log(i + 1 - 0.5) + Math.Log(dt);
                c += Math.Pow(Math.Log(i + 1 - 0.5) + Math.Log(dt), 2);
                e += Math.Log(CountLambda(i)) * (Math.Log(i + 1 - 0.5) + Math.Log(dt));
            }
            double x = (a  * c - b * e) / (K * c - b * b);
            double y = (a * b - K * e) / (b * b - K * c);
            double B = 1 + y;
            double A = Math.Exp(x) / B;
            var vejbZedGraph = (ZedGraphControl) VejbWFH.Child;
            GraphPane pane = vejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowVejb(A, B, t));
            }//
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a21", list1, System.Drawing.Color.Black, SymbolType.None);
            vejbZedGraph.AxisChange();
            // Обновляем график
            vejbZedGraph.Invalidate();

            double Tg = Math.Pow((1.0 - Pg) / A, (1.0 / B));

            if (!double.IsNaN(A))
            {
                Inf1Label.Content = "a э = " + A.ToString();
            }
            if (!double.IsNaN(B))
            {
                Inf2Label.Content = "b э = " + B.ToString();
            }
            if (!double.IsNaN(Tg))
            {
                Inf3Label.Content = "Т гэ = " + Tg.ToString();
                Inf4Label.Content = "σ тэ находится по фомулам 18, д-е.";
            }
        }

        private void DrawNormGraphics()
        {
            double b = 0, g = 1;
            double a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i <= K - 1; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            double Y = (K - 1) * (0.5 * K * a - g) / Math.Abs(a * a - b * (K - 1));
            double X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            double T4 = X * dt;
            double q4 = Math.Sqrt(Y ) * dt;
            var normZedGraph = (ZedGraphControl) NormWFH.Child;
            GraphPane pane = normZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list4 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list4.Add(t, FuncLowNorm(T4, q4, t));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a54", list4, System.Drawing.Color.Green, SymbolType.None);
            normZedGraph.AxisChange();
            // Обновляем график
            normZedGraph.Invalidate();
            if (!double.IsNaN(T4))
            {
                Inf1Label.Content = "T нэ = " + T4.ToString();
            }
            if (!double.IsNaN(q4))
            {
                Inf2Label.Content = "σ нэ = " + q4.ToString();
                Inf3Label.Content = "T гэ находится по формуле 19, ж.";
                Inf4Label.Content = "σ тэ = " + q4.ToString();
            }
        }

        private void DrawShortNormGraphics()
        {
            double b = 0, g = 1;
            double a = Math.Log(CountF(1) / CountF(K));
            for (int i = 1; i < K; i++)
            {
                b += Math.Pow(Math.Log(CountF(i) / CountF(i + 1)), 2);
                g *= Math.Pow(CountF(i) / CountF(i + 1), i);
            }
            g = Math.Log(g);
            double Y = (K - 1) * (0.5 * K * a - g) / Math.Abs(a * a - b * (K - 1));
            double X = (0.5 * K * b * (K - 1) - g * a) / (b * (K - 1) - a * a);
            double T = X * dt;
            double q = Math.Sqrt(Y) * dt;
            double C = 1.0 / (0.5 + FLaplas(T / q));
            var shortNormZedGraph = (ZedGraphControl) ShortNormWFH.Child;
            GraphPane pane = shortNormZedGraph.GraphPane;
            pane.CurveList.Clear();
            var list0 = new PointPairList();
            var list1 = new PointPairList();
            // Заполняем список точек
            for (int i = 0; i < K; i++)
            {
                double t = (_resCol[i].I - 0.5) * dt;
                list0.Add(t, P[i]);
                list1.Add(t, FuncLowShortNorm(T, q, t, C));
            }
            pane.AddCurve("Экспериментально", list0, System.Drawing.Color.Red, SymbolType.None);
            pane.AddCurve("a6", list1, System.Drawing.Color.Black, SymbolType.None);
            shortNormZedGraph.AxisChange();
            // Обновляем график
            shortNormZedGraph.Invalidate();

            double k = C / Math.Sqrt(2 * Math.PI) * Math.Exp(T * T / (-2 * q * q));
            double q1 = q * Math.Sqrt(1 - k * k + k * T / q);

            if (!double.IsNaN(T))
            {
                Inf1Label.Content = "T ун.э = " + T.ToString();
            }
            if (!double.IsNaN(q))
            {
                Inf2Label.Content = "σ ун.э = " + q.ToString();
            }
            if (!double.IsNaN(C))
            {
                Inf3Label.Content = "Cyэ = " + C.ToString();
            }
            if (!double.IsNaN(q1))
            {
                Inf5Label.Content = "σ тэ = " + q1.ToString();
                Inf4Label.Content = "Tгэ находится по формуле 20, е.";
            }
        }

        private double FuncLowExp(double a, double t)
        {
            return Math.Exp(-a * t);
        }

        private double FuncLowErl(double a, double t)
        {
            return (1 + a * t) * Math.Exp(-a * t);
        }

        private double FuncLowRel(double a, double t)
        {
            return Math.Exp(-a * t * t);
        }

        private double FuncLowVejb(double a, double b, double t)
        {
            return Math.Exp(-a * Math.Pow(t, b));
        }

        private double FuncNorm(double T, double q, double t)
        {
            return (1 / (q * Math.Sqrt(2 * Math.PI))) * Math.Exp(-Math.Pow(t - T, 2) / (2 * Math.Pow(q, 2)));
        }

        private double FuncLowNorm(double T, double q, double t)
        {
            return 1 - IntegralForNorm(0, t, T, q);
        }

        private double FuncLowShortNorm(double T, double q, double t, double C)
        {
            return C * (1 - IntegralForNorm(0, t, T, q));
        }

 /*       double FuncLowShortNorm(double T, double q, double t, double C)
        {
            double res = C * (0.5 - FLaplas((t - T) / q));
            return res;
        }
*/
        private double IntegralForNorm(double lim1, double lim2, double T, double q)
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

        private double CountDRel(double a)
        {
            /*
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
            for (int i = 1; i <= K; i++)
            {
                Double temp = Math.Exp(-a * (i - 0.5) * dt * (i - 0.5) * dt) - 0.5 * (CountProbability(i - 2) + CountProbability(i - 1));
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

        private double Factorial(int x)
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

        private double FLaplas(double x)
        {
            double res = 1.0 / Math.Sqrt(2 * Math.PI) * IntegralForLaplas(0, x, x);
            return res;
        }

        private double IntegralForLaplas(double lim1, double lim2, double x)
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

        private void Clear()
        {
            var expZedGraph = (ZedGraphControl) ExpWFH.Child;
            var erlZedGraph = (ZedGraphControl) ErlWFH.Child;
            var relZedGraph = (ZedGraphControl) RelWFH.Child;
            var vejbZedGraph = (ZedGraphControl) VejbWFH.Child;
            var normZedGraph = (ZedGraphControl) NormWFH.Child;
            var shortNormZedGraph = (ZedGraphControl) ShortNormWFH.Child;
            var pane = expZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = erlZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = relZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = vejbZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = normZedGraph.GraphPane;
            pane.CurveList.Clear();
            pane = shortNormZedGraph.GraphPane;
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
