using System.Windows;

namespace NumericsLabs
{
    public partial class MainWindow : Window
    {
        private Controller controller;
        public MainWindow()
        {
            InitializeComponent();

            controller = new Controller();

            CalcComboBox.ItemsSource = controller.CalculationNames;
            CalcComboBox.SelectedIndex = 0;

            RefreshIssue();

            TimeSlider.Minimum = 0;
            TimeSlider.Maximum = controller.Issue.Parameters.K * controller.Issue.Parameters.Tau;
            TimeSlider.TickFrequency = controller.Issue.Parameters.Tau;
            TimeSlider.Value = 0;
            TimeSlider.IsSnapToTickEnabled = true;

            RefreshErrorGraph();
            RefreshGraph();
        }

        private void RefreshIssue()
        {
            if (CalcComboBox.SelectedItem != null)
            {
                controller.SetCalculator(CalcComboBox.SelectedItem.ToString());
                controller.Refresh();
            }
        }

        private void RefreshGraph()
        {
            GraphVisualizer.Plot.Clear();

            int t = (int)(TimeSlider.Value / controller.Issue.Parameters.Tau);
            GraphVisualizer.Plot.AddScatter(controller.RangeOfX, controller.GetResultByTIndex(t), System.Drawing.Color.Red, 2, 2);
            GraphVisualizer.Plot.AddScatter(controller.RangeOfX, controller.GetAnalyticalResultByTIndex(t), System.Drawing.Color.Blue, 3, 3);
            GraphVisualizer.Plot.XLabel("x");
            GraphVisualizer.Plot.YLabel($"u(x, t={t})");

            GraphVisualizer.Refresh();
        }
        private void RefreshErrorGraph()
        {
            ErrorVisualizer.Plot.Clear();
            ErrorVisualizer.Plot.AddScatter(controller.RangeOfT, controller.GetError());
            ErrorVisualizer.Plot.XLabel("t");
            ErrorVisualizer.Plot.YLabel("ε(t)");
            ErrorVisualizer.Refresh();
        }

        private void Slider_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            RefreshGraph();
        }

        private void SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            if (controller?.Issue != null)
            {
                RefreshIssue();
                RefreshErrorGraph();
                RefreshGraph();
            }
        }
    }
}
