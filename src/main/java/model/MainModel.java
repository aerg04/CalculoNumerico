package model;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author DELL & Doztor
 */
public class MainModel {
    
        // Parámetros del circuito
    private static final double R_eq = 3 + 1 + 4; // Resistencia equivalente
    private static final double L_eq = 10e-3 + 10e-3; // Inductancia equivalente
    private static final double C = 100e-6; // Capacitancia
    
    //PARTE 1

    // Nuevo método para calcular la corriente i(t)
    public XYSeriesCollection createDatasetCurrent(double t0, double i0, double dt, double tEnd) {
        double[] y = new double[]{i0};
        double t = t0;

        XYSeries series = new XYSeries("Corriente i(t)");

        while (t <= tEnd) {
            series.add(t, y[0]);
            y = rungeKutta4thOrderCurrent(t, y, dt);
            t += dt;
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        return dataset;
    }

    // Nuevo método para calcular el voltaje V(t)
    public XYSeriesCollection createDatasetVoltage(double t0, double i0, double dt, double tEnd) {
        double[] y = new double[]{i0};
        double t = t0;

        XYSeries seriesCurrent = new XYSeries("Corriente i(t)");
        XYSeries seriesVoltage = new XYSeries("Voltaje V(t)");

        while (t <= tEnd) {
            seriesCurrent.add(t, y[0]);

            // Calcular voltaje en el circuito: V(t) = V_R + V_L + V_C
            double V_R = R_eq * y[0];
            double V_L = L_eq * dydtCurrent(t, y)[0];
            double V_C = (1 / C) * integrateCurrent(seriesCurrent, t0, t);
            seriesVoltage.add(t, V_R + V_L + V_C);

            y = rungeKutta4thOrderCurrent(t, y, dt);
            t += dt;
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(seriesVoltage);
        return dataset;
    }

 
    // Método Runge-Kutta para la corriente i(t)
    private double[] rungeKutta4thOrderCurrent(double t, double[] y, double dt) {
        double[] k1 = multiply(dydtCurrent(t, y), dt);
        double[] k2 = multiply(dydtCurrent(t + dt / 2, add(y, divide(k1, 2))), dt);
        double[] k3 = multiply(dydtCurrent(t + dt / 2, add(y, divide(k2, 2))), dt);
        double[] k4 = multiply(dydtCurrent(t + dt, add(y, k3)), dt);

        return add(y, divide(add(add(k1, multiply(add(k2, k3), 2)), k4), 6));
    }

    // Ecuación diferencial para la corriente i(t)
    private double[] dydtCurrent(double t, double[] y) {
        double i = y[0];
        double diDt = (120 * 60 * Math.cos(60 * t) - 1.001803e-4 * i) / 2;
        return new double[]{diDt};
    }

    // Método para integrar la corriente i(t) (usado para V_C)
    private double integrateCurrent(XYSeries series, double t0, double t) {
        double integral = 0;
        for (int i = 1; i < series.getItemCount(); i++) {
            double dt = series.getX(i).doubleValue() - series.getX(i - 1).doubleValue();
            double avgCurrent = (series.getY(i).doubleValue() + series.getY(i - 1).doubleValue()) / 2;
            integral += avgCurrent * dt;
        }
        return integral;
    }


//PARTE 2
    public XYSeriesCollection createDatasetPendulum(double lambdax,double omegax,double t0,double theta) {
        double lambda = lambdax;
        double omega = omegax;
        double h = 0.01;
        double[] y0 = {theta, t0};
        double a = 0;
        double b = 10.0;

        List<Double> tValues = new ArrayList<>();
        List<double[]> yValues = new ArrayList<>();
        rungeKutta4thOrder(tValues, yValues, a, b, y0, h, lambda, omega);

        XYSeries series = new XYSeries("Pendulum");

        for (int i = 0; i < tValues.size(); i++) {
            series.add((Number) tValues.get(i), (Number) yValues.get(i)[0]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        return dataset;
    }

    private void rungeKutta4thOrder(List<Double> tValues, List<double[]> yValues,
                                    double a, double b, double[] y0, double h, double lambda, double omega) {
        int N = (int) ((b - a) / h);
        double t = a;
        double[] w = y0.clone();
        tValues.add(t);
        yValues.add(w.clone());

        for (int i = 1; i <= N; i++) {
            double[] K1 = multiply(pendulum(t, w, lambda, omega), h);
            double[] K2 = multiply(pendulum(t + h / 2, add(w, divide(K1, 2)), lambda, omega), h);
            double[] K3 = multiply(pendulum(t + h / 2, add(w, divide(K2, 2)), lambda, omega), h);
            double[] K4 = multiply(pendulum(t + h, add(w, K3), lambda, omega), h);

            w = add(w, divide(add(add(K1, multiply(add(K2, K3), 2)), K4), 6));
            t = a + i * h;

            tValues.add(t);
            yValues.add(w.clone());
        }
    }

    private double[] pendulum(double t, double[] y, double lambda, double omega) {
        double theta = y[0];
        double thetaDot = y[1];
        double dThetaDt = thetaDot;
        double dThetaDotDt = -2 * lambda * thetaDot - omega * omega * Math.sin(theta);
        return new double[]{dThetaDt, dThetaDotDt};
    }

    private double[] add(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    private double[] multiply(double[] a, double b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] * b;
        }
        return result;
    }

    private double[] divide(double[] a, double b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] / b;
        }
        return result;
    }
}
