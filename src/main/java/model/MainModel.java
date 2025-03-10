package model;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.UnivariateFunction;


import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author DELL & Doztor
 */
public class MainModel {
    
    //PARTE 1
    
     // Método para calcular la corriente (RK4)
     public XYSeriesCollection createDatasetCorriente() {
        double t0 = 0;
        double i0 = 0;
        double dt = 0.01;
        double t_end = 1;

        XYSeries series = new XYSeries("Corriente");

        double t = t0;
        double i = i0;

        while (t <= t_end) {
            double[] rkResult = rungeKuttaCorriente(t, i, dt);
            t = rkResult[0];
            i = rkResult[1];
            series.add(t, i);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        return dataset;
    }

    private double[] rungeKuttaCorriente(double t, double i, double dt) {
        double k1 = dt * dydtCorriente(t, i);
        double k2 = dt * dydtCorriente(t + dt / 2.0, i + k1 / 2.0);
        double k3 = dt * dydtCorriente(t + dt / 2.0, i + k2 / 2.0);
        double k4 = dt * dydtCorriente(t + dt, i + k3);
        double i_new = i + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        double t_new = t + dt;
        return new double[]{t_new, i_new};
    }

    private double dydtCorriente(double t, double i) {
        return (120 * 60 * Math.cos(60 * t) - (1.001803e-4) * i);
    }

    // Método para calcular el voltaje
    public XYSeriesCollection createDatasetVoltaje() {
        double t0 = 0;
        double i0 = 0;
        double dt = 0.001;
        double t_end = 2;

        XYSeries series = new XYSeries("Voltaje");

        double t = t0;
        double i = i0;

        while (t <= t_end) {
            double[] rkResult = rungeKuttaCorriente(t, i, dt);
            t = rkResult[0];
            i = rkResult[1];
            double V = calcularVoltaje(t, i);
            series.add(t, V);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        return dataset;
    }

    private double calcularVoltaje(double t, double i) {
        double R_eq = 3 + 1 + 4; // Resistencia equivalente
        double L_eq = 10e-3 + 10e-3; // Inductancia equivalente
        double C = 100e-6; // Capacitancia

        double V_R = R_eq * i; // Voltaje en la resistencia
        double V_L = L_eq * dydtCorriente(t, i); // Voltaje en la inductancia
        double V_C = (1.0 / C) * integrarCorriente(t, i); // Voltaje en el capacitor
        return V_R + V_L + V_C;
    }

    private double integrarCorriente(double t, double i) {
        // Crear una instancia de TrapezoidIntegrator
        TrapezoidIntegrator integrator = new TrapezoidIntegrator();

        // Definir la función a integrar (en este caso, una función constante: f(x) = i)
        UnivariateFunction function = x -> i;

        // Integrar la función en el intervalo [0, t]
        return integrator.integrate(1000, function, 0, t);
    }


//PARTE 2
    public XYSeriesCollection createDatasetPendulum(double lambdax,double omegax,double t0,double theta) {
        double lambda = lambdax;
        double omega = omegax;
        double h = 0.01;
        double[] y0 = {theta, t0};
        double a = t0;
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
