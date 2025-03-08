/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package model;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author DELL
 */
public class MainModel {
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
