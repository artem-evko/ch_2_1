package org.artem.evko;

import java.io.*;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) {

        String inputFileName = "input.txt";
        String outputFileName = "output.txt";

        double[][] A = null;
        int N = 0;
        try (Scanner sc = new Scanner(new File(inputFileName))) {
            sc.useLocale(Locale.US);
            N = sc.nextInt();
            A = new double[N][N];
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    A[i][j] = sc.nextDouble();
                }
            }
        } catch (FileNotFoundException e) {
            System.err.println("Файл не найден: " + inputFileName);
            return;
        }

        double[] x = new double[N];
        Arrays.fill(x, 1.0);


        double eps = 1e-12;
        int maxIter = 1000;


        double lambdaOld = 0.0;
        double lambdaCurrent = 0.0;
        double lambdaOlder = 0.0;


        int iterCount = 0;
        boolean converged = false;


        double[] Ax = multiply(A, x);
        double normAx = norm(Ax);
        if (normAx < 1e-30) {
            System.err.println("Начальный вектор дает нулевое преобразование A*x");
            return;
        }

        for (int i = 0; i < N; i++) {
            x[i] = Ax[i] / normAx;
        }

        Ax = multiply(A, x);
        lambdaCurrent = dot(x, Ax);  // (x^T)(A x)


        while (iterCount < maxIter) {
            iterCount++;


            lambdaOlder = lambdaOld;
            lambdaOld = lambdaCurrent;


            Ax = multiply(A, x);
            normAx = norm(Ax);


            if (normAx < 1e-30) {
                System.err.println("Произошло вырождение на итерации: " + iterCount);
                break;
            }

            for (int i = 0; i < N; i++) {
                x[i] = Ax[i] / normAx;
            }


            Ax = multiply(A, x);
            lambdaCurrent = dot(x, Ax);


            if (Math.abs(lambdaCurrent - lambdaOld) < eps) {
                converged = true;
                break;
            }
        }


        double lambdaPower = lambdaCurrent;


        double lambdaAitken = lambdaPower;
        if (iterCount >= 2) {
            double d1 = (lambdaOld - lambdaOlder);
            double d2 = (lambdaCurrent - 2.0*lambdaOld + lambdaOlder);
            if (Math.abs(d2) > 1e-30) {

                lambdaAitken = lambdaOlder - (d1 * d1)/d2;
            }
        }


        try (PrintWriter out = new PrintWriter(new FileWriter(outputFileName))) {
            out.println("Число итераций: " + iterCount + (converged ? " (успешно)" : " (достигнут предел итераций)"));
            out.printf(Locale.US, "Приближение собств. значения (степенной метод): %.12f%n", lambdaPower);
            out.printf(Locale.US, "Приближение собств. значения (Эйткен):          %.12f%n", lambdaAitken);
            out.println("Нормированный собственный вектор (приближение):");
            for (int i = 0; i < N; i++) {
                out.printf(Locale.US, "%.12f ", x[i]);
            }
            out.println();
        } catch (IOException e) {
            System.err.println("Ошибка записи в файл: " + outputFileName);
        }
    }


    private static double[] multiply(double[][] A, double[] x) {
        int N = x.length;
        double[] result = new double[N];
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            result[i] = sum;
        }
        return result;
    }


    private static double dot(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }


    private static double norm(double[] v) {
        double sum = 0.0;
        for (double val : v) {
            sum += val * val;
        }
        return Math.sqrt(sum);
    }
}
