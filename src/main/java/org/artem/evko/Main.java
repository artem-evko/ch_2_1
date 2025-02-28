
import java.io.*;
import java.util.Arrays;

public class AitkenEigenvalues {

    private static final double EPS = 1e-6; // Точность вычислений

    // Метод для чтения матрицы из файла
    private static double[][] readMatrixFromFile(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = reader.readLine();
        int n = Integer.parseInt(line);
        double[][] matrix = new double[n][n];
        
        for (int i = 0; i < n; i++) {
            line = reader.readLine();
            String[] values = line.split(" ");
            for (int j = 0; j < n; j++) {
                matrix[i][j] = Double.parseDouble(values[j]);
            }
        }
        reader.close();
        return matrix;
    }

    // Метод для записи результата в файл
    private static void writeResultToFile(String filename, double eigenvalue, double[] eigenvector) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
        writer.write("Собственное значение: " + eigenvalue + "\n");
        writer.write("Собственный вектор:\n");
        for (double value : eigenvector) {
            writer.write(value + " ");
        }
        writer.close();
    }

    // Вычисление скалярного произведения векторов
    private static double dotProduct(double[] u, double[] v) {
        double result = 0;
        for (int i = 0; i < u.length; i++) {
            result += u[i] * v[i];
        }
        return result;
    }

    // Умножение матрицы на вектор
    private static double[] matrixVectorMultiply(double[][] A, double[] x) {
        int n = A.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i] += A[i][j] * x[j];
            }
        }
        return result;
    }

    // Нормализация вектора
    private static double[] normalize(double[] v) {
        double norm = 0;
        for (double value : v) {
            norm += value * value;
        }
        norm = Math.sqrt(norm);
        
        double[] result = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            result[i] = v[i] / norm;
        }
        return result;
    }

    // Процесс Эйткена для ускорения сходимости
    private static double aitkenAcceleration(double x0, double x1, double x2) {
        return x2 - ((x2 - x1) * (x2 - x1)) / (x2 - 2 * x1 + x0);
    }

    // Основной метод нахождения собственного значения
    public static void findEigenvalue(String inputFile, String outputFile) throws IOException {
        double[][] A = readMatrixFromFile(inputFile);
        int n = A.length;
        
        // Начальное приближение
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = 1;
        }
        
        double lambdaOld = 0;
        double lambdaNew = 0;
        double[] xOld = new double[n];
        double[] xNew = new double[n];
        
        // Итерационный процесс
        while (true) {
            // Сохраняем предыдущее приближение
            System.arraycopy(x, 0, xOld, 0, n);
            lambdaOld = lambdaNew;
            
            // Умножаем матрицу на вектор
            xNew = matrixVectorMultiply(A, x);
            
            // Вычисляем приближение собственного значения
            lambdaNew = dotProduct(x, xNew) / dotProduct(x, x);
            
            // Нормализуем вектор
            x = normalize(xNew);
            
            // Проверяем сходимость
            if (Math.abs(lambdaNew - lambdaOld) < EPS) {
                break;
