package sr1;

import Jama.Matrix;

/**
 * Created with IntelliJ IDEA.
 * User: Sergey Shilin
 * Mail: sergey.shilin@phystech.edu
 * Date: 17.05.13
 * Time: 19:57
 */
public class SR1Point extends Matrix {
    /**
     * Класс используется для описания точек приближения, градиента
     * и результата произведения матрицы на градиент - вектора.
     */

    public SR1Point(double... args) {
        super(args, args.length);
    }

    public SR1Point(int size) {
        super(size, 1);
    }

    public SR1Point(SR1Point point) {
        this(point.getPointArray());
    }

    public SR1Point(Matrix matrix) {
        super(matrix.getArray());
    }

    public Double get(int i) {
        return (Double) get(i, 0);
    }

    public int size() {
        return getRowDimension();
    }

    public void set(int position, double arg) {
        set(position, 0, arg);
    }

    public void print() {
        System.out.println();
        for (int i = 0; i < getRowDimension(); i++)
            System.out.print("  " + get(i));
        System.out.println();
    }

    public double[] getPointArray() {
        double[] args = new double[getRowDimension()];
        for (int i = 0; i < getRowDimension(); i++)
            args[i] = get(i);
        return args;
    }
}
