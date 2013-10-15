package sr1;

/**
 * Created with IntelliJ IDEA.
 * User: Sergey Shilin
 * Mail: sergey.shilin@phystech.edu
 * Date: 18.05.13
 * Time: 10:49
 */
public class SR1MathUtils {

    public static final double ARGUMENT_DELTA = 1e-6;
    public static final double EPSILON = 1e-3;
    public static final double DELTA_ONE = 1e-4;
    public static final double DELTA_TWO = 0.9;
    public static final double MAX_ALPHA = 1.0;
    public static final double MIN_ALPHA = 0.0;


    /**
     * Подсчет частной производной функции в точке.
     */
    public static SR1Function derivative(final SR1Function f, final int position) {
        return new SR1Function() {
            @Override
            public double calculate(SR1Point spoint) {
                SR1Point derivPoint = new SR1Point(spoint);
                derivPoint.set(position, spoint.get(position) + ARGUMENT_DELTA);
                return (f.calculate(derivPoint) - f.calculate(spoint)) / (ARGUMENT_DELTA);
            }
        };
    }

    /**
     * Подсчет градиента функции в точке
     */
    public static SR1Point grad(SR1Function f, SR1Point point) {
        SR1Point newpoint = new SR1Point(point.size());
        for (int i = 0; i < point.size(); i++) {
            SR1Function deriv = derivative(f, i);
            newpoint.set(i, deriv.calculate(point));
        }
        return newpoint;
    }

    /**
     * Подсчет нормы вектора(точки) как корень из суммы квадратов
     */
    public static double norm(SR1Point point) {
        return point.norm2();
    }
}
