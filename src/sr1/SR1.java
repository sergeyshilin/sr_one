package sr1;

/**
 * Created with IntelliJ IDEA.
 * User: Sergey Shilin
 * Mail: sergey.shilin@phystech.edu
 * Date: 16.05.13
 * Time: 15:49
 */
public class SR1 {
    public static void main(String[] args) {
        final SR1Point startPoint = new SR1Point(6, 4); // задаем начальное приближение

        SR1Function f = new SR1Function() {
            @Override
            public double calculate(SR1Point point) {
                /**
                 * Указываем, какие переменные будут присутствовать в функции
                 */
                double x = point.get(0);
                double y = point.get(1);

                return Math.exp(Math.pow((x-5), 2) + Math.pow((y-3), 2));  // указываем вид функции(тут все, что душе угодно)
            }
        };
        new SR1Optimizer().optimize(f, startPoint);  // Запускаем алгоритм оптимизации +MSR1
    }
}
