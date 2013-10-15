package sr1;

import Jama.Matrix;
import sr1.SR1MathUtils;

import java.util.Random;

import static sr1.SR1MathUtils.grad;
import static sr1.SR1MathUtils.norm;

/**
 * Created with IntelliJ IDEA.
 * User: Sergey Shilin
 * Mail: sergey.shilin@phystech.edu
 * Date: 18.05.13
 * Time: 10:48
 */
public class SR1Optimizer {

    private Matrix hessian;
    private int hessianDimention;
    private double alpha;
    private Matrix direction;

    public void optimize(SR1Function f, SR1Point startPoint) {
        setHessianDimention(startPoint.size()); // Задаем размерность матрицы Гессе
        setHessian();  // Изначально матрица Гессе есть единичная матрица

        SR1Point curpoint = new SR1Point(startPoint);  // Определяем текущую точку приближения
        while(!isStopCondition(f, curpoint)) {
            /**
             * Оптимизируем, пока не выполняются условия остановки программы.
             * 1) Получаем информацию о кривизне функци (Находим направление)
             */
            direction = (hessian.times(grad(f, curpoint))).uminus();
            SR1Point next = getNextPoint(f, curpoint);  // 2) получаем следующую точку (приближение)

            /**
             * 3)
             *     <-- Далее задаем необходимые для данного шага алгоритма переменные -->
             * Sk - разность между следующим и предыдущим приближением
             * Yk - разность градиентов функции для следующего и тукущего приближения
             * Psi - некоторое double число, вычисляется в шаге 7 pdf-файла
             * Ycap - вектор(точка), определяется в том же шаге
             * lambda - double число, вычисляется по ф-ле (11)
             */

            SR1Point sk = new SR1Point(next.minus(curpoint));
            SR1Point yk = new SR1Point(grad(f, next).minus(grad(f, curpoint)));
            double psi = getPsiFunction(f, curpoint, next, sk);
            SR1Point ycap = new SR1Point(getYCap(yk, psi, sk));
            double lambda = getLambda(sk, ycap);

            /**
             * 4) Пытаемся обновить матрицу Гессе, основываясь на условиях (3* - 5*)
             */
            tryUpdateHessian(f, curpoint, sk, yk, lambda);
            setNextHessian(sk, ycap);   // 5) Определяем матрицу Гессе для следующего шага. (формула 12)

            curpoint = next;
            System.out.println("current point is");
            curpoint.print();

        }

        System.out.print("final point is");
        curpoint.print();
        System.out.println("function value in this point is " + f.calculate(curpoint));
        System.out.print("Hessian matrix looks so: ");
        hessian.print(hessianDimention, hessianDimention);

    }

    private void setNextHessian(SR1Point sk, SR1Point ycap) {
        SR1Point a = new SR1Point(sk.minus(hessian.times(ycap)));
        double denominator = Math.pow(a.transpose().times(ycap).get(0, 0), -1);
        hessian = hessian.plus((a.times(a.transpose())).times(denominator));
    }

    private double getLambda(SR1Point sk, SR1Point ycap) {
        double stsk = sk.transpose().times(sk).get(0, 0);
        double ytsk = ycap.transpose().times(sk).get(0, 0);
        double ytyk = ycap.transpose().times(ycap).get(0, 0);
        double undersqrt = Math.pow((stsk / ytsk), 2) - (stsk / ytyk);
        double sqrt = (undersqrt > 0) ? Math.sqrt(undersqrt) : 0;
        /**
         * Вот тут корень часто хочет считаться из отрицательного значения, но
         * так как это значение достаточно мало, возвращаем ноль.
         */
        return ((stsk / ytsk) - sqrt);
    }

    private double getPsiFunction(SR1Function f, SR1Point curpoint, SR1Point next, SR1Point sk) {
        double first = 2 * (f.calculate(curpoint) - f.calculate(next));
        double second = (grad(f, next).plus(grad(f, curpoint))).transpose().times(sk).get(0, 0);
        return first + second;
    }

    private SR1Point getYCap(SR1Point yk, double psi, SR1Point sk) {
        double skfactor = Math.signum(psi) * ( psi / sk.transpose().times(sk).get(0, 0) );
        return new SR1Point(yk.plus(sk.times(skfactor)));
    }

    private void tryUpdateHessian(SR1Function f, SR1Point curpoint, SR1Point sk, SR1Point yk, double lambda) {
        if(canUpdateHessian(sk, yk)) {
            hessian = Matrix.identity(hessianDimention, hessianDimention).times(lambda);
            direction = grad(f, curpoint).times(-1 * lambda);
        }
    }

    private boolean canUpdateHessian(SR1Point sk, SR1Point yk) {
        boolean firstcond = sk.transpose().times(yk).get(0, 0) - (yk.transpose().times(hessian)).times(yk).get(0, 0) < 0;
//        boolean secondcond = ;
//        boolean thirdcond = ;
//
//        return firstcond || secondcond || thirdcond;
        /**
         * Нужно рассмотреть два других условия. Достаточно одного из них,
         * но для большей точности обычно рассматривают и другие
         */
        return firstcond;
    }

    private SR1Point getNextPoint(SR1Function f, SR1Point point) {
        setAlpha(1.0); // Изначально пробуем длину шага alpha = 1;
        while(!wolfeConditions(f, point)) {
            /**
             * После первого круга пробуем взять рандомное alpha
             * в границах [MIN_ALPHA, MAX_ALPHA]
             */
            Random r = new Random();
            alpha = SR1MathUtils.MIN_ALPHA + (SR1MathUtils.MAX_ALPHA - SR1MathUtils.MIN_ALPHA) * r.nextDouble();
        }
        SR1Point nextpoint = new SR1Point(point.plus(direction.times(alpha)));
        return nextpoint;
    }

    private boolean wolfeConditions(SR1Function f, SR1Point point) {
        /**
         *  Запускаем рандомизацию для alpha, чтобы найти такую длину шага, которая
         *  будет удовлетворять условиям Wolfe.
         */
        System.out.println("alpha is " + alpha);
        SR1Point nextpoint = new SR1Point(point.plus(direction.times(alpha))); // следующая точка (с учетом текущего шага)
        double matrixpoint = (grad(f, point).transpose().times(direction)).get(0, 0);

        /**
         * matrixpoint - результат перемножения матриц размерностей (1, 3)x(3, 1) к примеру.
         * Задаем левую и правую части первого условия Wolfe
         */
        double firstcondleft = f.calculate(nextpoint);
        double firstcondright = f.calculate(point) + SR1MathUtils.DELTA_ONE * alpha * matrixpoint;
        /**
         * Далее задаем левую и правую части второго условия Wolfe
         */
        double secondcondleft = (grad(f, nextpoint).transpose().times(direction)).get(0, 0);
        double secondcondright = SR1MathUtils.DELTA_TWO * matrixpoint;
        /**
         * Описываем сами условия
         */
        boolean firstcond = firstcondleft <= firstcondright;
        boolean secondcond = secondcondleft >= secondcondright;
        return firstcond && secondcond;
    }

    private boolean isStopCondition(SR1Function f, SR1Point point) {
        double leftcond = norm(grad(f, point));
        double rightcond = SR1MathUtils.EPSILON * Math.max(1, norm(point));
        return leftcond <= rightcond;
    }

    private void setHessianDimention(int dimention) {
        hessianDimention = dimention;
    }

    private void setHessian() {
        hessian = Matrix.identity(hessianDimention, hessianDimention);
    }


    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }
}
