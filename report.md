<!--- type: markdown; required-extensions: mathjax smarty --->
<!--- Created with ReText: http://retext.sourceforge.net/ --->

Отчёт о численном решении системы уравнений Навье--Стокса.

## Постановка задачи

Дана система уравнений, описывающая одномерное нестационарное движение вязкого
баротропного газа:

$$
\left\{
\begin{align*}
& \frac{\partial\rho}{\partial t} + \frac{\partial\rho u}{\partial x} = 0,\\
& \rho\frac{\partial u}{\partial t} + \rho u\frac{\partial u}{\partial x} +
  \frac{\partial p}{\partial x} = \mu\frac{\partial^2 u}{\partial x^2} + \rho f
\end{align*}
\right.
$$

### Обозначения

Переменная | Описание
---------- | --------
$x$        | пространственная координата
$t$        | время
$\rho$     | плотность газа
$p$        | давление газа, $p = p(\rho)$
$g$        | сокращение для $\ln p$
$u$        | скорость газа
$\mu$      | вязкость газа
$f$        | внешняя сила

### Линейная система

Линейное уравнение для $G$:

$$
- G_{m-1} \cdot \frac{\hat{V}_m + \hat{V}_{m-1}}{4h}
- V_{m-1} \cdot \frac{1}{2h}
+ G_{m}   \cdot \frac{1}{\tau}
+ G_{m+1} \cdot \frac{\hat{V}_{m-1} + \hat{V}_m}{4h}
+ V_{m+1} \cdot \frac{1}{2h}
= \frac{\hat{G}_m}{\tau} +
  \frac{\hat{G}_m \cdot \left( \hat{V}_{m+1} - \hat{V}_{m-1} \right)}{4h} +
  R_1(m)
$$

Начальное условие:

$$
  G_0 \cdot \left( \frac{1}{\tau} - \frac{\hat{V}_0}{2h} \right)
- V_0 \cdot \frac{1}{h}
+ G_1 \cdot \frac{V_1}{2h}
+ V_1 \cdot \frac{1}{h}
= \frac{\hat{G}_0}{\tau} +
  \frac{\hat{G}_0 \cdot \left( \hat{V}_1 - \hat{V}_0 \right)}{2h} +
  \frac{\hat{G}_2 \hat{V}_2 - 2 \hat{G}_1 V_1 + \hat{G}_0 \hat{V}_0 +
      \left( 2 - \hat{G}_0 \right) \left( \hat{V}_2 - 2 \hat{V}_1 + \hat{V}_0 \right)}{4h} +
  R_1(0)
$$

Конечное условие:

$$
- G_{M-1} \cdot \frac{\hat{V}_{M-1}}{2h}
- V_{M-1} \cdot \frac{1}{h}
+ G_{M}   \cdot \left( \frac{1}{\tau} + \frac{\hat{V}_M}{2h} \right)
+ V_{M}   \cdot \frac{1}{h}
= \frac{\hat{G}_M}{\tau} +
  \frac{\hat{G}_M \left( \hat{V}_M - \hat{V}_{M-1} \right) }{2h} +
  \frac{\hat{G}_M \hat{V}_M - 2 \hat{G}_{M-1} \hat{V}_{M-1} + \hat{G}_{M-2} \hat{V}_{M-2} +
      \left( 2 - \hat{G}_M \right) \left( \hat{V}_M - 2 \hat{V}_{M-1} + \hat{V}_{M-2} \right) }{4h} +
  R_1(M)
$$

Линейное уравнение для $V$ (в случае нулевой вязкости):

$$
- G_{m-1} \cdot \frac{\tilde{p}'}{2h}
- V_{m-1} \cdot \left( \frac{\hat{V}_{m-1} + \hat{V}_m}{6h} + \frac{4}{3}\tilde\mu \right)
+ V_{m}   \cdot \left( \frac{1}{\tau} + \frac{8}{3} \tilde\mu \right)
+ G_{m+1} \cdot \frac{\tilde{p}'}{2h}
+ V_{m+1} \cdot \left( \frac{\hat{V}_m + \hat{V}_{m+1}}{6h} - \frac{4}{3}\tilde\mu \right)
= \frac{\hat{V}_m}{\tau} +
  R_2(m)
$$

## Организация хранения матрицы и векторов

Значения $G$ и $V$ располагаются по очереди: $G_0 V_0 G_1 V_1 \dots G_M V_M$.

Поскольку на главной диагонали все элементы с индексами одинаковой чётности
равны между собой, то из диагональной части матрицы хранятся только два элемента.

Остальные элементы хранятся в формате Modified Compressed Sparse Row (MSR),
в виде массива индексов и массива значений.

## Методы решения разреженных линейных систем

В программе реализован метод BiCGSTAB (стабилизированный метод бисопряжённых
градиентов) с предобуславливателем.

Кроме того, возможно подключение библиотеки LASPack, в которой реализованы
другие итерационные методы и предобуславливатели.

### Алгоритм метода BiCGSTAB

Пусть $A$ --- матрица системы, $\mathbf b$ --- правая часть, $\mathbf x$ ---
начальное приближение системы, $K$ --- матрица предобуславливателя.

Начальные значения переменных:

- $\mathbf v = \mathbf p = \mathbf 0$,
- $\mathbf r = \mathbf{\hat r} = \mathbf b - A \mathbf x$.

Каждая итерация состоит из следующих шагов:

1. $\rho_\text{prev} = \rho$,
2. $\rho = (\mathbf{\hat r}, \mathbf r)$,
3. $\beta = (\rho \cdot \alpha) / (\rho_\text{prev} \cdot \omega)$,
4. $\mathbf p = \mathbf r + \beta \cdot (\mathbf p - \omega \mathbf v)$,
5. $\mathbf y = K \mathbf p$,
6. $\mathbf v = A \mathbf y$,
7. $\alpha = \rho / (\mathbf{\hat r}, \mathbf v)$,
8. $\mathbf s = \mathbf r - \alpha \mathbf v$,
9. $\mathbf z = K \mathbf s$,
10. $\mathbf t = A \mathbf z$,
11. $\omega = (\mathbf t, \mathbf s) / (\mathbf t, \mathbf t)$,
12. $\mathbf x = \mathbf x + \alpha \mathbf y + \omega \mathbf z$,
13. $d = \left\|A \mathbf x - \mathbf b\right\|$,
14. $\mathbf r = \mathbf s - \omega \mathbf t$.

Алгоритм заканчивается, когда $d$ будет меньше требуемой точности, или будет
достигнуто максимальное число итераций (1000).

В качестве предобуславливателя используется предобуславливатель Якоби
с матрицей $K = \omega_\text{prec} \cdot D^{-1}(A)$.

## Проведённые тесты

Были проведены испытания на сетке размера 400×400 при разных значениях $\tau$, $h$.
Результаты представлены в таблице:

Время    | $\omega$ | $\tau$      | $h$         | Невязка $V$ в $C$ | Невязка $V$ в $L_2$ | Невязка $G$ в $C$ | Невязка $G$ в $L_2$
-------- | -------- | ----------- | ----------- | ----------------- | ------------------- | ----------------- | -------------------
1.682 с  | 1.0      | 2.50627e-03 | 2.50627e-02 | 2.08565e-03       | 1.02849e-03         | 1.43495e-03       | 5.39350e-04
3.383 с  | 1.0      | 1.25156e-03 | 2.50627e-02 | 1.10104e-03       | 5.32001e-04         | 7.14881e-04       | 2.69397e-04
6.440 с  | 1.0      | 6.25391e-04 | 2.50627e-02 | 6.01320e-04       | 2.80338e-04         | 3.56732e-04       | 1.35181e-04
4.473 с  | 1.0      | 2.50627e-03 | 1.25156e-02 | 2.11002e-03       | 1.04097e-03         | 1.45481e-03       | 5.44843e-04
6.986 с  | 1.0      | 1.25156e-03 | 1.25156e-02 | 1.12532e-03       | 5.44123e-04         | 7.34616e-04       | 2.75298e-04
12.325 с | 1.0      | 6.25391e-04 | 1.25156e-02 | 6.25383e-04       | 2.92187e-04         | 3.76392e-04       | 1.41248e-04
10.467 с | 1.0      | 2.50627e-03 | 6.25391e-03 | 2.11612e-03       | 1.04425e-03         | 1.45975e-03       | 5.45978e-04
16.058 с | 1.0      | 1.25156e-03 | 6.25391e-03 | 1.13137e-03       | 5.47243e-04         | 7.39548e-04       | 2.76651e-04
25.220 с | 1.0      | 6.25391e-04 | 6.25391e-03 | 6.31391e-04       | 2.95209e-04         | 3.81302e-04       | 1.42702e-04
