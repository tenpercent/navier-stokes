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

**FIXME: write this section**

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
