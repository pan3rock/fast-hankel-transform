# fast-hankel-transform

## Fast Hankel transform through FFT<sup>[1](#footnote1),[2](#footnote2)</sup> 

定义0阶Hankel变换：

<a href="https://www.codecogs.com/eqnedit.php?latex=H(y)&space;=&space;\int_0^1&space;f(x)J_0(C\cdot&space;yx)xdx&space;\quad&space;(0&space;\le&space;y&space;\le&space;1)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?H(y)&space;=&space;\int_0^1&space;f(x)J_0(C\cdot&space;yx)xdx&space;\quad&space;(0&space;\le&space;y&space;\le&space;1)" title="H(y) = \int_0^1 f(x)J_0(C\cdot yx)xdx \quad (0 \le y \le 1)" /></a>

则在特定采样点上可以用FFT计算

<a href="https://www.codecogs.com/eqnedit.php?latex=H(y_n)&space;=&space;\frac{1}{C&space;y_n}&space;F\{F\{\phi_n\}\cdot&space;F^{-1}\{J_{0n}\}\}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?H(y_n)&space;=&space;\frac{1}{C&space;y_n}&space;F\{F\{\phi_n\}\cdot&space;F^{-1}\{J_{0n}\}\}" title="H(y_n) = \frac{1}{C y_n} F\{F\{\phi_n\}\cdot F^{-1}\{J_{0n}\}\}" /></a>

其中特定采样点取值如下

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;y_n&space;&=&space;x_n&space;=&space;x_0&space;exp(\alpha&space;n)&space;\qquad&space;\textrm{for}&space;\quad&space;n&space;=&space;0,&space;1,&space;\dots,&space;N-1\\&space;x_0&space;&=&space;[1&space;&plus;&space;exp(\alpha)]exp(-\alpha&space;N)&space;/&space;2&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;y_n&space;&=&space;x_n&space;=&space;x_0&space;exp(\alpha&space;n)&space;\qquad&space;\textrm{for}&space;\quad&space;n&space;=&space;0,&space;1,&space;\dots,&space;N-1\\&space;x_0&space;&=&space;[1&space;&plus;&space;exp(\alpha)]exp(-\alpha&space;N)&space;/&space;2&space;\end{align*}" title="\begin{align*} y_n &= x_n = x_0 exp(\alpha n) \qquad \textrm{for} \quad n = 0, 1, \dots, N-1\\ x_0 &= [1 + exp(\alpha)]exp(-\alpha N) / 2 \end{align*}" /></a>

那么，如果积分区间在[0, a]，而变量y的范围为[0, b]，对应的Hankel变换为

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;\hat{H}(y')&space;=&&space;\hat{H}(bx[n])\\&space;=&\int_0^{a}f(x')J_0(y'x')x'dx'\\&space;=&a^2\int_0^1&space;f(ax)&space;J_0(y'&space;ax)xdx\\&space;=&a^2\int_0^1&space;f(ax)&space;J_0(abyx)xdx\\&space;=&a^2H(f(x'),&space;C=ab,&space;y&space;=&space;x[n])&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\begin{align*}&space;\hat{H}(y')&space;=&&space;\hat{H}(bx[n])\\&space;=&\int_0^{a}f(x')J_0(y'x')x'dx'\\&space;=&a^2\int_0^1&space;f(ax)&space;J_0(y'&space;ax)xdx\\&space;=&a^2\int_0^1&space;f(ax)&space;J_0(abyx)xdx\\&space;=&a^2H(f(x'),&space;C=ab,&space;y&space;=&space;x[n])&space;\end{align*}" title="\begin{align*} \hat{H}(y') =& \hat{H}(bx[n])\\ =&\int_0^{a}f(x')J_0(y'x')x'dx'\\ =&a^2\int_0^1 f(ax) J_0(y' ax)xdx\\ =&a^2\int_0^1 f(ax) J_0(abyx)xdx\\ =&a^2H(f(x'), C=ab, y = x[n]) \end{align*}" /></a>

## ~~Less points to evaluate Hankel transform<sup>[3](#footnote3)</sup>~~

## References

<a name="footnote1">1</a>. Zhang, D., Yuan, X., Ngo, N., & Shum, P. (2002). Fast Hankel transform and its application for studying the propagation of cylindrical electromagnetic fields. Optics Express, 10(12), 521. https://doi.org/10.1364/OE.10.000521

<a name="footnote2">2</a>. Magni, V., Cerullo, G., & De Silvestri, S. (1992). High-accuracy fast Hankel transform for optical beam propagation. Journal of the Optical Society of America A, 9(11), 2031. https://doi.org/10.1364/JOSAA.9.002031

<a name="footnote3">3</a>. Levin, D. (1996). Fast integration of rapidly oscillatory functions. Journal of Computational and Applied Mathematics, 67(1), 95–101. https://doi.org/10.1016/0377-0427(94)00118-9