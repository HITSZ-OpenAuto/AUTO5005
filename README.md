# AUTO5005 - 非线性与自适应控制

![考查课](https://img.shields.io/badge/%E8%80%83%E6%9F%A5%E8%AF%BE-green)
![学分](https://img.shields.io/badge/%E5%AD%A6%E5%88%86-2-moccasin)
![本研共通](https://img.shields.io/badge/本研共通-lightskyblue)

![成绩构成](https://img.shields.io/badge/%E6%88%90%E7%BB%A9%E6%9E%84%E6%88%90-gold)
![平时成绩30%](https://img.shields.io/badge/平时成绩-30%25-wheat)
![期末考试70%](https://img.shields.io/badge/%E6%9C%9F%E6%9C%AB%E8%80%83%E8%AF%95-70%25-wheat)

## 课程基本信息

可见校内网盘 Lec1（Introduction）。

### 教材与参考书

以讲义为主。可参考：[Oliver Wu & Hye's《非线性与自适应控制》笔记](https://oliverwu.top/nac.html)。

{{% details title="参考书（主要是前五本）" closed="true" %}}

1. J-J E. Slotine and Weiping Li, Applied Nonlinear Control, English Edition, 机械工业出版社, 2004
2. Hassan K. Khalil, Nonlinear Systems, Third Eidition, 电子工业出版社, 2007
3. E. Lavertsky and K. Wise, Robust and Adaptive Control with Aerospace Applications, Springer-Verlag, 2013
4. K.J. Åström and B. Wittenmark, Adaptive Control, Second Edition, Dover Publications, INC. Mineola, New York, 2008
5. M. Krstic, I. Kanellakopoulos, and P. Kokotovic, Nonlinear and Adaptive Control Design, John Wiley and Sons, 1995.
6. P. Ioannou and B. Fidan, Adaptive Control  Tutorial, SIAM Press, Philadelphia, 2006
7. Shankar Sastry, Nonlinear Systems: Analysis, Stability, and Control, 世界图书出版社, 1999
8. K. Narendra and A. Annaswamy, Stable Adaptive Systems, Dover Publications, INC. Mineola, New York, 2005.

{{% /details %}}

### 主要内容

（认真看这一节，这可是一道考题）

{{% details title="第一、二讲 绪论" closed="true" %}}
主要内容是： 非线性控制系统概述(什么是非线性系统与自适应控制、为什么要研究它们)、非线性控制系统的建模、非线性常微分方程的解（存在、唯一性，主要是 Lipschitz 条件）。这两讲以了解为主。
{{% /details %}}

{{% details title="第三讲 自治系统的稳定性分析" closed="true" %}}
主要内容是：自治系统 Lyapunov 稳定性的概念（稳定性、渐近稳定、全局渐近稳定、指数稳定）、正定函数、自治系统的 Lyapunov 稳定性定理、LaSalle 不变集原理（分析原点与极限环的稳定性）、线性系统与线性化（利用局部线性化系统的 Jacobian）
{{% /details %}}

{{% details title="第四讲 非自治系统的稳定性分析" closed="true" %}}
主要内容是：比较函数、非自治系统稳定性的概念（除了上述稳定性之外，特别注意新引入的一致稳定）、时变正定函数与非自治系统稳定性定理（注意时变函数正定性的判别）、线性时变系统（了解为主） 、Barbalat 引理及其推论（重点，需掌握证明）、有界性和最终有界性（第六讲用到）、输入-状态稳定性（第六讲用到）
{{% /details %}}

{{% details title="第五讲 自适应控制" closed="true" %}}
主要内容：自适应控制的概念与类型（简答题考查）、各种各样的模型参考自适应控制（线性标量系统[直接型、间接型]、具有非线性项的标量系统、多输入多输出（MIMO）系统）、无参考模型的自适应控制设计、鲁棒自适应控制（死区修正、𝜎-修正、α-修正、自适应 𝜎-修正，需要会推广至多输入多输出（MIMO）系统与时变参数情形）、预设性能控制（2024 新加入）
{{% /details %}}

{{% details title="第六讲 非线性控制系统设计" closed="true" %}}
主要内容是：反馈控制问题的类别（状态反馈镇定问题、输出反馈镇定问题、跟踪问题）、反馈线性化、各种各样的滑模控制（镇定问题、跟踪问题、有不确定性的跟踪问题、有外部干扰）、反步法（基本形式、自适应反步法及减少过参数化、有调节函数（tuning function）的自适应反步法）
{{% /details %}}

{{% details title="第七讲 实例：机械臂控制" closed="true" %}}
主要内容是：将前面讲解过的方法应用到机械臂的控制上，包括基于 Lyapunov 分析的位置控制（主要依靠反馈线性化） 、跟踪控制（使用滑模控制、自适应滑模控制、反步法滑模控制与预设性能控制）
{{% /details %}}

> 文 / [Oliver Wu](https://github.com/OliverWu515), 2025.1

## 授课教师

- 梅杰
  - 上课使用 iPad 手写投屏，节奏飞快。
  - 课后会更新本年度的手写讲义，但并不及时。建议同学们看往年讲义预习。
  - 偶尔以提问方式点名。

## 关于考试

期末考试为闭卷、不允许携带计算器（closed-book, closed-note, closed-calculator exam）。实际上也用不到计算器。
题量很大，不过基本上是上课讲过的例子或者提问的变种（讲义上标有 Q 字样且不给出明确解答的问题）。给分非常宽松（2024 年，98 分排名 7/53），大家即使不太笃定也要把大致思路写上。

## 学习建议
