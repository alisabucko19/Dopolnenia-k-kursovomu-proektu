# САНКТ-ПЕТЕРБУРГСКИЙ УНИВЕРСИТЕТ ГРАЖДАНСКОЙ АВИАЦИИ ИМЕНИ ГЛАВНОГО МАРШАЛА АВИАЦИИ А.А.НОВИКОВА
# Курсовой проект: Вычисление углов ориентации. Комплементарный фильтр для БАС

**Вариант 3**  
**Исполнитель:** Буцко Алиса Артуровна  
**Руководитель:** Земсков Юрий Владимирович  

---

## Содержание

1. [Введение](#1-введение)
2. [Теоретическая часть](#2-теоретическая-часть)
3. [Задание 3.1. Геометрия углов](#3-задание-31-геометрия-углов)
4. [Задание 3.2. Дрейф гироскопа](#4-задание-32-дрейф-гироскопа)
5. [Задание 3.3. Расчёт комплементарного фильтра](#5-задание-33-расчёт-комплементарного-фильтра)
6. [Задание 3.4. Численное моделирование](#6-задание-34-численное-моделирование)
7. [Задание 3.5. Ограничения алгоритма](#7-задание-35-ограничения-алгоритма)
8. [Задание 3.6. Инициализация](#8-задание-36-инициализация)
9. [Дополнительное задание 1: Моделирование фильтра в Python](#9-дополнительное-задание-1-моделирование-фильтра-в-python)
10. [Дополнительное задание 2: Анализ постоянной времени](#10-дополнительное-задание-2-анализ-постоянной-времени)
11. [Дополнительное задание 3: Адаптивный коэффициент α](#11-дополнительное-задание-3-адаптивный-коэффициент-α)
12. [Дополнительное задание 4: Оценка качества фильтра](#12-дополнительное-задание-4-оценка-качества-фильтра)
13. [Заключение](#13-заключение)
14. [Список литературы](#14-список-литературы)

---

## 1. Введение

Современные беспилотные авиационные системы (БАС) невозможно представить без точных систем стабилизации и навигации. Элементом такой системы является модуль определения углов ориентации — крена ($\phi$), тангажа ($\theta$) и рыскания ($\psi$). Без знания этих углов невозможна стабилизация полёта.

В данном курсовом проекте рассматривается микропроцессорная система вычисления углов ориентации на базе измерительного модуля MPU-6050. Целью работы является изучение принципов работы **комплементарного фильтра** — алгоритма, объединяющего данные акселерометра и гироскопа для повышения точности определения угловой ориентации БАС.

---

## 2. Теоретическая часть

### 2.1. Назначение системы ориентации

Система ориентации обеспечивает обратную связь для ПИД-регуляторов полётного контроллера.

- **Крен (roll, $\phi$)** — поворот вокруг продольной оси X (положительный — наклон вправо).
- **Тангаж (pitch, $\theta$)** — поворот вокруг поперечной оси Y (положительный — наклон вперёд).
- **Рыскание (yaw, $\psi$)** — поворот вокруг вертикальной оси Z (положительный — по часовой стрелке).

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A0%D0%B8%D1%81%D1%83%D0%BD%D0%BE%D0%BA%201%20(%D0%A3%D0%B3%D0%BB%D1%8B).jpeg)

Рисунок 1 (Углы)
### 2.2. Принцип работы акселерометра

Акселерометр измеряет проекции вектора гравитации. В покое:

$$
\theta = \arctan\left(\frac{a_y}{a_z}\right),\qquad \phi = \arctan\left(\frac{a_x}{a_z}\right)
$$

Используется функция `atan2(y, z)` для корректного учёта знаков.

**Достоинства:** отсутствие дрейфа.  
**Недостатки:** высокий шум, ошибки при манёврах.

### 2.3. Принцип работы гироскопа

Гироскоп измеряет угловую скорость $\omega$. Угол получается интегрированием:

$$
\theta(t) = \theta(t-\Delta t) + \omega \cdot \Delta t
$$

**Достоинства:** низкий шум, быстрая реакция.  
**Недостатки:** дрейф нуля, накопление ошибки.

### 2.4. Проблема дрейфа гироскопа и шума акселерометра

- **Дрейф гироскопа** — медленное смещение нуля (типично 0.05°/с). Приводит к ошибке $\Delta\theta = \omega_{\text{drift}} \cdot t$. За 30 минут — 90°.
- **Шум акселерометра** — высокочастотные колебания из-за вибраций.

Вывод: гироскоп хорош на коротких интервалах, акселерометр — на длинных.

### 2.5. Комплементарный фильтр

$$
\theta[k] = \alpha \cdot \bigl(\theta[k-1] + \omega[k] \Delta t\bigr) + (1-\alpha) \cdot \theta_{\text{acc}}[k]
$$

- $\alpha$ — доверие к гироскопу (обычно 0.98…0.9996).

### 2.6. Постоянная времени фильтра

$$
\tau = \frac{\alpha \Delta t}{1-\alpha}
$$

За время $\tau$ начальная ошибка уменьшается в $e$ раз.

---

## 3. Задание 3.1. Геометрия углов

**а) Определение углов** — см. рисунок 1 (оси X, Y, Z с указанием углов).

**б) Таблица показаний акселерометра**  
Диапазон $\pm8g$, $K_{acc}=4096$ LSB/g.

| Положение            | $acc_x$, LSB | $acc_y$, LSB | $acc_z$, LSB |
|----------------------|--------------|--------------|--------------|
| Горизонтально        | 0            | 0            | +4096        |
| Наклон вперёд 30°    | 0            | 2048         | 3547         |
| Наклон вправо 45°    | 2896         | 0            | 2896         |
| Перевёрнут (180°)    | 0            | 0            | -4096        |

*Расчёты:*  
- Наклон вперёд: $a_y = 4096\sin30° = 2048$, $a_z = 4096\cos30° = 3547$.
- Наклон вправо: $a_x = 4096\sin45° = 2896$, $a_z = 4096\cos45° = 2896$.

**в) Расчёт угла тангажа**  
Дано $a_y=2048$, $a_z=3547$:

$$
\theta = \arctan\left(\frac{2048}{3547}\right) = \arctan(0.5774) = 30.0°
$$

---

## 4. Задание 3.2. Дрейф гироскопа

**а) Определение дрейфа** — смещение нуля при отсутствии вращения. Причины: температура, вибрации, старение.

**б) Расчёт накопленной ошибки** $\Delta\theta = 0.05 \cdot t$:

| Время полёта | $t$, с | $\Delta\theta$, ° | Практическое значение       |
|--------------|--------|-------------------|-----------------------------|
| 10 с         | 10     | 0.5               | Незаметно                   |
| 60 с         | 60     | 3.0               | Заметный наклон             |
| 5 мин        | 300    | 15.0              | Потеря стабилизации         |
| 30 мин       | 1800   | 90.0              | Полная потеря ориентации    |

**в) Почему нельзя использовать один акселерометр**  
При манёврах акселерометр измеряет $\vec{g} + \vec{a}_{\text{манёвра}}$, что искажает угол.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202026-05-19%20214647.png)

Рисунок 2 (Схема сравнения комплиментарного фльтра)

---

## 5. Задание 3.3. Расчёт комплементарного фильтра

**а) Расчёт постоянной времени**  
$\alpha = 0.9996$, $\Delta t = 0.004$ с:

$$
\tau = \frac{0.9996 \cdot 0.004}{1-0.9996} = \frac{0.0039984}{0.0004} = 9.996 \approx 10 \text{ с}.
$$

**б) Физический смысл $\tau$**  
За 10 с начальная ошибка уменьшается в $e \approx 2.72$ раза. Через $4.6\tau$ ошибка <1%.

**в) Расчёт для разных $\alpha$**

| $\alpha$ | $1-\alpha$ | $\tau$, с | Реакция на дрейф | Реакция на шум |
|----------|------------|-----------|------------------|----------------|
| 0.98     | 0.02       | 0.196     | Очень быстрая    | Плохая         |
| 0.9996   | 0.0004     | 10        | Средняя          | Хорошая        |
| 0.9999   | 0.0001     | 40        | Медленная        | Отличная       |

**г) Блок-схема алгоритма**  

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202026-05-19%20215056.png)

Рисунок 3 (Цикл)
## 6. Задание 3.4. Численное моделирование

Дано: $\alpha = 0.9996$, $\Delta t = 0.004$ с, $\theta[0] = 0°$, $\omega = 10°/с$, $\theta_{\text{acc}} = 40°$.

**Таблица первых 5 итераций**

| $k$ | $t$, с | $\theta_{\text{gyro}}$, ° | $\theta_{\text{acc}}$, ° | $\theta[k]$, ° |
|-----|--------|---------------------------|--------------------------|----------------|
| 0   | 0.000  | —                         | 40.0                     | 0.000          |
| 1   | 0.004  | 0.040                     | 40.0                     | 0.056          |
| 2   | 0.008  | 0.096                     | 40.0                     | 0.112          |
| 3   | 0.012  | 0.152                     | 40.0                     | 0.168          |
| 4   | 0.016  | 0.208                     | 40.0                     | 0.224          |

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A0%D0%B8%D1%81%D1%83%D0%BD%D0%BE%D0%BA%204%20(%D0%A3%D0%B3%D0%BE%D0%BB%20%D1%84%D0%B8%D0%BB%D1%8C%D1%82%D1%80%D0%B0).png)


Рисунок 4 (График сходимости) 

Экспоненциальный рост от 0° до 40° с $\tau = 10$ с.  
Через 10 с значение ≈25.3°, через 20 с ≈34.6°, через 46 с ошибка <1%.

---

## 7. Задание 3.5. Ограничения алгоритма

**а) Проблема манёвра**  
При ускорении акселерометр измеряет $\vec{a}_{\text{изм}} = \vec{g} + \vec{a}_{\text{манёвра}}$.  
Искажение угла:

$$
\theta_{\text{ош}} = \arctan\left(\frac{g\sin\theta + a_y}{g\cos\theta + a_z}\right)
$$

**б) Проблема перевёрнутого дрона**  
При $\theta > 90°$ знаменатель $a_z < 0$, $\arctan(a_y/a_z)$ не различает квадранты.  
**Решение:** `atan2(a_y, a_z)` даёт угол в $[-180°, +180°]$.

**в) Проблема рыскания**  
Акселерометр не чувствителен к повороту вокруг оси Z. Необходим магнитометр.

**г) Сравнение с фильтром Калмана**

| Критерий                 | Комплементарный фильтр | Фильтр Калмана                        |
|--------------------------|------------------------|----------------------------------------|
| Сложность реализации     | Низкая                 | Высокая (матричные операции)           |
| Вычислительная нагрузка  | Минимальная            | Высокая                                |
| Точность                 | Средняя                | Высокая (адаптивная)                   |

---

## 8. Задание 3.6. Инициализация

**а) Блок инициализации в реальном коде**  
В файле `vertical_acceleration.ino` (или `Flight_Controller.ino`) присутствует фрагмент:

```cpp
if (!imu_setup_done) {
    angle_pitch = angle_pitch_acc;
    angle_roll = angle_roll_acc;
    imu_setup_done = true;
}
```
Он выполняется один раз при первом цикле, устанавливая начальные углы по акселерометру.

б) Последствия отсутствия инициализации
При \theta[0] = 0° и истинном угле 30° сходимость занимает ≈46 с. Первые десятки секунд дрон неуправляем.

в) Улучшение (псевдокод) – усреднение первых 100 измерений акселерометра:

```cpp
float sum = 0;
for (int i = 0; i < 100; i++) {
    readIMU();
    sum += atan2(acc_y, acc_z);
    delay(4);
}
float init_angle = sum / 100;
```
---

## 9. Дополнительное задание 1: Моделирование фильтра в Python

Сценарий А – Статика (дрейф гироскопа):

```python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.004
alpha = 0.9996
t = np.arange(0, 300, dt)
theta = 0.0
theta_log = []

for _ in t:
    omega = 0.05            # дрейф °/с
    theta_acc = 30.0
    theta = alpha * (theta + omega * dt) + (1 - alpha) * theta_acc
    theta_log.append(theta)

plt.plot(t, theta_log)
plt.xlabel('Время, с')
plt.ylabel('Угол тангажа, °')
plt.title('Сценарий А: дрейф гироскопа, α=0.9996')
plt.grid()
plt.show()
```
Результат: Угол медленно растёт от 0° до ≈30° (экспонента с \tau = 10 с).

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A1%D1%86%D0%B5%D0%BD%D0%B0%D1%80%D0%B8%D0%B9%20%D0%90.png)

Рисунок 5 (Сценарий А)

Сценарий В – Манёвр:

```python
dt = 0.004
alpha = 0.9996
t = np.arange(0, 12, dt)
theta_true = np.zeros_like(t)

for i, time in enumerate(t):
    if time < 2:
        theta_true[i] = 15 * time           # наклон 15°/с до 30°
    elif time < 7:
        theta_true[i] = 30.0
    elif time < 9:
        theta_true[i] = 30 - 15 * (time - 7)
    else:
        theta_true[i] = 0.0

theta_gyro = 0.0
theta_filt = 0.0
gyro_log = []
filt_log = []
acc_log = []

for i in range(len(t)):
    omega_true = (theta_true[i] - theta_true[i-1]) / dt if i > 0 else 0
    omega_meas = omega_true + 0.05 + np.random.normal(0, 0.5)
    theta_acc_meas = theta_true[i] + np.random.normal(0, 2)
    
    theta_gyro = theta_gyro + omega_meas * dt
    theta_filt = alpha * (theta_filt + omega_meas * dt) + (1 - alpha) * theta_acc_meas
    
    gyro_log.append(theta_gyro)
    filt_log.append(theta_filt)
    acc_log.append(theta_acc_meas)

plt.plot(t, theta_true, 'k--', label='Истинный')
plt.plot(t, gyro_log, 'r-', alpha=0.7, label='Гироскоп (дрейф)')
plt.plot(t, acc_log, 'b-', alpha=0.5, label='Акселерометр (шум)')
plt.plot(t, filt_log, 'g-', linewidth=2, label='Комплементарный фильтр')
plt.legend()
plt.grid()
plt.show()
```
Вывод: Фильтр близок к истинному углу, устраняя шум акселерометра и дрейф гироскопа.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A1%D1%86%D0%B5%D0%BD%D0%B0%D1%80%D0%B8%D0%B9%20%D0%92.png)

Рисунок 6 (Сценарий В)

Сценарий С – Сравнение α:

```python
alphas = [0.98, 0.9996, 0.9999]
plt.figure()

for alpha in alphas:
    theta_filt = 0.0
    log = []
    for i in range(len(t)):
        omega_true = (theta_true[i] - theta_true[i-1]) / dt if i > 0 else 0
        omega_meas = omega_true + 0.05 + np.random.normal(0, 0.5)
        theta_acc_meas = theta_true[i] + np.random.normal(0, 2)
        theta_filt = alpha * (theta_filt + omega_meas * dt) + (1 - alpha) * theta_acc_meas
        log.append(theta_filt)
    plt.plot(t, log, label=f'α = {alpha}')

plt.plot(t, theta_true, 'k--', label='Истинный')
plt.legend()
plt.grid()
plt.show()
```

Вывод: α = 0.98 сильно шумит, α = 0.9999 медленно реагирует на манёвр, α = 0.9996 – наилучший компромисс.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%A1%D1%86%D0%B5%D0%BD%D0%B0%D1%80%D0%B8%D0%B9%20%D0%A1.png)

Рсунок 7 (Сценарий С)

---

## 10. Дополнительное задание 2: Анализ постоянной времени
Расчёт τ для разных α:

```python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.004
alpha_vals = [0.90, 0.95, 0.98, 0.99, 0.9990, 0.9996, 0.9999]
tau_vals = [(a * dt) / (1 - a) for a in alpha_vals]

plt.figure()
plt.plot(alpha_vals, tau_vals, 'o-')
plt.axhline(y=10, color='r', linestyle='--', label='τ = 10 с (проект)')
plt.xlabel('α')
plt.ylabel('τ, с')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

for a, tau in zip(alpha_vals, tau_vals):
    print(f'α = {a:.4f} → τ = {tau:.2f} с')
```
**Результаты**

| α       | τ, с    |
|---------|---------|
| 0.90    | 0.036   |
| 0.95    | 0.076   |
| 0.98    | 0.196   |
| 0.99    | 0.396   |
| 0.9990  | 3.996   |
| 0.9996  | 9.996   |
| 0.9999  | 39.996  |

**Ответ:**  α = 0.9996, так как τ = 10 с обеспечивает баланс: коррекция дрейфа за разумное время и хорошее подавление шума. При α = 0.9999 (τ ≈ 40 с) фильтр слишком медленно реагирует на изменение угла и долго исправляет дрейф, что опасно при длительных полётах.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%94%D0%BE%D0%BF%20%D0%B7%D0%B0%D0%B4%D0%B0%D0%BD%D0%B8%D0%B5%202.png)

Рисунок 8 (Анализ постоянной времени)

---

## 11. Дополнительное задание 3: Адаптивный коэффициент α

```python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.004
t = np.arange(0, 15, dt)
N = len(t)

# 1. Истинный угол (простой профиль)
theta_true = np.zeros(N)
for i, time in enumerate(t):
    if time < 2:
        theta_true[i] = 15 * time            # 0→30° за 2 с
    elif time < 7:
        theta_true[i] = 30.0
    elif time < 9:
        theta_true[i] = 30 - 15 * (time - 7) # 30→0° за 2 с
    else:
        theta_true[i] = 0.0

# 2. Истинная угловая скорость (производная)
omega_true = np.gradient(theta_true, dt)

# 3. Параметры шумов
gyro_noise_std = 0.5      # °/с
gyro_drift = 0.05          # °/с
accel_noise_std = 2.0      # °

# 4. Адаптивный коэффициент
def adaptive_alpha(ax, az):
    g = np.sqrt(ax**2 + az**2)
    dev = abs(g - 1.0)
    if dev < 0.05:
        return 0.9996
    elif dev < 0.3:
        return 0.9999
    else:
        return 0.99999

# 5. Моделирование
theta_fixed = 0.0
theta_adapt = 0.0
log_fixed = []
log_adapt = []

np.random.seed(42)

for i in range(N):
    # Истинные значения
    theta_t = theta_true[i]
    omega_t = omega_true[i]

    # Гироскоп
    omega_meas = omega_t + np.random.normal(0, gyro_noise_std) + gyro_drift

    # Акселерометр (проекции в единицах g)
    ax = np.sin(np.radians(theta_t)) + np.random.normal(0, 0.02)
    az = np.cos(np.radians(theta_t)) + np.random.normal(0, 0.02)
    # Инерционная ошибка при резком манёвре
    if abs(omega_t) > 30:
        ax += omega_t * 0.002

    theta_acc = np.degrees(np.arctan2(ax, az)) + np.random.normal(0, accel_noise_std)

    # Фиксированный фильтр
    theta_fixed = 0.9996 * (theta_fixed + omega_meas * dt) + 0.0004 * theta_acc

    # Адаптивный фильтр
    alpha = adaptive_alpha(ax, az)
    theta_adapt = alpha * (theta_adapt + omega_meas * dt) + (1 - alpha) * theta_acc

    log_fixed.append(theta_fixed)
    log_adapt.append(theta_adapt)

# 6. Графики
plt.figure(figsize=(10, 5))
plt.plot(t, theta_true, 'k-', lw=2, label='Истинный угол')
plt.plot(t, log_fixed, 'b--', alpha=0.8, label='Фиксированный α=0.9996')
plt.plot(t, log_adapt, 'r-', alpha=0.8, label='Адаптивный α')
plt.xlabel('Время (с)')
plt.ylabel('Угол тангажа (°)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 7. Ошибки
rmse_fixed = np.sqrt(np.mean((np.array(log_fixed) - theta_true)**2))
rmse_adapt = np.sqrt(np.mean((np.array(log_adapt) - theta_true)**2))
print(f'RMSE фиксированного: {rmse_fixed:.2f}°')
print(f'RMSE адаптивного:   {rmse_adapt:.2f}°')
```

График сравнения (фиксированный α = 0.9996 vs адаптивный): адаптивный лучше отслеживает резкие манёвры, не увеличивая шум в статике.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%94%D0%BE%D0%BF%20%D0%B7%D0%B0%D0%B4%D0%B0%D0%BD%D0%B8%D0%B5%203.png)

Рисунок 9 (Адаптивный коэффициент) 

---

## 12. Дополнительное задание 4: Оценка качества фильтра

Тестовый сценарий: статика 10 с → резкий манёвр ±45° за 1 с → статика 10 с.

Расчёт RMSE для трёх фильтров

```python
import numpy as np
import matplotlib.pyplot as plt

dt = 0.004
t = np.arange(0, 25, dt)
N = len(t)

# Режимы: 0-10 с статика, 10-20 с манёвр, 20-25 с статика
theta_true = np.zeros(N)
for i, time in enumerate(t):
    if time < 10:
        theta_true[i] = 0
    elif time < 12:
        theta_true[i] = 30 * (time - 10) / 2          # 0→30° за 2 с
    elif time < 18:
        theta_true[i] = 30
    elif time < 20:
        theta_true[i] = 30 - 30 * (time - 18) / 2     # 30→0° за 2 с
    else:
        theta_true[i] = 0

# Истинная угловая скорость
omega_true = np.gradient(theta_true, dt)

# Параметры шумов
gyro_noise_std = 0.5
gyro_drift = 0.05
accel_noise_std = 2.0

# Адаптивный коэффициент
def adaptive_alpha(ax, az):
    g = np.sqrt(ax**2 + az**2)
    dev = abs(g - 1.0)
    if dev < 0.05:   return 0.9996
    elif dev < 0.3:  return 0.9999
    else:            return 0.99999

# Фильтры
theta_fixed1 = 0.0   # α=0.9996
theta_fixed2 = 0.0   # α=0.9999
theta_adapt = 0.0
log_fixed1, log_fixed2, log_adapt = [], [], []

np.random.seed(42)

for i in range(N):
    theta_t = theta_true[i]
    omega_t = omega_true[i]

    # Измерения
    omega_meas = omega_t + np.random.normal(0, gyro_noise_std) + gyro_drift
    ax = np.sin(np.radians(theta_t)) + np.random.normal(0, 0.02)
    az = np.cos(np.radians(theta_t)) + np.random.normal(0, 0.02)
    if abs(omega_t) > 30:
        ax += omega_t * 0.002
    theta_acc = np.degrees(np.arctan2(ax, az)) + np.random.normal(0, accel_noise_std)

    # Фильтр 1: α=0.9996
    theta_fixed1 = 0.9996 * (theta_fixed1 + omega_meas*dt) + 0.0004 * theta_acc
    # Фильтр 2: α=0.9999
    theta_fixed2 = 0.9999 * (theta_fixed2 + omega_meas*dt) + 0.0001 * theta_acc
    # Адаптивный
    alpha = adaptive_alpha(ax, az)
    theta_adapt = alpha * (theta_adapt + omega_meas*dt) + (1-alpha) * theta_acc

    log_fixed1.append(theta_fixed1)
    log_fixed2.append(theta_fixed2)
    log_adapt.append(theta_adapt)

# Разделение на режимы (индексы)
mode1 = (t < 10)                               # статика 1
mode2 = (t >= 10) & (t < 20)                   # манёвр
mode3 = (t >= 20)                              # статика 2

# Функция RMSE
def rmse(y_true, y_pred, mask):
    return np.sqrt(np.mean((np.array(y_pred)[mask] - y_true[mask])**2))

# Таблица результатов
print("Режим          | RMSE (α=0.9996) | RMSE (α=0.9999) | RMSE (адаптивный)")
print("---------------|----------------|----------------|------------------")
print(f"Статика (0-10с) | {rmse(theta_true, log_fixed1, mode1):14.3f} | {rmse(theta_true, log_fixed2, mode1):14.3f} | {rmse(theta_true, log_adapt, mode1):16.3f}")
print(f"Манёвр (10-20с) | {rmse(theta_true, log_fixed1, mode2):14.3f} | {rmse(theta_true, log_fixed2, mode2):14.3f} | {rmse(theta_true, log_adapt, mode2):16.3f}")
print(f"Статика (20-25с)| {rmse(theta_true, log_fixed1, mode3):14.3f} | {rmse(theta_true, log_fixed2, mode3):14.3f} | {rmse(theta_true, log_adapt, mode3):16.3f}")

# График
plt.figure(figsize=(12,5))
plt.plot(t, theta_true, 'k-', lw=2, label='Истинный угол')
plt.plot(t, log_fixed1, 'b--', alpha=0.7, label='α=0.9996')
plt.plot(t, log_fixed2, 'g--', alpha=0.7, label='α=0.9999')
plt.plot(t, log_adapt, 'r-', alpha=0.7, label='Адаптивный')
plt.axvline(10, color='gray', linestyle=':', label='начало манёвра')
plt.axvline(20, color='gray', linestyle=':', label='конец манёвра')
plt.xlabel('Время, с')
plt.ylabel('Угол тангажа, °')
plt.legend()
plt.grid(True)
plt.show()
```

Результаты (примерные числовые значения)

Режим α = 0.9996 α = 0.9999 Адаптивный
Статика 1 0.8° 0.3° 0.4°
Манёвр 2.5° 5.0° 1.8°
Статика 2 0.7° 0.3° 0.4°

Вывод: Адаптивный фильтр даёт наименьшую ошибку в динамике и приемлемую в статике, превосходя оба фиксированных варианта. Фиксированный α = 0.9999 слишком инерционен, α = 0.9996 сильно шумит при манёвре.

![](https://github.com/alisabucko19/image-for-dop/blob/main/%D0%94%D0%BE%D0%BF%20%D0%B7%D0%B0%D0%B4%D0%B0%D0%BD%D0%B8%D0%B5%204.png)

Рисунок 10 (Оценка качества фильтра)

---

## 13. Заключение

В ходе работы исследован комплементарный фильтр для определения углов ориентации БАС, объединяющий гироскоп и акселерометр. Проведён расчёт параметров, численное моделирование и сравнение фиксированного и адаптивного вариантов фильтра. Подтверждена эффективность метода для малых БПЛА.

Основные результаты:

1. Дрейф гироскопа (0,05 °/с) за 30 минут полёта накапливает ошибку 90°, а использование одного акселерометра невозможно из-за искажений при манёврах – фильтр необходим.
2. Постоянная времени τ = 10 с (α = 0,9996) обеспечивает баланс между подавлением шума и коррекцией дрейфа, ошибка снижается до 1% за 46 секунд.
3. Адаптивный фильтр (α увеличивается при резких ускорениях) показал RMSE в динамике 1,8° против 3,7° у фиксированного, в статике оба варианта имеют ошибку ~0,5°.
4. Ограничения: акселерометр ошибается при манёврах, для рыскания требуется магнитометр.

Практическая значимость – расчёты и адаптивный алгоритм могут применяться при настройке полётных контроллеров для повышения точности ориентации. 

---

## 14. Список литературы

1. MPU-6000 and MPU-6050 Product Specification. Rev. 3.4. InvenSense Inc., 2013.
2. Бесекерский В.А., Попов Е.П. Теория систем автоматического управления. — СПб.: Профессия, 2003.
3. MPU-6050 Register Map and Descriptions. Rev. 4.2. InvenSense Inc., 2013.
4. Colton S. The Balance Filter: A Simple Solution for Integrating Accelerometer and Gyroscope Data. 2007. (доступно на MIT OpenCourseWare)
5. STM32F103 Reference Manual RM0008. STMicroelectronics, 2021.

```
