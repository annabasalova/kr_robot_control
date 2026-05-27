import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ahrs.filters import Madgwick, Mahony, EKF
from ahrs.common.orientation import q2euler

# ==============================================================================
# 1. Загрузка и подготовка данных
# ==============================================================================
FILE_PATH = "data"
SAMPLE_FREQ = 50.0


def load_imu_data(filepath, freq):
    """Загружает 9 столбцов: магнитометр, акселерометр, гироскоп."""
    data = pd.read_csv(filepath, header=None, sep=None, engine='python')

    if data.shape[1] != 9:
        raise ValueError(f"Файл должен содержать 9 столбцов, найдено {data.shape[1]}")

    mag = data.iloc[:, 0:3].values.astype(np.float64)  # Магнитометр
    acc = data.iloc[:, 3:6].values.astype(np.float64)  # Акселерометр
    gyro = data.iloc[:, 6:9].values.astype(np.float64)  # Гироскоп

    dt = 1.0 / freq
    time = np.arange(len(data)) * dt

    acc_mean = acc.mean(axis=0)
    gravity_lsb = acc_mean[2]
    acc_bias = np.array([acc_mean[0], acc_mean[1], 0.0])
    acc = (acc - acc_bias) / gravity_lsb
    gyro_mean = gyro.mean(axis=0)
    gyro = np.radians(gyro - gyro_mean)
    mag_mean = mag.mean(axis=0)
    mag = mag - mag_mean

    return gyro, acc, mag, time, dt


print("Загрузка данных...")
gyro, acc, mag, time, dt = load_imu_data(FILE_PATH, SAMPLE_FREQ)
n_samples = len(time)
print(f"Загружено {n_samples} отсчётов, длительность {time[-1]:.2f} сек")

# ==============================================================================
# 2. Инициализация фильтров
# ==============================================================================
madgwick_filter = Madgwick(frequency=SAMPLE_FREQ, beta=0.1)
mahony_filter = Mahony(frequency=SAMPLE_FREQ, K_P=0.5, K_I=0.05)
ekf_filter = EKF(frequency=SAMPLE_FREQ)

Q_madgwick = np.zeros((n_samples, 4))
Q_mahony = np.zeros((n_samples, 4))
Q_ekf = np.zeros((n_samples, 4))

Q_madgwick[0] = np.array([1.0, 0.0, 0.0, 0.0])
Q_mahony[0] = np.array([1.0, 0.0, 0.0, 0.0])
Q_ekf[0] = np.array([1.0, 0.0, 0.0, 0.0])

# ==============================================================================
# 3. Основной цикл фильтрации (обработка данных по шагам)
# ==============================================================================
print("Выполняется фильтрация данных...")
for i in range(1, n_samples):
    g = gyro[i]  # Гироскоп
    a = acc[i]  # Акселерометр
    m = mag[i]  # Магнитометр

    #Q_madgwick[i] = madgwick_filter.updateMARG(Q_madgwick[i - 1], gyr=g, acc=a, mag = m)
    Q_madgwick[i] = madgwick_filter.updateIMU(Q_madgwick[i - 1], gyr=g, acc=a)
    #Q_mahony[i] = mahony_filter.updateMARG(Q_mahony[i - 1], gyr=g, acc=a, mag = m)
    Q_mahony[i] = mahony_filter.updateIMU(Q_mahony[i - 1], gyr=g, acc=a)
    Q_ekf[i] = ekf_filter.update(Q_ekf[i - 1], gyr=g, acc=a)

print("Фильтрация завершена.")

# ==============================================================================
# 4. Преобразование кватернионов в углы Эйлера (в градусы)
# ==============================================================================
euler_madgwick = np.degrees(np.array([q2euler(q) for q in Q_madgwick]))
euler_mahony = np.degrees(np.array([q2euler(q) for q in Q_mahony]))
euler_ekf = np.degrees(np.array([q2euler(q) for q in Q_ekf]))

# ==============================================================================
# 5. Визуализация результатов
# ==============================================================================
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
titles = ['Крен (Roll)', 'Тангаж (Pitch)', 'Рыскание (Yaw)']

for i, (title, ax) in enumerate(zip(titles, axes)):
    ax.plot(time, euler_madgwick[:, i], label='Madgwick (MARG)', linewidth=1.5)
    ax.plot(time, euler_mahony[:, i], label='Mahony (MARG)', linewidth=1.5)
    ax.plot(time, euler_ekf[:, i], label='EKF', linewidth=1.5)
    ax.set_ylabel(f'{title} [°]')
    ax.legend(loc='best')
    ax.grid(True)

axes[-1].set_xlabel('Время [с]')
plt.suptitle('Сравнение алгоритмов оценки ориентации (дрон неподвижен)')
plt.tight_layout()
plt.show()

# ==============================================================================
# 6. Оценка стабильности (для неподвижного дрона)
# ==============================================================================
print("\n=== Статистика углов (градусы) ===")
for name, euler in zip(['Madgwick', 'Mahony', 'EKF'],
                       [euler_madgwick, euler_mahony, euler_ekf]):
    print(f"\n{name}:")
    print(f"  Roll:  mean = {np.mean(euler[:, 0]):6.2f}°,  std = {np.std(euler[:, 0]):.3f}°")
    print(f"  Pitch: mean = {np.mean(euler[:, 1]):6.2f}°,  std = {np.std(euler[:, 1]):.3f}°")
    print(f"  Yaw:   mean = {np.mean(euler[:, 2]):6.2f}°,  std = {np.std(euler[:, 2]):.3f}°")

    print("Madgwick Yaw: min =", euler_madgwick[:, 2].min(), "max =", euler_madgwick[:, 2].max())
    print("Mahony Yaw:   min =", euler_mahony[:, 2].min(), "max =", euler_mahony[:, 2].max())
    print("EKF Yaw:      min =", euler_ekf[:, 2].min(), "max =", euler_ekf[:, 2].max())


