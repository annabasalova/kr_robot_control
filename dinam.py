import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ahrs.filters import Madgwick
from ahrs.common.orientation import q2euler
from scipy.interpolate import interp1d
import os
import glob

# ==============================================================================
# НАСТРОЙКИ
# ==============================================================================
DATA_FOLDER = r"C:\Users\anna_\PycharmProjects\pythonProject\curs"   # папка с файлами
IMU_FREQ = 400.0
DT_IMU = 1.0 / IMU_FREQ
USE_MAGNETOMETER = True          # рекомендуется True для коррекции рыскания
CALIB_MODE = 'auto'              # 'auto', 'force', 'none'

GPS_ERROR_M = 5.0                # СКО GPS (м)
PROCESS_NOISE = 0.01

SAVE_PLOTS = True
PLOTS_FOLDER = os.path.join(DATA_FOLDER, "plots")
os.makedirs(PLOTS_FOLDER, exist_ok=True)

# ==============================================================================
# ФУНКЦИИ
# ==============================================================================
def need_calibration(acc_vals, gyro_vals):
    """Определяет, нужна ли калибровка (если значения >1000 – значит LSB)"""
    acc_max = np.max(np.abs(acc_vals[:1000]))
    gyro_max = np.max(np.abs(gyro_vals[:1000]))
    need_acc = acc_max > 1000
    need_gyro = gyro_max > 10
    return need_acc, need_gyro

def rotate_vector(v, q):
    """Поворачивает вектор v кватернионом q (из связанной в глобальную)"""
    v_q = np.array([0, v[0], v[1], v[2]])
    qc = q * [1, -1, -1, -1]
    def qmul(a, b):
        w1,x1,y1,z1 = a
        w2,x2,y2,z2 = b
        return np.array([
            w1*w2 - x1*x2 - y1*y2 - z1*z2,
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2
        ])
    return qmul(qmul(q, v_q), qc)[1:]

def plot_results(time_imu, pos_imu, pos_corr, gps_t, gps_x, gps_y, gps_z, quat, filename):
    """Строит и сохраняет графики"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # XY
    axes[0,0].plot(pos_imu[:,0], pos_imu[:,1], 'b-', alpha=0.7, label='IMU only')
    axes[0,0].plot(pos_corr[:,0], pos_corr[:,1], 'r-', label='IMU+GPS Kalman')
    axes[0,0].plot(gps_x, gps_y, 'go', markersize=4, label='GPS')
    axes[0,0].set_xlabel('X (м)'); axes[0,0].set_ylabel('Y (м)')
    axes[0,0].set_title('Горизонтальная траектория')
    axes[0,0].legend(); axes[0,0].grid(True); axes[0,0].axis('equal')
    # XZ
    axes[0,1].plot(pos_imu[:,0], pos_imu[:,2], 'b-', alpha=0.7, label='IMU')
    axes[0,1].plot(pos_corr[:,0], pos_corr[:,2], 'r-', label='IMU+GPS')
    axes[0,1].plot(gps_x, gps_z, 'go', markersize=4, label='GPS')
    axes[0,1].set_xlabel('X (м)'); axes[0,1].set_ylabel('Z (м)')
    axes[0,1].set_title('Вертикальное сечение'); axes[0,1].legend(); axes[0,1].grid(True)
    # X(t)
    axes[1,0].plot(time_imu, pos_imu[:,0], 'b-', alpha=0.7, label='IMU')
    axes[1,0].plot(time_imu, pos_corr[:,0], 'r-', label='IMU+GPS')
    axes[1,0].plot(gps_t, gps_x, 'go', markersize=4, label='GPS')
    axes[1,0].set_xlabel('Время (с)'); axes[1,0].set_ylabel('X (м)')
    axes[1,0].set_title('X(t)'); axes[1,0].legend(); axes[1,0].grid(True)
    # Y(t)
    axes[1,1].plot(time_imu, pos_imu[:,1], 'b-', alpha=0.7, label='IMU')
    axes[1,1].plot(time_imu, pos_corr[:,1], 'r-', label='IMU+GPS')
    axes[1,1].plot(gps_t, gps_y, 'go', markersize=4, label='GPS')
    axes[1,1].set_xlabel('Время (с)'); axes[1,1].set_ylabel('Y (м)')
    axes[1,1].set_title('Y(t)'); axes[1,1].legend(); axes[1,1].grid(True)

    plt.suptitle(f'Результаты фильтрации – {os.path.basename(filename)}')
    plt.tight_layout()
    out_path = os.path.join(PLOTS_FOLDER, f"{os.path.basename(filename)}.png")
    plt.savefig(out_path, dpi=150)
    print(f"График сохранён: {out_path}")
    plt.close()

def safe_float_convert(df, columns):
    """Приводит колонки к float, заменяя нечисловые значения на NaN"""
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def process_file(filepath, calib_mode='auto', use_mag=True, save_plot=False):
    print(f"\n{'='*70}\nОбработка: {os.path.basename(filepath)}\n{'='*70}")
    df = pd.read_excel(filepath, engine='openpyxl')
    df.columns = df.columns.str.strip()

    # Определение временной колонки
    if 'TimeUS' in df.columns:
        df['Time_sec'] = df['TimeUS'] / 1_000_000.0
        time_col = 'TimeUS'
    elif 'TimeS' in df.columns:
        df['Time_sec'] = df['TimeS']
        time_col = 'TimeS'
    else:
        raise ValueError(f"Нет TimeUS или TimeS в {filepath}")
    print(f"Используется временная колонка: {time_col}")
    df = df.sort_values('Time_sec').reset_index(drop=True)

    # Преобразование всех нужных колонок в числа
    numeric_cols = ['ACC.AccX','ACC.AccY','ACC.AccZ',
                    'GYR.GyrX','GYR.GyrY','GYR.GyrZ',
                    'MAG.MagX','MAG.MagY','MAG.MagZ',
                    'GPS.Lat','GPS.Lng','GPS.Alt']
    df = safe_float_convert(df, numeric_cols)

    # Выделение датчиков
    acc_df = df[df['ACC.AccX'].notna()][['Time_sec','ACC.AccX','ACC.AccY','ACC.AccZ']].copy()
    acc_df.rename(columns={'ACC.AccX':'ax','ACC.AccY':'ay','ACC.AccZ':'az'}, inplace=True)
    gyro_df = df[df['GYR.GyrX'].notna()][['Time_sec','GYR.GyrX','GYR.GyrY','GYR.GyrZ']].copy()
    gyro_df.rename(columns={'GYR.GyrX':'gx','GYR.GyrY':'gy','GYR.GyrZ':'gz'}, inplace=True)
    mag_df = df[df['MAG.MagX'].notna()][['Time_sec','MAG.MagX','MAG.MagY','MAG.MagZ']].copy()
    mag_df.rename(columns={'MAG.MagX':'mx','MAG.MagY':'my','MAG.MagZ':'mz'}, inplace=True)
    gps_df = df[df['GPS.Lat'].notna()][['Time_sec','GPS.Lat','GPS.Lng','GPS.Alt']].copy()
    gps_df.rename(columns={'GPS.Lat':'lat','GPS.Lng':'lon','GPS.Alt':'alt'}, inplace=True)

    # Слияние датчиков по времени (допуск 1 мс)
    merged = pd.merge_asof(acc_df, gyro_df, on='Time_sec', direction='nearest', tolerance=0.001)
    if len(mag_df) > 0:
        merged = pd.merge_asof(merged, mag_df, on='Time_sec', direction='nearest', tolerance=0.001)
    else:
        merged['mx'] = merged['my'] = merged['mz'] = 0.0
    merged = merged.dropna(subset=['ax','ay','az','gx','gy','gz'])
    merged = merged.sort_values('Time_sec').reset_index(drop=True)

    # Калибровка (если нужно)
    calib_len = min(1000, len(merged))
    acc_raw = merged[['ax','ay','az']].values.astype(float)
    gyro_raw = merged[['gx','gy','gz']].values.astype(float)
    mag_raw = merged[['mx','my','mz']].values.astype(float)

    if calib_mode == 'auto':
        need_acc, need_gyro = need_calibration(acc_raw, gyro_raw)
        print(f"Авто: акселерометр калибровать? {need_acc}, гироскоп? {need_gyro}")
    elif calib_mode == 'force':
        need_acc, need_gyro = True, True
    else:
        need_acc, need_gyro = False, False

    if calib_len > 0:
        if need_acc:
            acc_bias = np.nanmean(acc_raw[:calib_len], axis=0)
            # Если значения по Z около 9.8, то уже м/с², просто вычитаем bias
            if abs(acc_bias[2]) > 5:
                acc_bias[2] = 0.0
                acc = acc_raw - acc_bias
            else:
                gravity = acc_bias[2]
                acc_bias[2] = 0.0
                acc = (acc_raw - acc_bias) / gravity if gravity != 0 else acc_raw
        else:
            acc = acc_raw

        if need_gyro:
            gyro_bias = np.nanmean(gyro_raw[:calib_len], axis=0)
            gyro = gyro_raw - gyro_bias
            # Если значения большие (>1) – возможно, градусы/с -> переводим в рад/с
            if np.max(np.abs(gyro[:100])) > 1.0:
                gyro = np.radians(gyro)
        else:
            gyro = gyro_raw

        # Магнитометр: всегда вычитаем среднее за период покоя
        mag_bias = np.nanmean(mag_raw[:calib_len], axis=0)
        mag = mag_raw - mag_bias
        mag = np.nan_to_num(mag)
    else:
        acc, gyro, mag = acc_raw, gyro_raw, np.nan_to_num(mag_raw)

    # Интерполяция на равномерную сетку (400 Гц)
    t_min, t_max = merged['Time_sec'].min(), merged['Time_sec'].max()
    time_grid = np.arange(t_min, t_max, DT_IMU)
    def interp3(data, times, tgrid):
        res = np.zeros((len(tgrid),3))
        for i in range(3):
            f = interp1d(times, data[:,i], kind='linear', fill_value='extrapolate')
            res[:,i] = f(tgrid)
        return res
    orig_t = merged['Time_sec'].values
    acc_i = interp3(acc, orig_t, time_grid)
    gyro_i = interp3(gyro, orig_t, time_grid)
    mag_i = interp3(mag, orig_t, time_grid)

    # Фильтр Маджвика
    mad = Madgwick(frequency=IMU_FREQ, beta=0.1)
    n = len(time_grid)
    quat = np.zeros((n,4))
    quat[0] = [1,0,0,0]
    for i in range(1,n):
        if use_mag:
            quat[i] = mad.updateMARG(quat[i-1], gyr=gyro_i[i], acc=acc_i[i], mag=mag_i[i])
        else:
            quat[i] = mad.updateIMU(quat[i-1], gyr=gyro_i[i], acc=acc_i[i])

    # Интегрирование ускорения
    pos = np.zeros((n,3))
    vel = np.zeros((n,3))
    g = np.array([0,0,9.80665])
    for i in range(1,n):
        a_glob = rotate_vector(acc_i[i], quat[i]) - g
        vel[i] = vel[i-1] + a_glob * DT_IMU
        pos[i] = pos[i-1] + vel[i-1] * DT_IMU + 0.5 * a_glob * DT_IMU**2

    # GPS в метры
    if len(gps_df) < 2:
        print("Недостаточно GPS-точек")
        return None
    gps_df = gps_df.dropna(subset=['lat','lon','alt']).sort_values('Time_sec')
    if len(gps_df) == 0:
        return None
    lat0, lon0, alt0 = gps_df.iloc[0][['lat','lon','alt']]
    def gps2m(lat, lon, alt):
        dy = (lat - lat0) * 111320.0
        dx = (lon - lon0) * 111320.0 * np.cos(np.radians((lat+lat0)/2))
        dz = alt - alt0
        return dx, dy, dz
    gps_xyz = np.array([gps2m(row.lat, row.lon, row.alt) for _, row in gps_df.iterrows()])
    gps_t = gps_df['Time_sec'].values
    gps_x, gps_y, gps_z = gps_xyz[:,0], gps_xyz[:,1], gps_xyz[:,2]

    # Фильтр Калмана (прогноз с ускорением)
    def kf_predict(x, P, dt, acc):
        F = np.eye(6)
        F[0,3]=dt; F[1,4]=dt; F[2,5]=dt
        B = np.zeros((6,3))
        B[0,0]=0.5*dt**2; B[1,1]=0.5*dt**2; B[2,2]=0.5*dt**2
        B[3,0]=dt; B[4,1]=dt; B[5,2]=dt
        x = F @ x + B @ acc
        Q = np.eye(6) * PROCESS_NOISE
        P = F @ P @ F.T + Q
        return x, P
    def kf_update(x, P, z, R):
        H = np.zeros((3,6)); H[0,0]=1; H[1,1]=1; H[2,2]=1
        y = z - H @ x
        S = H @ P @ H.T + R
        K = P @ H.T @ np.linalg.inv(S)
        x = x + K @ y
        P = (np.eye(6) - K @ H) @ P
        return x, P

    init_pos = np.mean(pos[:100], axis=0)
    xk = np.array([init_pos[0], init_pos[1], init_pos[2], 0,0,0])
    Pk = np.eye(6) * 10.0
    Rk = np.eye(3) * (GPS_ERROR_M**2)
    pos_corr = np.zeros((n,3))
    pos_corr[0] = init_pos
    gps_idx = 0
    for i in range(1,n):
        a_glob = rotate_vector(acc_i[i], quat[i]) - g
        xk, Pk = kf_predict(xk, Pk, DT_IMU, a_glob)
        if gps_idx < len(gps_t) and abs(time_grid[i] - gps_t[gps_idx]) < 0.5*DT_IMU:
            z_gps = np.array([gps_x[gps_idx], gps_y[gps_idx], gps_z[gps_idx]])
            xk, Pk = kf_update(xk, Pk, z_gps, Rk)
            gps_idx += 1
        pos_corr[i] = xk[:3]

    # Ошибки в точках GPS
    errors = []
    for j, t in enumerate(gps_t):
        idx = np.argmin(np.abs(time_grid - t))
        if abs(time_grid[idx] - t) < 0.05:
            e_imu = np.linalg.norm(pos[idx] - [gps_x[j], gps_y[j], gps_z[j]])
            e_corr = np.linalg.norm(pos_corr[idx] - [gps_x[j], gps_y[j], gps_z[j]])
            errors.append((e_imu, e_corr))
    if not errors:
        return None
    e_imu_arr = [e[0] for e in errors]
    e_corr_arr = [e[1] for e in errors]

    # Дрейф рыскания
    yaw = np.degrees([q2euler(q)[2] for q in quat])
    yaw_start = np.mean(yaw[:int(IMU_FREQ)])
    yaw_end = np.mean(yaw[-int(IMU_FREQ):])
    gps_displacement = np.linalg.norm([gps_x[-1]-gps_x[0], gps_y[-1]-gps_y[0], gps_z[-1]-gps_z[0]])

    result = {
        'file': os.path.basename(filepath),
        'duration': time_grid[-1] - time_grid[0],
        'gps_points': len(errors),
        'mean_err_imu': np.mean(e_imu_arr),
        'max_err_imu': np.max(e_imu_arr),
        'mean_err_corr': np.mean(e_corr_arr),
        'max_err_corr': np.max(e_corr_arr),
        'yaw_drift': yaw_end - yaw_start,
        'gps_displacement': gps_displacement
    }

    print(f"Длит: {result['duration']:.1f} с, GPS-точек: {result['gps_points']}")
    print(f"Ошибка IMU: средняя {result['mean_err_imu']:.1f} м, макс {result['max_err_imu']:.0f} м")
    print(f"Ошибка коррекции: средняя {result['mean_err_corr']:.2f} м, макс {result['max_err_corr']:.2f} м")
    print(f"Дрейф рыскания: {result['yaw_drift']:.1f}°, GPS перемещение: {result['gps_displacement']:.1f} м")

    if save_plot:
        plot_results(time_grid, pos, pos_corr, gps_t, gps_x, gps_y, gps_z, quat, filepath)

    return result

# ==============================================================================
# ОСНОВНОЙ ЗАПУСК
# ==============================================================================
if __name__ == "__main__":
    if not os.path.exists(DATA_FOLDER):
        print(f"Папка не найдена: {DATA_FOLDER}")
        exit(1)

    # Найти все файлы .xlsx, исключая raw_data
    all_files = glob.glob(os.path.join(DATA_FOLDER, "*.xlsx"))
    file_list = [f for f in all_files if "raw_data" not in os.path.basename(f).lower()]
    print(f"Найдено файлов: {len(file_list)}")
    if not file_list:
        print("Нет файлов для обработки.")
        exit(0)

    all_res = []
    for f in file_list:
        # Быстрая проверка наличия временной колонки
        try:
            temp = pd.read_excel(f, nrows=0, engine='openpyxl')
            temp.columns = temp.columns.str.strip()
            if not ('TimeUS' in temp.columns or 'TimeS' in temp.columns):
                print(f"Файл {f} не содержит TimeUS или TimeS, пропускаем")
                continue
        except Exception as e:
            print(f"Ошибка при проверке {f}: {e}")
            continue

        res = process_file(f, calib_mode=CALIB_MODE, use_mag=USE_MAGNETOMETER, save_plot=SAVE_PLOTS)
        if res:
            all_res.append(res)

    if all_res:
        print("\n" + "="*70)
        print("СРЕДНИЕ РЕЗУЛЬТАТЫ ПО ВСЕМ ФАЙЛАМ")
        print("="*70)
        metrics = ['mean_err_imu', 'max_err_imu', 'mean_err_corr', 'max_err_corr', 'yaw_drift', 'gps_displacement']
        for m in metrics:
            vals = [r[m] for r in all_res]
            print(f"{m}: среднее = {np.mean(vals):.2f} ± {np.std(vals):.2f}")
    else:
        print("Нет обработанных файлов")
