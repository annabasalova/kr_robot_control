#include <bits/stdc++.h>
#include "quaternion.h"
#include "matrix.h"

using namespace std;

const double g = 9.80666;
const double PI = 3.14159265358979323846;

// Коэффициенты преобразования из LSB в физические единицы
const double ACCEL_LSB_TO_MS2 = 9.80665 / 16384.0; // для акселерометра ±2g
const double GYRO_LSB_TO_RADS = (250.0 * PI / 180.0) / 32768.0; // для гироскопа ±250°/с

// Функции преобразования из LSB в физические единицы
pure convertAccelerometerLSBToMS2(pure acc_lsb) {
    return pure(acc_lsb.x * ACCEL_LSB_TO_MS2,
                acc_lsb.y * ACCEL_LSB_TO_MS2,
                acc_lsb.z * ACCEL_LSB_TO_MS2);
}

pure convertGyroscopeLSBToRadS(pure gyro_lsb) {
    return pure(gyro_lsb.x * GYRO_LSB_TO_RADS,
                gyro_lsb.y * GYRO_LSB_TO_RADS,
                gyro_lsb.z * GYRO_LSB_TO_RADS);
}

class ExtendedKalmanFilter {
private:
    Matrix P; // Ковариационная матрица ошибки
    Matrix Q; // Ковариация шума процесса
    Matrix R; // Ковариация шума измерений
    Matrix x; // Вектор состояния [px, py, pz, vx, vy, vz, q0, q1, q2, q3]

    double dt;

public:
    ExtendedKalmanFilter(double _dt) :
            P(10, 10), Q(10, 10), R(9, 9), x(10, 1), dt(_dt) {

        // Инициализация ковариационных матриц
        for(int i = 0; i < 10; i++) P.data[i][i] = 10.0;
        for(int i = 0; i < 10; i++) Q.data[i][i] = 0.1;
        for(int i = 0; i < 9; i++) R.data[i][i] = 1.0;

        // Начальное состояние: позиция (0,0,0), скорость (0,0,0), ориентация (1,0,0,0)
        for(int i = 0; i < 10; i++) x.data[i][0] = 0.0;
        x.data[6][0] = 1.0; // q0 = 1
    }

    // Функция предсказания (модель движения)
    void predict(pure acceleration, pure gyro) {
        // Извлекаем текущее состояние
        double px = x.data[0][0], py = x.data[1][0], pz = x.data[2][0];
        double vx = x.data[3][0], vy = x.data[4][0], vz = x.data[5][0];
        double q0 = x.data[6][0], q1 = x.data[7][0], q2 = x.data[8][0], q3 = x.data[9][0];

        // Преобразуем ускорение из локальной системы в глобальную с учетом ориентации
        pure acc_global = rotateVectorByQuaternion(acceleration,
                                                   quat(q0, q1, q2, q3));

        // Вычитаем гравитацию
        acc_global.z -= g;

        // Обновляем позицию и скорость (линейная часть)
        x.data[0][0] = px + vx * dt + 0.5 * acc_global.x * dt * dt;
        x.data[1][0] = py + vy * dt + 0.5 * acc_global.y * dt * dt;
        x.data[2][0] = pz + vz * dt + 0.5 * acc_global.z * dt * dt;

        x.data[3][0] = vx + acc_global.x * dt;
        x.data[4][0] = vy + acc_global.y * dt;
        x.data[5][0] = vz + acc_global.z * dt;

        // Обновляем ориентацию (кватернион) с учетом гироскопа
        updateOrientation(gyro);

        // Обновляем ковариационную матрицу (упрощенно)
        P = P + Q;
    }

    // Функция обновления (коррекция по измерениям)
    void update(pure accelerometer, pure gyroscope, pure magnetometer) {
        // Вектор измерений [ax, ay, az, gx, gy, gz, mx, my, mz]
        Matrix z(9, 1);
        z.data[0][0] = accelerometer.x;
        z.data[1][0] = accelerometer.y;
        z.data[2][0] = accelerometer.z;
        z.data[3][0] = gyroscope.x;
        z.data[4][0] = gyroscope.y;
        z.data[5][0] = gyroscope.z;
        z.data[6][0] = magnetometer.x;
        z.data[7][0] = magnetometer.y;
        z.data[8][0] = magnetometer.z;

        // Предсказанные измерения (упрощенно)
        Matrix h(9, 1);
        // Для акселерометра - ускорение в локальной системе с учетом гравитации
        quat current_orientation = getOrientation();
        pure gravity_local = rotateVectorByQuaternion(pure(0, 0, g), current_orientation.conjugate());
        h.data[0][0] = gravity_local.x;
        h.data[1][0] = gravity_local.y;
        h.data[2][0] = gravity_local.z;

        // Для гироскопа - угловая скорость (предсказание равно нулю в стационарном состоянии)
        h.data[3][0] = 0;
        h.data[4][0] = 0;
        h.data[5][0] = 0;

        // Для магнитометра - магнитное поле (упрощенно)
        pure north(1, 0, 0); // предполагаем, что север по оси X
        pure mag_local = rotateVectorByQuaternion(north, current_orientation.conjugate());
        h.data[6][0] = mag_local.x;
        h.data[7][0] = mag_local.y;
        h.data[8][0] = mag_local.z;

        // Матрица Якоби H (упрощенно)
        Matrix H(9, 10);
        for(int i = 0; i < 9; i++) {
            if(i < 6) H.data[i][i] = 1.0; // Простая диагональная матрица
        }

        // Вычисление коэффициента Калмана
        Matrix Ht = H.transpose();
        Matrix S = (H * P * Ht) + R;
        Matrix K = P * Ht * S.inverse();

        // Обновление состояния
        Matrix innovation = z - h;
        x = x + (K * innovation);

        // Обновление ковариации
        Matrix I = Matrix::identity(10);
        P = (I - (K * H)) * P;

        // Нормализуем кватернион после обновления
        normalizeQuaternion();
    }

    vector<double> getState() const {
        vector<double> result(10);
        for(int i = 0; i < 10; i++) result[i] = x.data[i][0];
        return result;
    }

    pure getPosition() const {
        return pure(x.data[0][0], x.data[1][0], x.data[2][0]);
    }

    pure getVelocity() const {
        return pure(x.data[3][0], x.data[4][0], x.data[5][0]);
    }

    quat getOrientation() const {
        return quat(x.data[6][0], x.data[7][0], x.data[8][0], x.data[9][0]);
    }

private:
    void updateOrientation(pure gyro) {
        double q0 = x.data[6][0], q1 = x.data[7][0], q2 = x.data[8][0], q3 = x.data[9][0];
        double wx = gyro.x, wy = gyro.y, wz = gyro.z;

        // Производная кватерниона от угловой скорости
        double dq0 = 0.5 * (-q1*wx - q2*wy - q3*wz);
        double dq1 = 0.5 * ( q0*wx + q2*wz - q3*wy);
        double dq2 = 0.5 * ( q0*wy - q1*wz + q3*wx);
        double dq3 = 0.5 * ( q0*wz + q1*wy - q2*wx);

        // Интегрирование
        x.data[6][0] = q0 + dq0 * dt;
        x.data[7][0] = q1 + dq1 * dt;
        x.data[8][0] = q2 + dq2 * dt;
        x.data[9][0] = q3 + dq3 * dt;

        // Нормализация кватерниона
        normalizeQuaternion();
    }

    void normalizeQuaternion() {
        double q0 = x.data[6][0], q1 = x.data[7][0], q2 = x.data[8][0], q3 = x.data[9][0];
        double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);

        if(norm > 1e-10) {
            x.data[6][0] /= norm;
            x.data[7][0] /= norm;
            x.data[8][0] /= norm;
            x.data[9][0] /= norm;
        }
    }

    pure rotateVectorByQuaternion(pure v, quat q) {
        // Поворот вектора v кватернионом q
        quat p(0, v.x, v.y, v.z);
        quat q_conj = q.conjugate();
        quat rotated = q * p * q_conj;
        return pure(rotated.q1, rotated.q2, rotated.q3);
    }
};

struct state {
    pure a; // ускорение
    pure v; // скорость
    pure p; // позиция
    pure w; // с магнитометра
    pure g; // с гироскопа

    state() : a(0), v(0), p(0) {}
    state(pure _a, pure _v, pure _p) : a(_a), v(_v), p(_p) {}
};

std::ostream& operator << (std::ostream &os, const state &st) {
    return os << st.a << " " << st.v << " " << st.p;
}

std::istream& operator >> (std::istream& in, state& st) {
    in >> st.a >> st.v >> st.p;
    return in;
}

struct position {
    state cur_state;
    ExtendedKalmanFilter ekf;
    double dt = 10; // время в мс
    bool first;

    position() : ekf(10.0/1000.0) { // преобразуем в секунды
        cur_state.a = 0;
        cur_state.v = 0;
        cur_state.p = 0;
        first = true;
    }

    void update(pure acc_lsb, pure gyro_lsb, pure mag_lsb) {
        if(first) {
            first = false;
            // Инициализация фильтра первыми измерениями
            return;
        }

        // Преобразование из LSB в физические единицы
        pure accelerometer = convertAccelerometerLSBToMS2(acc_lsb);
        pure gyroscope = convertGyroscopeLSBToRadS(gyro_lsb);
        pure magnetometer = mag_lsb; // магнитометр оставляем как есть (предполагаем уже в правильных единицах)

        // Вызов функции предсказания фильтра Калмана
        ekf.predict(accelerometer, gyroscope);

        // Вызов функции обновления фильтра Калмана
        ekf.update(accelerometer, gyroscope, magnetometer);

        // Обновление текущего состояния из фильтра
        cur_state.p = ekf.getPosition();
        cur_state.v = ekf.getVelocity();
        cur_state.a = accelerometer; // Используем преобразованные измерения ускорения
    }
};

class Matrix2D {
public:
    vector<vector<double>> data;
    int rows, cols;

    Matrix2D(int r, int c) : rows(r), cols(c) {
        data.resize(r, vector<double>(c, 0.0));
    }

    Matrix2D operator*(const Matrix2D& other) const {
        Matrix2D result(rows, other.cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < other.cols; j++) {
                for(int k = 0; k < cols; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix2D transpose() const {
        Matrix2D result(cols, rows);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    Matrix2D operator+(const Matrix2D& other) const {
        Matrix2D result(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix2D operator-(const Matrix2D& other) const {
        Matrix2D result(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix2D inverse() const {
        // Упрощенная реализация для диагональных матриц
        Matrix2D result(rows, cols);
        for(int i = 0; i < min(rows, cols); i++) {
            if(fabs(data[i][i]) > 1e-10) {
                result.data[i][i] = 1.0 / data[i][i];
            }
        }
        return result;
    }

    static Matrix2D identity(int n) {
        Matrix2D result(n, n);
        for(int i = 0; i < n; i++) {
            result.data[i][i] = 1.0;
        }
        return result;
    }
};

class KalmanFilter2D {
private:
    Matrix2D P; // Ковариационная матрица ошибки
    Matrix2D Q; // Ковариация шума процесса
    Matrix2D R; // Ковариация шума измерений
    Matrix2D x; // Вектор состояния [px, py, vx, vy]

    double dt;

public:
    KalmanFilter2D(double _dt) :
            P(4, 4), Q(4, 4), R(2, 2), x(4, 1), dt(_dt) {

        // Инициализация ковариационных матриц
        for(int i = 0; i < 4; i++) P.data[i][i] = 10.0;
        for(int i = 0; i < 4; i++) Q.data[i][i] = 0.1;
        for(int i = 0; i < 2; i++) R.data[i][i] = 1.0;

        // Начальное состояние: позиция (0,0), скорость (0,0)
        for(int i = 0; i < 4; i++) x.data[i][0] = 0.0;
    }

    // Функция предсказания (модель движения)
    void predict(pure acceleration) {
        // Матрица состояния A для 2D модели [px, py, vx, vy]
        Matrix2D A(4, 4);
        for(int i = 0; i < 4; i++) A.data[i][i] = 1.0;
        A.data[0][2] = dt;
        A.data[1][3] = dt;

        // Матрица управления B
        Matrix2D B(4, 2);
        B.data[0][0] = 0.5 * dt * dt;
        B.data[1][1] = 0.5 * dt * dt;
        B.data[2][0] = dt;
        B.data[3][1] = dt;

        // Вектор управления u [ax, ay]
        Matrix2D u(2, 1);
        u.data[0][0] = acceleration.x;
        u.data[1][0] = acceleration.y;

        // Предсказание состояния: x = A*x + B*u
        x = (A * x) + (B * u);

        // Обновление ковариационной матрицы: P = A*P*A^T + Q
        P = (A * P * A.transpose()) + Q;
    }

    // Функция обновления (коррекция по измерениям позиции)
    void update(pure position_measurement) {
        // Матрица измерений H (измеряем только позицию x,y)
        Matrix2D H(2, 4);
        H.data[0][0] = 1.0;
        H.data[1][1] = 1.0;

        // Вектор измерений [px, py]
        Matrix2D z(2, 1);
        z.data[0][0] = position_measurement.x;
        z.data[1][0] = position_measurement.y;

        // Вычисление коэффициента Калмана
        Matrix2D Ht = H.transpose();
        Matrix2D S = (H * P * Ht) + R;
        Matrix2D K = P * Ht * S.inverse();

        // Обновление состояния
        Matrix2D innovation = z - (H * x);
        x = x + (K * innovation);

        // Обновление ковариации
        Matrix2D I = Matrix2D::identity(4);
        P = (I - (K * H)) * P;
    }

    vector<double> getState() const {
        vector<double> result(4);
        for(int i = 0; i < 4; i++) result[i] = x.data[i][0];
        return result;
    }

    pure getPosition() const {
        return pure(x.data[0][0], x.data[1][0], 0);
    }

    pure getVelocity() const {
        return pure(x.data[2][0], x.data[3][0], 0);
    }
};

// Простой комплементарный фильтр для ориентации
class OrientationFilter2D {
private:
    double yaw; // угол поворота вокруг оси Z (в радианах)
    double dt;

public:
    OrientationFilter2D(double _dt) : yaw(0), dt(_dt) {}

    void update(pure gyro, pure mag) {
        // Интегрируем гироскоп для получения угла
        yaw += gyro.z * dt;

        // Корректируем угол по магнитометру (простой комплементарный фильтр)
        double mag_yaw = atan2(mag.y, mag.x);
        yaw = 0.98 * yaw + 0.02 * mag_yaw;
    }

    double getYaw() const {
        return yaw;
    }

    // Поворот вектора из локальной системы в глобальную
    pure rotateToGlobal(pure local) const {
        double cos_yaw = cos(yaw);
        double sin_yaw = sin(yaw);

        pure global;
        global.x = local.x * cos_yaw - local.y * sin_yaw;
        global.y = local.x * sin_yaw + local.y * cos_yaw;
        global.z = 0;
        return global;
    }
};

struct state_2d {
    pure a; // ускорение в плоскости
    pure v; // скорость в плоскости
    pure p; // позиция в плоскости

    state_2d() : a(0), v(0), p(0) {}
    state_2d(pure _a, pure _v, pure _p) : a(_a), v(_v), p(_p) {}
};

std::ostream& operator << (std::ostream &os, const state_2d &st) {
    return os << st.p.x << " " << st.p.y;
}


struct position_2d {
    state_2d cur_state;
    KalmanFilter2D kf;
    OrientationFilter2D orientation_filter;
    double dt = 20; // время в мс
    bool first;

    position_2d() : kf(20.0/1000.0), orientation_filter(20.0/1000.0) {
        cur_state.a = 0;
        cur_state.v = 0;
        cur_state.p = 0;
        first = true;
    }

    void update(pure acc_lsb, pure gyro_lsb, pure mag_lsb) {
        if(first) {
            first = false;
            return;
        }

        // Преобразование из LSB в физические единицы
        pure accelerometer = convertAccelerometerLSBToMS2(acc_lsb);
        pure gyroscope = convertGyroscopeLSBToRadS(gyro_lsb);

        // Обновление ориентации
        orientation_filter.update(gyroscope, mag_lsb);

        // Поворот ускорения из локальной системы в глобальную
        pure acc_global = orientation_filter.rotateToGlobal(accelerometer);

        // Компенсация гравитации (предполагаем, что устройство ориентировано горизонтально)
        // В 2D мы не учитываем наклоны, поэтому просто используем акселерометр для движения

        // Предсказание шага фильтра Калмана
        kf.predict(acc_global);

        // Создаем измерение позиции через простое интегрирование (для демонстрации)
        pure measured_position = cur_state.p;
        measured_position.x += cur_state.v.x * dt + 0.5 * acc_global.x * dt * dt;
        measured_position.y += cur_state.v.y * dt + 0.5 * acc_global.y * dt * dt;

        // Коррекция фильтра Калмана
        kf.update(measured_position);

        // Обновление текущего состояния из фильтра
        cur_state.p = kf.getPosition();
        cur_state.v = kf.getVelocity();
        cur_state.a = acc_global;
    }
};


signed main() {
    ifstream in("data.txt");
    ofstream out("out.txt");
    position cur_pos;

    while(!in.eof()) {
        pure accelerometer_lsb, gyroscope_lsb, magnetometer_lsb;

        // Чтение данных акселерометра в LSB
        in >> accelerometer_lsb;

        // Чтение данных гироскопа в LSB
        in >> gyroscope_lsb;

        // Чтение данных магнитометра
        in >> magnetometer_lsb;

        // Вызов функции update с данными в LSB
        cur_pos.update(accelerometer_lsb, gyroscope_lsb, magnetometer_lsb);

        // Вывод отфильтрованной позиции
        cout << fixed << setprecision(10) << cur_pos.cur_state.p << "\n";
    }

    in.close();
    out.close();
    return 0;
}
