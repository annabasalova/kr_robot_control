//
// Created by anna_ on 20.10.2025.
//

#ifndef KR_3_QUATERNION_H
#define KR_3_QUATERNION_H

#endif //KR_3_QUATERNION_H

#define double long double

struct pure{
    double x,y,z;
    pure(){
        x = 0;
        y = 0;
        z = 0;
    }
    pure(int a){
        x= a;
        y = a;
        x = a;
    }
    pure(double _x, double _y, double _z){
        x = _x;
        y = _y;
        z = _z;
    }
    pure operator =(pure pp){
        x = pp.x;
        y = pp.y;
        z = pp.z;
        return *this;
    }
    pure operator =(double a){
        x = a;
        y = a;
        z = a;
        return *this;
    }
};

struct quat{
    double q0,q1,q2,q3;//-,x,y,z
    pure p;//чистый кватернион
    quat(){
        q0=0;q1=0;q2=0;q3=0;
        p = {0,0,0};
    }
    quat(double aa, double bb, double cc, double dd){
        q0=aa;q1=bb;q2=cc;q3=dd;
        p = {bb,cc,dd};
    }
    quat(pure qq){
        q0 = 0;
        q1 = qq.x;
        q2 = qq.y;
        q3 = qq.z;
    }
    double scal(){
        return q0;
    }
    quat pure_quat(){
        return {0,q1,q2,q3};
    }
    double norm(){
        return q0*q0+q1*q1+q2*q2+q3*q3;
    }
    double abs(){
        return sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
    }
    double abs_vect(){
        return sqrt(q1*q1+q2*q2+q3*q3);
    }
    double arg(){
        return atan2(q0,abs_vect());
    }
    quat vers(){
        return{q0/abs(), q1/abs(), q2/abs(), q3/abs()};
    }
    quat conjugate(){
        return {q0,-q1,-q2,-q3};
    }
    quat operator + (quat q){
        q0 += q.q0;
        q1 += q.q1;
        q2 += q.q2;
        q3 += q.q3;
        return *this;
    }
    quat operator - (quat q){
        q0 -= q.q0;
        q1 -= q.q1;
        q2 -= q.q2;
        q3 -= q.q3;
        return *this;
    }
    quat operator +(double d){
        q0 += d;
        return *this;
    }
    quat operator -(double d){
        q0 -= d;
        return *this;
    }
    quat operator /(double d){
        q0/=d;
        q1/=d;
        q2/=d;
        q3/=d;
        return *this;
    }
    quat operator *(double d){
        q0 *=d;
        q1 *=d;
        q2 *=d;
        q3 *=d;
        return *this;
    }
    quat operator *(quat p){
        double a0,a1,a2,a3;
        a0 = q0*p.q0 - q1*p.q1 -  q2*p.q2 - q3*p.q3;
        a1 = q0*p.q1 + q1*p.q0 + q2*p.q3 - q3*p.q2;
        a2 = q0*p.q2 - q1*p.q3 + q2*p.q0 + q3*p.q1;
        a3 = q0*p.q3 + q1*p.q2 + q2*p.q1 + q3*p.q0;
        return {a0,a1,a2,a3};
    }
    quat scalar_product(quat p){
        quat q = *this;
        return ((p.conjugate()*q + q.conjugate()*p)/2);
    }
    quat vector_product(quat p){
        quat q = *this;
        return ((p*q - q*p)/2);
    }
    quat reciprocal(){
        quat q = *this;
        return (q.conjugate()/(q.norm()));
    }
    quat operator =(quat p){
        q0 = p.q0;
        q1 = p.q1;
        q2 = p.q2;
        q3 = p.q3;
        return *this;
    }
    quat operator =(double a){
        q0 = a;
        q1 = 0;
        q2 = 0;
        q3 = 0;
        return *this;
    }
};

quat operator +(double d, quat q){
    q.q0 += d;
    return q;
}

quat operator -(double d, quat q){
    q.q0 -= d;
    return q;
}

quat operator *(double d, quat q){
    q = q*d;
    return q;
}

std::ostream& operator << (std::ostream &os, const quat &q)
{
    return os  << q.q0 << " " << q.q1 << " " << q.q2 << " " << q.q3;
}

std::istream& operator >> (std::istream& in, quat& q)
{
    in >> q.q0 >> q.q1 >> q.q2 >> q.q3;
    return in;
}

std::ostream& operator << (std::ostream &os, const pure &q)
{
    return os  << q.x << " " << q.y << " " << q.z;
}

std::istream& operator >> (std::istream& in, pure& q)
{
    in >> q.x >> q.y >> q.z;
    return in;
}


