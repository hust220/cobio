#pragma once

#include "jnc_core.h"

namespace jnc {
namespace geom {

template <typename NumType> using MatX = Eigen::Matrix<NumType, -1, -1>;
using Matd = MatX<double>;

template <typename NumType> using VecX = Eigen::Matrix<NumType, -1, 1>;
using Vecd = VecX<double>;

template <typename T, typename U> void translate(T &&t, const U &u) {
  for (int i = 0; i < 3; i++)
    t[i] += u[i];
}

/**
 * Rotate fixed in the origin
 */
template <typename T, typename _Mat> void rotate(T &&t, const _Mat &mat) {
  auto x = t[0] * mat(0, 0) + t[1] * mat(1, 0) + t[2] * mat(2, 0);
  auto y = t[0] * mat(0, 1) + t[1] * mat(1, 1) + t[2] * mat(2, 1);
  auto z = t[0] * mat(0, 2) + t[1] * mat(1, 2) + t[2] * mat(2, 2);
  t[0] = x;
  t[1] = y;
  t[2] = z;
}

/**
 * Rotate fixed in an specified point
 */
template <typename T, typename U>
void rotate(T &&t, const U &origin, const MatX<double> &mat) {
  for (int i = 0; i < 3; i++)
    t[i] -= origin[i];
  rotate(t, mat);
  for (int i = 0; i < 3; i++)
    t[i] += origin[i];
}

template <class C1, class C2> MatX<double> x_rot_mat(C1 c, C2 s) {
  MatX<double> rot_mat(3, 3);
  rot_mat << 1, 0, 0, 0, c, s, 0, -s, c;
  return rot_mat;
}

template <class C1, class C2> MatX<double> y_rot_mat(C1 c, C2 s) {
  MatX<double> rot_mat(3, 3);
  rot_mat << c, 0, -s, 0, 1, 0, s, 0, c;
  return rot_mat;
}

template <class C1, class C2> MatX<double> z_rot_mat(C1 c, C2 s) {
  MatX<double> rot_mat(3, 3);
  rot_mat << c, s, 0, -s, c, 0, 0, 0, 1;
  return rot_mat;
}

template <typename T> MatX<double> rot_mat(int i, T &&v) {
  double c = std::cos(v);
  double s = std::sin(v);
  assert(i >= 0 && i < 3);
  return (i == 0 ? x_rot_mat(c, s)
                 : (i == 1 ? y_rot_mat(c, s) : z_rot_mat(c, s)));
}

// Rotate along with an axis
class RotateAlong {
public:
  using mat_t = MatX<double>;
  using vec_t = VecX<double>;

  mat_t m_rm;
  vec_t m_beg, m_end;

  RotateAlong() = default;

  template <typename Begin, typename End>
  RotateAlong(const Begin &begin, const End &end, double angle) {
    init(begin, end, angle);
  }

  template <typename Begin, typename End>
  void init(const Begin &begin, const End &end, double t) {
    double r1, r2, c1, c2, s1, s2;
    // L l = end - begin;
    vec_t l(3);
    for (int i = 0; i < 3; i++)
      l[i] = end[i] - begin[i];

    m_beg.resize(3);
    m_end.resize(3);
    for (int i = 0; i < 3; i++) {
      m_beg[i] = begin[i];
      m_end[i] = end[i];
    }

    r1 = std::sqrt(l[0] * l[0] + l[1] * l[1]);
    c1 = l[1] / r1;
    s1 = l[0] / r1;

    r2 = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    c2 = l[2] / r2;
    s2 = std::sqrt(l[0] * l[0] + l[1] * l[1]) / r2;

    m_rm = mat_t::Identity(3, 3);
    if (r1 != 0)
      m_rm = m_rm * z_rot_mat(c1, s1);
    if (r2 != 0)
      m_rm = m_rm * x_rot_mat(c2, s2);
    m_rm = m_rm * z_rot_mat(std::cos(t), std::sin(t));
    if (r2 != 0)
      m_rm = m_rm * x_rot_mat(c2, -s2);
    if (r1 != 0)
      m_rm = m_rm * z_rot_mat(c1, -s1);
  }

  template <typename T> RotateAlong &operator()(T &&t) {
    rotate(t, m_beg, m_rm);
    return *this;
  }
};

template <typename T> T square(T n) { return n * n; }

template <class P1, class P2> double distance(const P1 &p1, const P2 &p2) {
  return std::sqrt(square(p1[0] - p2[0]) + square(p1[1] - p2[1]) +
                   square(p1[2] - p2[2]));
}

template <class P1, class P2> double dist2(const P1 &p1, const P2 &p2) {
  return square(p1[0] - p2[0]) + square(p1[1] - p2[1]) + square(p1[2] - p2[2]);
}

template <typename P1, typename P2>
inline double norm(P1 &&p1, P2 &&p2, int n) {
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  return std::sqrt(sum);
}

template <class P1, class P2, class P3>
inline double angle(const P1 &p1, const P2 &p2, const P3 &p3) {
  double a1 = p1[0] - p2[0], a2 = p1[1] - p2[1], a3 = p1[2] - p2[2];
  double b1 = p3[0] - p2[0], b2 = p3[1] - p2[1], b3 = p3[2] - p2[2];
  double ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
  double rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
  return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb));
}

template <class T1, class T2, class T3, class T4>
inline double chirality(const T1 &p1, const T2 &p2, const T3 &p3,
                        const T4 &p4) {
  double a[3] = {p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]};
  double b[3] = {p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]};
  double c[3] = {p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]};
  double d[3] = {b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2],
                 b[0] * c[1] - b[1] * c[0]};
  return a[0] * d[0] + a[1] * d[1] + a[2] * d[2];
}

template <class P = Eigen::Vector3d, class P1, class P2, class P3>
P normal_vector(const P1 &p1, const P2 &p2, const P3 &p3) {
  double a1, a2, a3, b1, b2, b3, r;
  P p;
  a1 = p2[0] - p1[0];
  a2 = p2[1] - p1[1];
  a3 = p2[2] - p1[2];
  b1 = p3[0] - p2[0];
  b2 = p3[1] - p2[1];
  b3 = p3[2] - p2[2];
  p[0] = a2 * b3 - a3 * b2;
  p[1] = b1 * a3 - a1 * b3;
  p[2] = a1 * b2 - a2 * b1;
  r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  p[0] /= r;
  p[1] /= r;
  p[2] /= r;
  return p;
}

template <class P1, class P2, class P3, class P4>
inline double dihedral(const P1 &p1, const P2 &p2, const P3 &p3, const P4 &p4) {
  double a1, a2, a3, b1, b2, b3, c1, c2, c3;
  double a1_, a2_, a3_, b1_, b2_, b3_, c1_, c2_, c3_;
  double r, c, s;

  a1 = p1[0] - p2[0];
  a2 = p1[1] - p2[1];
  a3 = p1[2] - p2[2];
  b1 = p3[0] - p2[0];
  b2 = p3[1] - p2[1];
  b3 = p3[2] - p2[2];
  c1 = p4[0] - p2[0];
  c2 = p4[1] - p2[1];
  c3 = p4[2] - p2[2];
  if (b1 * b1 + b2 * b2 != 0) {
    r = sqrt(b1 * b1 + b2 * b2);
    c = b1 / r;
    s = -b2 / r;
    a1_ = c * a1 - s * a2;
    a2_ = s * a1 + c * a2;
    a3_ = a3;
    b1_ = c * b1 - s * b2;
    b2_ = s * b1 + c * b2;
    b3_ = b3;
    c1_ = c * c1 - s * c2;
    c2_ = s * c1 + c * c2;
    c3_ = c3;
  } else {
    a1_ = a1;
    a2_ = a2;
    a3_ = a3;
    b1_ = b1;
    b2_ = b2;
    b3_ = b3;
    c1_ = c1;
    c2_ = c2;
    c3_ = c3;
  }
  if (b1_ * b1_ + b3_ * b3_ != 0) {
    r = sqrt(b1_ * b1_ + b3_ * b3_);
    c = b1_ / r;
    s = b3_ / r;
    a1 = c * a1_ + s * a3_;
    a2 = a2_;
    a3 = -s * a1_ + c * a3_;
    b1 = c * b1_ + s * b3_;
    b2 = b2_;
    b3 = -s * b1_ + c * b3_;
    c1 = c * c1_ + s * c3_;
    c2 = c2_;
    c3 = -s * c1_ + c * c3_;
  } else {
    a1 = a1_;
    a2 = a2_;
    a3 = a3_;
    b1 = b1_;
    b2 = b2_;
    b3 = b3_;
    c1 = c1_;
    c2 = c2_;
    c3 = c3_;
  }
  if (a2 * a2 + a3 * a3 != 0) {
    r = sqrt(a2 * a2 + a3 * a3);
    c = a3 / r;
    s = a2 / r;
    a1_ = a1;
    a2_ = c * a2 - s * a3;
    a3_ = s * a2 + c * a3;
    b1_ = b1;
    b2_ = c * b2 - s * b3;
    b3_ = s * b2 + c * b3;
    c1_ = c1;
    c2_ = c * c2 - s * c3;
    c3_ = s * c2 + c * c3;
  }
  if (c2_ * c2_ + c3_ * c3_ == 0)
    return 0;
  double temp = acos(c3_ / sqrt(c2_ * c2_ + c3_ * c3_));
  if (c2_ > 0) {
    temp = -temp;
  }
  return temp;
}

template <typename NumType> class SupPos {
public:
  using mat_t = MatX<NumType>;
  using vec_t = VecX<NumType>;

  mat_t rot;
  vec_t c1;
  vec_t c2;
  NumType rmsd;

  SupPos() = default;

  SupPos(const mat_t &m, const mat_t &n) { init(m, n); }

  void init(const mat_t &m, const mat_t &n) {
    int i, j, len;
    mat_t x, y, g, u, v, I, d;
    std::ostringstream stream;

    if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
      stream << "jian::geom::suppos error! (" << m.rows() << ' ' << m.cols()
             << ") (" << n.rows() << ' ' << n.cols() << ")\n";
      throw stream.str();
    }

    len = m.rows();
    x = m;
    y = n;
    c1 = vec_t::Zero(3);
    c2 = vec_t::Zero(3);
    for (i = 0; i < len; i++)
      for (j = 0; j < 3; j++) {
        c1[j] += x(i, j);
        c2[j] += y(i, j);
      }
    for (i = 0; i < 3; i++) {
      c1[i] = c1[i] / len;
      c2[i] = c2[i] / len;
    }
    for (i = 0; i < len; i++)
      for (j = 0; j < 3; j++) {
        x(i, j) -= c1[j];
        y(i, j) -= c2[j];
      }

    g = x.transpose() * y;
    Eigen::JacobiSVD<mat_t> svd(g, Eigen::ComputeFullU | Eigen::ComputeFullV);
    u = svd.matrixU();
    v = svd.matrixV();

    I = mat_t::Identity(3, 3);
    if (g.determinant() < 0)
      I(2, 2) = -1;
    rot = u * I * v.transpose();

    d = x * rot - y;
    rmsd = 0;
    for (int i = 0; i < len; i++)
      for (int j = 0; j < 3; j++)
        rmsd += d(i, j) * d(i, j);
    rmsd = std::sqrt(rmsd / len);

    c1 = -c1;
  }

  template <typename T> void apply(T &&t) {
    translate(t, c1);
    rotate(t, rot);
    translate(t, c2);
  }

  template <typename M> void apply_m(M &m) {
    int i, j;
    for (i = 0; i < m.rows(); i++) {
      for (j = 0; j < 3; j++) {
        m(i, j) += c1[j];
      }
    }
    m = m * rot;
    for (i = 0; i < m.rows(); i++) {
      for (j = 0; j < 3; j++) {
        m(i, j) += c2[j];
      }
    }
  }
};

template <typename NumType>
SupPos<NumType> suppos(const MatX<NumType> &m, const MatX<NumType> &n) {
  return SupPos<NumType>(m, n);
}

template <typename T, typename NumType>
SupPos<NumType> suppos(T &t, const MatX<NumType> &m, const MatX<NumType> &n) {
  SupPos<NumType> sp(m, n);
  sp.apply_m(t);
  return sp;
}

template <typename NumType, typename U, typename F, typename V>
SupPos<NumType> suppos(MatX<NumType> &src, const U &src_indices, const F &tgt,
                       const V &tgt_indices) {
  MatX<NumType> m(src_indices.size(), 3), n(tgt_indices.size(), 3);
  for (int i = 0; i < src_indices.size(); i++)
    for (int j = 0; j < 3; j++) {
      m(i, j) = src(src_indices[i], j);
      n(i, j) = tgt(tgt_indices[i], j);
    }

  auto sp = suppos(m, n);

  sp.apply_m(src);
  return sp;
}

template <typename NumType1, typename NumType2>
NumType1 rmsd(const MatX<NumType1> &x, const MatX<NumType2> &y) {
  return SupPos<NumType1>(x, y).rmsd;
}

template <typename NumType, typename A, typename B, typename C, typename D>
MatX<NumType> suppos_axis_polar(A theta_o, B phi_o, C theta_n, D phi_n) {
  MatX<NumType> m = MatX<NumType>::Identity(3, 3);
  double ang = -phi_o;
  if (ang != 0)
    m *= z_rot_mat(std::cos(ang), std::sin(ang));
  ang = -theta_o + theta_n;
  if (ang != 0)
    m *= y_rot_mat(std::cos(ang), std::sin(ang));
  ang = phi_n;
  if (ang != 0)
    m *= z_rot_mat(std::cos(ang), std::sin(ang));
  return m;
}

template <typename NumType, typename O, typename N>
MatX<NumType> suppos_axis_xyz(const O &o, const N &n) {
  MatX<NumType> m = MatX<NumType>::Identity(3, 3);
  double r, x, y, r1, x1, y1, r2, x2, y2, c, s, c1, c2, s1, s2;
  r = std::sqrt(o[0] * o[0] + o[1] * o[1]);
  x = o[0];
  y = o[1];
  if (r != 0) {
    c = y / r;
    s = x / r;
    if (s != 0)
      m *= z_rot_mat(c, s);
  }
  r1 = std::sqrt(o[0] * o[0] + o[1] * o[1] + o[2] * o[2]);
  x1 = std::sqrt(o[0] * o[0] + o[1] * o[1]);
  y1 = o[2];
  r2 = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  x2 = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  y2 = n[2];
  if (r1 != 0 && r2 != 0) {
    c1 = y1 / r1;
    s1 = x1 / r1;
    c2 = y2 / r2;
    s2 = x2 / r2;
    c = c1 * c2 + s1 * s2;
    s = s1 * c2 - s2 * c1;
    if (s != 0)
      m *= x_rot_mat(c, s);
  }
  r = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  x = n[0];
  y = n[1];
  if (r != 0) {
    c = y / r;
    s = -x / r;
    if (s != 0)
      m *= z_rot_mat(c, s);
  }
  return m;
}

} // namespace geom
} // namespace jnc

