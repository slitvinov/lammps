#ifndef LMP_MATH_EIGEN_IMPL_H
#define LMP_MATH_EIGEN_IMPL_H 
#include <numeric>
#include <complex>
#include <limits>
#include <cmath>
#include <vector>
#include <random>
#include <functional>
namespace MathEigen {
  template<typename Entry>
  void Alloc2D(size_t nrows,
               size_t ncols,
               Entry ***paaX);
  template<typename Entry>
  void Dealloc2D(Entry ***paaX);
  template <typename T>
  struct realTypeMap {
    typedef T type;
  };
  template <typename T>
  struct realTypeMap<std::complex<T>> {
    typedef T type;
  };
  template <typename T>
  using real_t = typename realTypeMap<T>::type;
  template <typename T>
  T inner_prod(const std::vector<T>& v1, const std::vector<T>& v2);
  template <typename T>
  real_t<T> l1_norm(const std::vector<T>& v);
  template <typename T>
  real_t<T> l2_norm(const std::vector<T>& v);
  template <typename T1, typename T2>
  void scalar_mul(T1 c, std::vector<T2>& v);
  template <typename T>
  void normalize(std::vector<T>& v);
  template<typename Scalar,
           typename Vector,
           typename Matrix,
           typename ConstMatrix=Matrix>
  class Jacobi
  {
    int n;
    Scalar **M;
    Scalar c;
    Scalar s;
    Scalar t;
    int *max_idx_row;
  public:
    Jacobi(int n);
    ~Jacobi();
    void SetSize(int n);
    typedef enum eSortCriteria {
      DO_NOT_SORT,
      SORT_DECREASING_EVALS,
      SORT_INCREASING_EVALS,
      SORT_DECREASING_ABS_EVALS,
      SORT_INCREASING_ABS_EVALS
    } SortCriteria;
    int
    Diagonalize(ConstMatrix mat,
                Vector eval,
                Matrix evec,
                SortCriteria sort_criteria=SORT_DECREASING_EVALS,
                bool calc_evecs=true,
                int max_num_sweeps=50
                );
    Jacobi(int n, Scalar **M, int *max_idx_row);
  private:
    bool is_preallocated;
    void _Jacobi(int n, Scalar **M, int *max_idx_row);
    void CalcRot(Scalar const *const *M, int i, int j);
    void ApplyRot(Scalar **M, int i, int j);
    void ApplyRotLeft(Matrix E, int i, int j);
    int MaxEntryRow(Scalar const *const *M, int i) const;
    void MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const;
    void SortRows(Vector v, Matrix M, int n, SortCriteria s=SORT_DECREASING_EVALS) const;
    void Init();
    void Alloc(int n);
    void Dealloc();
  public:
    Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source);
    Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other);
    void swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other);
    Jacobi<Scalar, Vector, Matrix, ConstMatrix>& operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source);
  };
  template <typename T>
  struct VectorRandomInitializer {
  public:
    static void init(std::vector<T>&);
  };
  template <typename T>
  struct VectorRandomInitializer<std::complex<T>> {
  public:
    static void init(std::vector<std::complex<T>>&);
  };
  template <typename T>
  inline constexpr int sig_decimal_digit() {
    return (int)(std::numeric_limits<T>::digits *
                 std::log10(std::numeric_limits<T>::radix));
  }
  template <typename T>
  inline constexpr T minimum_effective_decimal() {
    return std::pow(10, -sig_decimal_digit<T>());
  }
  template <typename T>
  class LambdaLanczos {
  public:
    LambdaLanczos();
    LambdaLanczos(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul, int matrix_size, bool find_maximum);
    LambdaLanczos(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul, int matrix_size) : LambdaLanczos(mv_mul, matrix_size, true) {}
    int run(real_t<T>&, std::vector<T>&) const;
    int matrix_size;
    std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul;
    bool find_maximum = false;
    real_t<T> eigenvalue_offset = 0.0;
    void ChooseOffset();
    int max_iteration;
    real_t<T> eps = minimum_effective_decimal<real_t<T>>() * 1e3;
    real_t<T> tridiag_eps_ratio = 1e-1;
    int initial_vector_size = 200;
    std::function<void(std::vector<T>&)> init_vector =
      VectorRandomInitializer<T>::init;
    int SetSize(int matrix_size);
    void SetMul(std::function<void(const std::vector<T>&,
                                   std::vector<T>&)> mv_mul);
    void SetInitVec(std::function<void(std::vector<T>&)> init_vector);
    void SetFindMax(bool find_maximum);
    void SetEvalOffset(T eigenvalue_offset);
    void SetEpsilon(T eps);
    void SetTriEpsRatio(T tridiag_eps_ratio);
  private:
    static void schmidt_orth(std::vector<T>&, const std::vector<std::vector<T>>&);
    real_t<T> find_minimum_eigenvalue(const std::vector<real_t<T>>&,
                                      const std::vector<real_t<T>>&) const;
    real_t<T> find_maximum_eigenvalue(const std::vector<real_t<T>>&,
                                      const std::vector<real_t<T>>&) const;
    static real_t<T> tridiagonal_eigen_limit(const std::vector<real_t<T>>&,
                                             const std::vector<real_t<T>>&);
    static int num_of_eigs_smaller_than(real_t<T>,
                                        const std::vector<real_t<T>>&,
                                        const std::vector<real_t<T>>&);
    real_t<T> UpperBoundEvals() const;
  };
  template<typename Scalar, typename Vector, typename ConstMatrix>
  class PEigenDense
  {
    size_t n;
    std::vector<Scalar> evec;
  public:
    PEigenDense(int matrix_size=0);
    real_t<Scalar>
    PrincipalEigen(ConstMatrix matrix,
                   Vector evector,
                   bool find_max=false);
    void SetSize(int matrix_size);
  };
template<typename Entry>
void Alloc2D(size_t nrows,
             size_t ncols,
             Entry ***paaX)
{
  *paaX = new Entry* [nrows];
  (*paaX)[0] = new Entry [nrows * ncols];
  for (size_t iy=0; iy<nrows; iy++)
    (*paaX)[iy] = (*paaX)[0] + iy*ncols;
}
template<typename Entry>
void Dealloc2D(Entry ***paaX)
{
  if (paaX && *paaX) {
    delete [] (*paaX)[0];
    delete [] (*paaX);
    *paaX = nullptr;
  }
}
template <typename T>
struct ConjugateProduct {
public:
  static T prod(T a, T b) { return a*b; }
};
template <typename T>
struct ConjugateProduct<std::complex<T>> {
public:
  static std::complex<T> prod(std::complex<T> a, std::complex<T> b) {
    return std::conj(a)*b;
  }
};
template <typename T>
inline T inner_prod(const std::vector<T>& v1, const std::vector<T>& v2) {
  return std::inner_product(std::begin(v1), std::end(v1),
                            std::begin(v2), T(),
                            [](T a, T b) -> T { return a+b; },
                            ConjugateProduct<T>::prod);
}
template <typename T>
inline real_t<T> l2_norm(const std::vector<T>& vec) {
  return std::sqrt(std::real(inner_prod(vec, vec)));
}
template <typename T1, typename T2>
inline void scalar_mul(T1 a, std::vector<T2>& vec) {
  int n = vec.size();
  for (int i = 0;i < n;i++)
    vec[i] *= a;
}
template <typename T>
inline void normalize(std::vector<T>& vec) {
  scalar_mul(1.0/l2_norm(vec), vec);
}
template <typename T>
inline real_t<T> l1_norm(const std::vector<T>& vec) {
  real_t<T> norm = real_t<T>();
  for (const T& element : vec)
    norm += std::abs(element);
  return norm;
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(int n) {
  _Jacobi(n, nullptr, nullptr);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(int n, Scalar **M, int *max_idx_row) {
  _Jacobi(n, M, max_idx_row);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
_Jacobi(int n, Scalar **M, int *max_idx_row) {
  Init();
  if (M) {
    is_preallocated = true;
    this->n = n;
    this->M = M;
    this->max_idx_row = max_idx_row;
  } else {
    is_preallocated = false;
    SetSize(n);
  }
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
~Jacobi() {
  if (! is_preallocated)
    Dealloc();
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Diagonalize(ConstMatrix mat,
            Vector eval,
            Matrix evec,
            SortCriteria sort_criteria,
            bool calc_evec,
            int max_num_sweeps)
{
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      M[i][j] = mat[i][j];
  if (calc_evec)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        evec[i][j] = (i==j) ? 1.0 : 0.0;
  for (int i = 0; i < n-1; i++)
    max_idx_row[i] = MaxEntryRow(M, i);
  int n_iters;
  int max_num_iters = max_num_sweeps*n*(n-1)/2;
  for (n_iters=0; n_iters < max_num_iters; n_iters++) {
    int i,j;
    MaxEntry(M, i, j);
    if ((M[i][i] + M[i][j] == M[i][i]) && (M[j][j] + M[i][j] == M[j][j])) {
      M[i][j] = 0.0;
      max_idx_row[i] = MaxEntryRow(M,i);
    }
    if (M[i][j] == 0.0)
      break;
    CalcRot(M, i, j);
    ApplyRot(M, i, j);
    if (calc_evec)
      ApplyRotLeft(evec,i,j);
  }
  for (int i = 0; i < n; i++)
    eval[i] = M[i][i];
  SortRows(eval, evec, n, sort_criteria);
  return (n_iters == max_num_iters);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
CalcRot(Scalar const *const *M,
        int i,
        int j)
{
  t = 1.0;
  Scalar M_jj_ii = (M[j][j] - M[i][i]);
  if (M_jj_ii != 0.0) {
    Scalar kappa = M_jj_ii;
    t = 0.0;
    Scalar M_ij = M[i][j];
    if (M_ij != 0.0) {
      kappa /= (2.0*M_ij);
      t = 1.0 / (std::sqrt(1 + kappa*kappa) + std::abs(kappa));
      if (kappa < 0.0)
        t = -t;
    }
  }
  c = 1.0 / std::sqrt(1 + t*t);
  s = c*t;
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRot(Scalar **M,
         int i,
         int j)
{
  M[i][i] -= t * M[i][j];
  M[j][j] += t * M[i][j];
  M[i][j] = 0.0;
  for (int w=0; w < i; w++) {
    M[i][w] = M[w][i];
    M[w][i] = c*M[w][i] - s*M[w][j];
    if (i == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][i])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=i;
  }
  for (int w=i+1; w < j; w++) {
    M[w][i] = M[i][w];
    M[i][w] = c*M[i][w] - s*M[w][j];
  }
  for (int w=j+1; w < n; w++) {
    M[w][i] = M[i][w];
    M[i][w] = c*M[i][w] - s*M[j][w];
  }
  max_idx_row[i] = MaxEntryRow(M, i);
  for (int w=0; w < i; w++) {
    M[w][j] = s*M[i][w] + c*M[w][j];
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
  }
  for (int w=i+1; w < j; w++) {
    M[w][j] = s*M[w][i] + c*M[w][j];
    if (j == max_idx_row[w]) max_idx_row[w] = MaxEntryRow(M, w);
    else if (std::abs(M[w][j])>std::abs(M[w][max_idx_row[w]])) max_idx_row[w]=j;
  }
  for (int w=j+1; w < n; w++) {
    M[j][w] = s*M[w][i] + c*M[j][w];
  }
  max_idx_row[j] = MaxEntryRow(M, j);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
ApplyRotLeft(Matrix E,
             int i,
             int j)
{
  for (int v = 0; v < n; v++) {
    Scalar Eiv = E[i][v];
    E[i][v] = c*E[i][v] - s*E[j][v];
    E[j][v] = s*Eiv + c*E[j][v];
  }
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
int Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntryRow(Scalar const *const *M, int i) const {
  int j_max = i+1;
  for (int j = i+2; j < n; j++)
    if (std::abs(M[i][j]) > std::abs(M[i][j_max]))
      j_max = j;
  return j_max;
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
MaxEntry(Scalar const *const *M, int& i_max, int& j_max) const {
  i_max = 0;
  j_max = max_idx_row[i_max];
  Scalar max_entry = std::abs(M[i_max][j_max]);
  int nm1 = n-1;
  for (int i=1; i < nm1; i++) {
    int j = max_idx_row[i];
    if (std::abs(M[i][j]) > max_entry) {
      max_entry = std::abs(M[i][j]);
      i_max = i;
      j_max = j;
    }
  }
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SortRows(Vector eval,
         Matrix evec,
         int n,
         SortCriteria sort_criteria) const
{
  for (int i = 0; i < n-1; i++) {
    int i_max = i;
    for (int j = i+1; j < n; j++) {
      switch (sort_criteria) {
      case SORT_DECREASING_EVALS:
        if (eval[j] > eval[i_max])
          i_max = j;
        break;
      case SORT_INCREASING_EVALS:
        if (eval[j] < eval[i_max])
          i_max = j;
        break;
      case SORT_DECREASING_ABS_EVALS:
        if (std::abs(eval[j]) > std::abs(eval[i_max]))
          i_max = j;
        break;
      case SORT_INCREASING_ABS_EVALS:
        if (std::abs(eval[j]) < std::abs(eval[i_max]))
          i_max = j;
        break;
      default:
        break;
      }
    }
    std::swap(eval[i], eval[i_max]);
    for (int k = 0; k < n; k++)
      std::swap(evec[i][k], evec[i_max][k]);
  }
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Init() {
  n = 0;
  M = nullptr;
  max_idx_row = nullptr;
  is_preallocated = false;
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
SetSize(int n) {
  Dealloc();
  Alloc(n);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Alloc(int n) {
  this->n = n;
  if (n > 0) {
    max_idx_row = new int[n];
    Alloc2D(n, n, &M);
  }
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Dealloc() {
  Dealloc2D(&M);
  delete[] max_idx_row;
  max_idx_row = nullptr;
  Init();
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(const Jacobi<Scalar, Vector, Matrix, ConstMatrix>& source)
{
  Init();
  SetSize(source.n);
  std::copy(source.max_idx_row,
            source.max_idx_row + n,
            max_idx_row);
  for (int i = 0; i < n; i++)
    std::copy(source.M[i],
              source.M[i] + n,
              M[i]);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
void Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
swap(Jacobi<Scalar, Vector, Matrix, ConstMatrix> &other) {
  std::swap(n, other.n);
  std::swap(is_preallocated, other.is_preallocated);
  std::swap(max_idx_row, other.max_idx_row);
  std::swap(M, other.M);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
Jacobi(Jacobi<Scalar, Vector, Matrix, ConstMatrix>&& other) {
  Init();
  this->swap(other);
}
template<typename Scalar,typename Vector,typename Matrix,typename ConstMatrix>
Jacobi<Scalar, Vector, Matrix, ConstMatrix>&
Jacobi<Scalar, Vector, Matrix, ConstMatrix>::
operator = (Jacobi<Scalar, Vector, Matrix, ConstMatrix> source) {
  this->swap(source);
  return *this;
}
template <typename T>
inline LambdaLanczos<T>::LambdaLanczos() {
  this->matrix_size = 0;
  this->max_iteration = 0;
  this->find_maximum = false;
}
template <typename T>
inline LambdaLanczos<T>::
LambdaLanczos(std::function<void(const std::vector<T>&,
                                 std::vector<T>&)> mv_mul,
                                 int matrix_size,
              bool find_maximum)
{
  this->mv_mul = mv_mul;
  this->matrix_size = matrix_size;
  this->max_iteration = matrix_size;
  this->find_maximum = find_maximum;
}
template <typename T>
inline int LambdaLanczos<T>::
run(real_t<T>& eigvalue, std::vector<T>& eigvec) const
{
  std::vector<std::vector<T>> u;
  std::vector<real_t<T>> alpha;
  std::vector<real_t<T>> beta;
  const int n = this->matrix_size;
  u.reserve(this->initial_vector_size);
  alpha.reserve(this->initial_vector_size);
  beta.reserve(this->initial_vector_size);
  u.emplace_back(n, 0.0);
  std::vector<T> vk(n, 0.0);
  real_t<T> alphak = 0.0;
  alpha.push_back(alphak);
  real_t<T> betak = 0.0;
  beta.push_back(betak);
  std::vector<T> uk(n);
  this->init_vector(uk);
  normalize(uk);
  u.push_back(uk);
  real_t<T> ev, pev;
  pev = std::numeric_limits<real_t<T>>::max();
  int itern = this->max_iteration;
  for (int k = 1;k <= this->max_iteration;k++) {
    for (int i = 0;i < n;i++) {
      vk[i] = uk[i]*this->eigenvalue_offset;
    }
    this->mv_mul(uk, vk);
    alphak = std::real(inner_prod(u.back(), vk));
    alpha.push_back(alphak);
    for (int i = 0;i < n; i++) {
      uk[i] = vk[i] - betak*u[k-1][i] - alphak*u[k][i];
    }
    schmidt_orth(uk, u);
    betak = l2_norm(uk);
    beta.push_back(betak);
    if (this->find_maximum) {
      ev = find_maximum_eigenvalue(alpha, beta);
    } else {
      ev = find_minimum_eigenvalue(alpha, beta);
    }
    const real_t<T> zero_threshold = minimum_effective_decimal<real_t<T>>()*1e-1;
    if (betak < zero_threshold) {
      u.push_back(uk);
      itern = k;
      break;
    }
    normalize(uk);
    u.push_back(uk);
    if (abs(ev-pev) < std::min(abs(ev), abs(pev))*this->eps) {
      itern = k;
      break;
    } else {
      pev = ev;
    }
  }
  eigvalue = ev - this->eigenvalue_offset;
  int m = alpha.size();
  std::vector<T> cv(m+1);
  cv[0] = 0.0;
  cv[m] = 0.0;
  cv[m-1] = 1.0;
  beta[m-1] = 0.0;
  if (eigvec.size() < n) {
    eigvec.resize(n);
  }
  for (int i = 0;i < n;i++) {
    eigvec[i] = cv[m-1]*u[m-1][i];
  }
  for (int k = m-2;k >= 1;k--) {
    cv[k] = ((ev - alpha[k+1])*cv[k+1] - beta[k+1]*cv[k+2])/beta[k];
    for (int i = 0;i < n;i++) {
      eigvec[i] += cv[k]*u[k][i];
    }
  }
  normalize(eigvec);
  return itern;
}
template <typename T>
inline void LambdaLanczos<T>::
schmidt_orth(std::vector<T>& uorth, const std::vector<std::vector<T>>& u)
{
  int n = uorth.size();
  for (int k = 0;k < u.size();k++) {
    T innprod = inner_prod(uorth, u[k]);
    for (int i = 0;i < n;i++)
      uorth[i] -= innprod * u[k][i];
  }
}
template <typename T>
inline real_t<T> LambdaLanczos<T>::
find_minimum_eigenvalue(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta) const
{
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid;
  while (upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    if (nmid >= 1) {
      upper = mid;
    } else {
      lower = mid;
    }
    if (mid == pmid) {
      break;
    }
    pmid = mid;
  }
  return lower;
}
template <typename T>
inline real_t<T> LambdaLanczos<T>::
find_maximum_eigenvalue(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta) const
{
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid;
  int m = alpha.size() - 1;
  while (upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    if (nmid < m) {
      lower = mid;
    } else {
      upper = mid;
    }
    if (mid == pmid) {
      break;
    }
    pmid = mid;
  }
  return lower;
}
template <typename T>
inline real_t<T> LambdaLanczos<T>::
tridiagonal_eigen_limit(const std::vector<real_t<T>>& alpha,
                        const std::vector<real_t<T>>& beta)
{
  real_t<T> r = l1_norm(alpha);
  r += 2*l1_norm(beta);
  return r;
}
template <typename T>
inline int LambdaLanczos<T>::
num_of_eigs_smaller_than(real_t<T> c,
                         const std::vector<real_t<T>>& alpha,
                         const std::vector<real_t<T>>& beta)
{
  real_t<T> q_i = 1.0;
  int count = 0;
  int m = alpha.size();
  for (int i = 1;i < m;i++){
    q_i = alpha[i] - c - beta[i-1]*beta[i-1]/q_i;
    if (q_i < 0){
      count++;
    }
    if (q_i == 0){
      q_i = minimum_effective_decimal<real_t<T>>();
    }
  }
  return count;
}
template <typename T>
inline void LambdaLanczos<T>::ChooseOffset() {
  const auto n = this->matrix_size;
  std::vector<T> unit_vec_j(n);
  std::vector<T> matrix_column_j(n);
  real_t<T> eval_upper_bound = 0.0;
  for (int j = 0; j < n; j++) {
    std::fill(unit_vec_j.begin(), unit_vec_j.end(), 0);
    unit_vec_j[j] = 1.0;
    this->mv_mul(unit_vec_j, matrix_column_j);
    real_t<T> sum_column = 0.0;
    for (int i = 0; i < n; i++)
      sum_column += std::abs(matrix_column_j[i]);
    if (eval_upper_bound < sum_column)
      eval_upper_bound = sum_column;
  }
  if (find_maximum)
    this->eigenvalue_offset = eval_upper_bound;
  else
    this->eigenvalue_offset = -eval_upper_bound;
}
template <typename T>
inline int LambdaLanczos<T>::SetSize(int matrix_size)
{
  this->matrix_size = matrix_size;
  this->max_iteration = matrix_size;
  return matrix_size;
}
template <typename T>
inline void LambdaLanczos<T>::SetMul(std::function<void(const std::vector<T>&, std::vector<T>&)> mv_mul)
{
  this->mv_mul = mv_mul;
}
template <typename T>
inline void LambdaLanczos<T>::SetInitVec(std::function<void(std::vector<T>&)> init_vector)
{
  this->init_vector = init_vector;
}
template <typename T>
inline void LambdaLanczos<T>::SetFindMax(bool find_maximum) {
  this->find_maximum = find_maximum;
}
template <typename T>
inline void LambdaLanczos<T>::SetEvalOffset(T offset)
{
  this->eigenvalue_offset = offset;
}
template <typename T>
inline void LambdaLanczos<T>::SetEpsilon(T epsilon)
{
  this->eps = epsilon;
}
template <typename T>
inline void LambdaLanczos<T>::SetTriEpsRatio(T tri_eps_ratio)
{
  this->tridiag_eps_ratio = tri_eps_ratio;
}
template <typename T>
inline void VectorRandomInitializer<T>::
init(std::vector<T>& v)
{
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));
  int n = v.size();
  for (int i = 0;i < n;i++) {
    v[i] = rand(mt);
  }
  normalize(v);
}
template <typename T>
inline void VectorRandomInitializer<std::complex<T>>::
init(std::vector<std::complex<T>>& v)
{
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));
  int n = v.size();
  for (int i = 0;i < n;i++) {
    v[i] = std::complex<T>(rand(mt), rand(mt));
  }
  normalize(v);
}
template<typename Scalar, typename Vector, typename ConstMatrix>
void PEigenDense<Scalar, Vector, ConstMatrix>::
SetSize(int matrix_size) {
  n = matrix_size;
  evec.resize(n);
}
template<typename Scalar, typename Vector, typename ConstMatrix>
PEigenDense<Scalar, Vector, ConstMatrix>::
PEigenDense(int matrix_size):evec(matrix_size) {
  SetSize(matrix_size);
}
template<typename Scalar, typename Vector, typename ConstMatrix>
real_t<Scalar> PEigenDense<Scalar, Vector, ConstMatrix>::
PrincipalEigen(ConstMatrix matrix,
               Vector eigenvector,
               bool find_max)
{
  auto matmul = [&](const std::vector<Scalar>& in, std::vector<Scalar>& out) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        out[i] += matrix[i][j]*in[j];
      }
    }
  };
  auto init_vec = [&](std::vector<Scalar>& vec) {
    for (int i = 0; i < n; i++)
      vec[i] = 0.0;
    vec[0] = 1.0;
  };
  LambdaLanczos<Scalar> ll_engine(matmul, n, find_max);
  Scalar eval_upper_bound = 0.0;
  for (int i = 0; i < n; i++) {
    Scalar sum_row = 0.0;
    for (int j = 0; j < n; i++)
      sum_row += std::abs(matrix[i][j]);
    if (eval_upper_bound < sum_row)
      eval_upper_bound = sum_row;
  }
  if (find_max)
    ll_engine.eigenvalue_offset = eval_upper_bound;
  else
    ll_engine.eigenvalue_offset = -eval_upper_bound;
  ll_engine.init_vector = init_vec;
  Scalar eval;
  ll_engine.run(eval, evec);
  for (int i = 0; i < n; i++)
    eigenvector[i] = evec[i];
  return eval;
}
}
#endif
