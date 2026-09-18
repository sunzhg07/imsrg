///////////////////////////////////////////////////////////////////////////////////
// C++-only tensor leftover lock: m-unfactorized, m-factorized, J-scheme.
// No Python expansion of MEs / CGs / packaging scans.
//
//   ./build/test_tensor_factorized
//   ./build/test_tensor_factorized 2 2 8 1
//   ./build/test_tensor_factorized emax=2 lambda=2 max_m_cmp=8 diagrams=1
//   ./build/test_tensor_factorized emax=2 lambda=2 step=1 diagrams=0
// step: 0=all, 1=m-unfact≡J-unfact, 2=m-unfact≡m-fact, 3=diagrams only,
//       4=J-unfact≡J-fact
///////////////////////////////////////////////////////////////////////////////////

#include "ModelSpace.hh"
#include "UnitTest.hh"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <string>

static void usage(const char *prog)
{
  std::cerr
      << "Usage:\n"
      << "  " << prog << " [emax] [lambda] [max_m_cmp] [diagrams]\n"
      << "  " << prog << " emax=2 lambda=2 T=0 max_m_cmp=8 diagrams=1 step=1\n"
      << "Defaults: emax=1 lambda=2 T=0 max_m_cmp=8 diagrams=1 step=0\n"
      << "step: 0=all 1=m-unfact≡J-nested 2=m-unfact≡m-fact 3=diagrams 4=J-nested≡J-fact\n";
}

static bool starts_with(const char *s, const char *pre)
{
  return std::strncmp(s, pre, std::strlen(pre)) == 0;
}

static int parse_kv_int(const char *arg, const char *key)
{
  const size_t n = std::strlen(key);
  return std::atoi(arg + n);
}

int main(int argc, char **argv)
{
  int emax = 1;
  int lam = 2;
  int trank = 0;
  int max_m = 8;
  int diagrams_i = 1;
  int step = 0;
  bool diagrams_set = false;
  const uint64_t seed = 17;
  int positional = 0;

  for (int i = 1; i < argc; ++i)
  {
    const char *a = argv[i];
    if (std::strcmp(a, "-h") == 0 or std::strcmp(a, "--help") == 0)
    {
      usage(argv[0]);
      return 0;
    }
    if (starts_with(a, "emax="))
      emax = parse_kv_int(a, "emax=");
    else if (starts_with(a, "lambda="))
      lam = parse_kv_int(a, "lambda=");
    else if (starts_with(a, "lam="))
      lam = parse_kv_int(a, "lam=");
    else if (starts_with(a, "T="))
      trank = parse_kv_int(a, "T=");
    else if (starts_with(a, "trank="))
      trank = parse_kv_int(a, "trank=");
    else if (starts_with(a, "max_m_cmp="))
      max_m = parse_kv_int(a, "max_m_cmp=");
    else if (starts_with(a, "diagrams="))
    {
      diagrams_i = parse_kv_int(a, "diagrams=");
      diagrams_set = true;
    }
    else if (starts_with(a, "step="))
      step = parse_kv_int(a, "step=");
    else if (a[0] == '-' and (a[1] < '0' or a[1] > '9'))
    {
      std::cerr << "Unknown option: " << a << std::endl;
      usage(argv[0]);
      return 1;
    }
    else
    {
      const int v = std::atoi(a);
      if (positional == 0)
        emax = v;
      else if (positional == 1)
        lam = v;
      else if (positional == 2)
        max_m = v;
      else if (positional == 3)
        diagrams_i = v;
      else
      {
        std::cerr << "Unexpected extra argument: " << a << std::endl;
        usage(argv[0]);
        return 1;
      }
      ++positional;
    }
  }

  if (not diagrams_set and step != 0 and step != 3)
    diagrams_i = 0;
  const bool diagrams = (step == 3) or (diagrams_i != 0 and step != 1 and step != 2 and step != 4);
  const int threeway_step = (step == 3) ? -1 : step;

  std::cout << "test_tensor_factorized  emax=" << emax << "  λ=" << lam
            << "  T=" << trank << "  max_m_cmp=" << max_m << "  diagrams=" << diagrams
            << "  step=" << step << "  seed=" << seed << std::endl;

  // GetSixJ is not thread-safe on first fill; gold benches run single-threaded.
  omp_set_num_threads(1);

  ModelSpace ms(emax, "He4", "He4");
  ms.SetHbarOmega(20.0);
  ms.PreCalculateSixJ();
  ms.PreCalculateNineJ();

  UnitTest ut(ms);
  ut.SetRandomSeed(seed);

  bool ok = true;
  if (threeway_step >= 0)
    ok = ut.TestTensorFactorizedThreeway(lam, max_m, threeway_step, trank);
  if (diagrams)
    ok = ut.TestTensorFactorizedDiagrams(lam, std::min(max_m, 4)) and ok;

  std::cout << "\ntest_tensor_factorized  OVERALL: "
            << (ok ? "PASS" : "FAIL") << std::endl;
  return ok ? 0 : 1;
}
