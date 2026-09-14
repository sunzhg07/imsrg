#ifndef evc_hh
#define evc_hh 1

#include "Operator.hh"

#include <array>
#include <cstddef>

class EVC
{
  public:
    struct ClusterAmplitudes
    {
      double z0 = 0.0;
      Operator z;

      ClusterAmplitudes(ModelSpace &ms)
          : z(ms, 0, 0, 0, 2)
      {
        z.SetNonHermitian();
      }
    };

    struct RHS
    {
      double dz0 = 0.0;
      Operator dz;
      Operator transformed;

      RHS(ModelSpace &ms)
          : dz(ms, 0, 0, 0, 2), transformed(ms, 0, 0, 0, 2)
      {
        dz.SetNonHermitian();
        transformed.SetNonHermitian();
      }
    };

    struct RHSDiagrams
    {
      double dz0 = 0.0;
      Operator c1;
      Operator c2;
      Operator c3;
      Operator c4;

      RHSDiagrams(ModelSpace &ms)
          : c1(ms, 0, 0, 0, 2), c2(ms, 0, 0, 0, 2), c3(ms, 0, 0, 0, 2), c4(ms, 0, 0, 0, 2)
      {
        c1.SetNonHermitian();
        c2.SetNonHermitian();
        c3.SetNonHermitian();
        c4.SetNonHermitian();
      }
    };

    explicit EVC(ModelSpace &ms);

    void SetUseIMSRG3Intermediates(bool tf) { use_imsrg3_intermediates_ = tf; }
    void SetEulerSteps(std::size_t nsteps) { euler_steps_ = nsteps; }
    void SetUseRK4(bool tf) { use_rk4_ = tf; }

    Operator ExtractExcitationPart(const Operator &op) const;
    Operator ExtractDeexcitationPart(const Operator &op) const;
    Operator ExtractTFromOmega(const Operator &omega) const;
    Operator ExtractTDaggerFromOmega(const Operator &omega) const;

    Operator BuildClusterOperator(const ClusterAmplitudes &amps) const;
    // J-scheme CC Z1 residual d z_i^a / dλ = -R_i^a, stored in non-Hermitian OneBody(a,i).
    // AMC: run/evc/amc_eqs/output/z1_terms.tex . m-scheme gold: ReferenceImplementations::evc_z1_mscheme.
    Operator BuildZ1RHS(const Operator &t, const Operator &z) const;
    // J-scheme CC Z2 residual d z_{ij}^{ab}/dλ = -R. AMC: run/evc/amc_eqs/output/z2_terms.tex
    // A,B FoldI6I7 (T1); C,D FoldCD; E FoldE; F,I,K,L FoldI9(τ=Z+W); G,H FoldGH; J FoldJ.
    // m-scheme gold: ReferenceImplementations::evc_z2_mscheme.
    Operator BuildZ2RHS(const Operator &t, const Operator &z) const;
    // CC Z0: dZ_0/dλ = -R_0. AMC: run/evc/amc_eqs/output/z0_terms.tex
    double BuildZ0RHS(const Operator &t, const Operator &z) const;
    Operator SimilarityTransformCCSD(const Operator &tdagger, const Operator &z) const;
    RHSDiagrams BuildRHSDiagrams(const Operator &tdagger, const Operator &z) const;
    // dZ/dλ from locked CC residuals. T is excitation-only (T1,T2); T† ME are those of T.
    RHS BuildRHS(const Operator &t, const ClusterAmplitudes &amps) const;
    // Integrate dZ/dλ = -R(T†,Z) from λ=0 to λ_max. IC: Z0=0, Z1=T1, Z2=T2 (evc.tex).
    ClusterAmplitudes Solve(const Operator &t, double lambda_max = 1.0) const;
    ClusterAmplitudes SolveEuler(const Operator &t, const Operator &tdagger, double lambda_max = 1.0) const;
    ClusterAmplitudes SolveEulerFromOmega(const Operator &omega, double lambda_max = 1.0) const;
    ClusterAmplitudes SolveFromOmega(const Operator &omega, double lambda_max = 1.0) const;

    // CC connected vacuum pieces of a general H (0-body + OV + OOVV) with excitation Z.
    double CoupledClusterEnergy(const Operator &h, const Operator &z) const;
    // ⟨0|H e^Z|0⟩ = e^{Z_0} (H_0 + Δ_CC). H_0 additive from unitary BCH; Z_0 multiplicative from the ODE.
    double HamiltonianKernel(const Operator &h, const ClusterAmplitudes &amps) const;
    double NormKernel(const ClusterAmplitudes &amps) const;
    // H_ik = ⟨0|e^{-Ω_i} H e^{Ω_k}|0⟩, N_ik = ⟨0|e^{-Ω_i} e^{Ω_k}|0⟩.
    // H^i = e^{-Ω_i} H e^{Ω_i} (keep H^i_0); Ω_D = BCH_Product(-Ω_i, Ω_k) with Ω_D,0 erased;
    // e^{Ω_D}|0⟩ = e^Z|0⟩ from the ODE (Z_0 in amps.z0 only).
    std::array<double, 2> EvaluateKernels(const Operator &h0, const Operator &omega_i, const Operator &omega_k) const;

  private:
    ModelSpace *modelspace_;
    bool use_imsrg3_intermediates_ = true;
    std::size_t euler_steps_ = 100;
    bool use_rk4_ = true;

    Operator MakeExcitationOperatorLike(const Operator &op) const;
    void PromoteToThreeBody(Operator &op) const;
    void KeepExcitationBlocks(Operator &op, bool keep_zero_body) const;
    void KeepDeexcitationBlocks(Operator &op) const;
    void AddRHS(ClusterAmplitudes &amps, const RHS &rhs, double scale) const;
    ClusterAmplitudes Integrate(const Operator &t, double lambda_max, bool rk4) const;
    double ConnectedOV_OOVV(const Operator &op, const Operator &z) const;
};

#endif