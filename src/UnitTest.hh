#ifndef UnitTest_h
#define UnitTest_h

#include "ModelSpace.hh"
#include "Operator.hh"
#include "evc.hh"
//#include "Commutator.hh"

class UnitTest
{

 static uint64_t random_seed;

 public:

  ModelSpace* modelspace;

  UnitTest(ModelSpace&);

  void SetRandomSeed( uint64_t s ){ random_seed = s;};

  Operator RandomOp( ModelSpace& modelspace, int jrank, int tz, int parity, int particle_rank, int hermitian);

  Operator RandomDaggerOp( ModelSpace& modelspace, index_t Q);


  double GetMschemeMatrixElement_1b( const Operator& Op, int a, int ma, int b, int mb );
  double GetMschemeMatrixElement_2b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md );
  double GetMschemeMatrixElement_3b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md, int e, int me, int f, int mf );

  /// Gold A: unfactorized 3-operator m-scheme Γ^{IV_b}. No J-scheme χ.
  /// W = Σ_abcd w(c,d,b) [Ω_dibc Γ_acdk Ω_jbla − Ω_dkbc Γ_acdi Ω_jalb],
  /// w = n̄_c n_d n̄_b + n_c n̄_d n_b, then Z(m) = (1−P_ij)(1−P_kl) W.
  /// Ω may be tensor (WE-reduced) or scalar; Γ unreduced scalar.
  /// which_term: 0 both (W1−W2), 1 only W1, 2 only −W2. Then (1−P_ij)(1−P_kl).
  double Mscheme_comm223_232_GIVb(const Operator &Eta, const Operator &Gamma,
                                 int i, int mi, int j, int mj, int k, int mk, int l, int ml,
                                 int which_term = 0);

  // Factorized m-scheme gold (Python run/test_chi_*_mscheme.py + eq:fact folds).
  // Ω×Ω → scalar χ uses CG(λ μ, λ −μ; 00) as in those scripts. Ω×Γ χ has no extra CG.
  double Mscheme_chi_alpha(const Operator &Eta, int i, int mi, int j, int mj);
  double Mscheme_chi_beta(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  double Mscheme_chi_gamma(const Operator &Eta, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_delta(const Operator &Eta, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_epsilon(const Operator &Eta, int i, int mi, int j, int mj);
  double Mscheme_chi_zeta(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  /// G^{II} bra: χ^{ΩΓ}_{ij} = 1/2 Σ w Ω_ciab Γ_abcj. Same occ as χ^ζ. T×S, no CG inside.
  double Mscheme_chi_OmegaGamma(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  double Mscheme_chi_eta(const Operator &Eta, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_theta(const Operator &Eta, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_iota(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_kappa(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_chi_lambda(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);

  double Mscheme_fact_fI(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  double Mscheme_fact_fII(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  double Mscheme_fact_fIIIa(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);
  double Mscheme_fact_fIIIb(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj);

  double Mscheme_fact_GI(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GII(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GIIIa(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  /// G^{III_b} = IIb + IId as in comm223_232_BruteForce (distinct occ).
  /// IIb: −Σ (n̄_b n_c n_d + n_b n̄_c n̄_d) Ω_dcbk Ω_biac Γ_jald
  /// IId: −Σ (n̄_c n_b n_d + n_c n̄_b n̄_d) Ω_jcbd Ω_balc Γ_diak
  /// then Z = (1−P_ij)(1−P_kl) W. Not the shared-occ χ^η fold.
  double Mscheme_fact_GIIIb(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GIIIc(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GIVa(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GIVb_chi(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_fact_GIVc(const Operator &Eta, const Operator &Gamma, int i, int mi, int j, int mj, int k, int mk, int l, int ml);


  double GetMschemeMatrixElement_1leg( const Operator& Op, int a, int ma );
  double GetMschemeMatrixElement_3leg( const Operator& Op, int a, int ma, int b, int mb, int c, int mc );

//  void Test3BodyAntisymmetry();
  void Test3BodyAntisymmetry(Operator& Y);
  void Test3BodyHermiticity(Operator& Y);

  bool TestNormalOrdering(Operator& Op);

//  void Test3BodySetGet(Operator& Y);

  // test strategy: Fill two random operators, calculate a specific commutator term
  // using the J-coupled expression, and in m-scheme (where the formula is simpler)
  // and make sure that they give the same answer
//  void TestCommutators();
  bool TestCommutators(Operator& X, Operator& Y);
  bool TestCommutators_Tensor(Operator& X, Operator& Y);
//  bool TestCommutators_Tensor();
  bool TestCommutators_IsospinChanging();
  bool TestCommutators_ParityChanging();
//  void TestCommutators3();
//  void TestCommutators3(Operator& X, Operator& Y); 
  void TestCommutators3(Operator& X, Operator& Y, std::vector<std::string>& skiplist ); 

  void TestDaggerCommutators(index_t Q);
  void TestDaggerCommutatorsAlln(index_t Q);

  typedef void commutator_func (const Operator&,const Operator&,Operator&) ;
//  bool Test_against_ref_impl(const Operator& X, const Operator& Y, void (*ComRef)(const Operator&,const Operator&), void (*ComOpt)(const Operator&,const Operator&), std::string output_tag="" );
  bool Test_against_ref_impl(const Operator& X, const Operator& Y, commutator_func ComOpt, commutator_func ComRef, std::string output_tag="" );

  bool Mscheme_Test_comm110ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm220ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm111ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm121ss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm221ss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm122ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_pp_hhss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_phss( const Operator& X, const Operator& Y ) ;

  bool Mscheme_Test_comm121st( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm221st( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm222_pp_hhst( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_phst( const Operator& X, const Operator& Y ) ;

  bool Test_comm110ss( const Operator& X, const Operator& Y );
  bool Test_comm220ss( const Operator& X, const Operator& Y );
  bool Test_comm111ss( const Operator& X, const Operator& Y );
  bool Test_comm121ss( const Operator& X, const Operator& Y ); 
  bool Test_comm221ss( const Operator& X, const Operator& Y ); 
  bool Test_comm122ss( const Operator& X, const Operator& Y );
  bool Test_comm222_pp_hhss( const Operator& X, const Operator& Y );
  bool Test_comm222_phss( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_pp_hh_221ss( const Operator& X, const Operator& Y );


  // Tensor
  bool Test_comm111st( const Operator& X, const Operator& Y );
  bool Test_comm121st( const Operator& X, const Operator& Y );
  bool Test_comm122st( const Operator& X, const Operator& Y );
  bool Test_comm221st( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_pp_hhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_phst( const Operator& X, const Operator& Y ) ;
  // Tensor 3b
  bool Test_comm331st( const Operator& X, const Operator& Y ) ;
  bool Test_comm231st( const Operator& X, const Operator& Y ) ;
  bool Test_comm232st( const Operator& X, const Operator& Y ) ;
  bool Test_comm132st( const Operator& X, const Operator& Y ) ;
  bool Test_comm231tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm110tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm220tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm111tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm121tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm122tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm221tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_pp_hhtts( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_phtts( const Operator& X, const Operator& Y ) ;
  bool Test_comm132tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm232tts( const Operator& X, const Operator& Y ) ;
  bool Test_comm223tts( const Operator& X, const Operator& Y ) ;
  bool Test_scalar_tts_matches_ss( const Operator& X, const Operator& Y ) ;
  bool Test_comm223st( const Operator& X, const Operator& Y ) ;
  bool Test_comm133st( const Operator& X, const Operator& Y ) ;

  bool Test_comm332_pphhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm332_ppph_hhhpst( const Operator& X, const Operator& Y ) ;
  bool Test_comm233_pp_hhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm233_phst( const Operator& X, const Operator& Y ) ;
  bool Test_comm333_ppp_hhhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm333_pph_hhpst( const Operator& X, const Operator& Y ) ;

  bool Test_evc_rhs_ccsd( const Operator& Tdagger, const Operator& Z );
  bool Test_evc_z1_jscheme( const Operator& T, const Operator& Z );
  bool Test_evc_z2_jscheme( const Operator& T, const Operator& Z );
  bool Test_evc_z0_jscheme( const Operator& T, const Operator& Z );
  bool Test_evc_ode( const Operator& T );
  bool Test_evc_kernels( const Operator& H );

  bool Mscheme_Test_comm330ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm331ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm231ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm132ss( const Operator& X, const Operator& Y );
//
  bool Mscheme_Test_comm332_ppph_hhhpss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm332_pphhss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm133ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm223ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm233_pp_hhss( const Operator& X, const Operator& Y );    
  bool Mscheme_Test_comm233_phss( const Operator& X, const Operator& Y );    
  bool Mscheme_Test_comm232ss( const Operator& X, const Operator& Y );

  bool Mscheme_Test_comm333_ppp_hhhss( const Operator& X, const Operator& Y );  
  bool Mscheme_Test_comm333_pph_hhpss( const Operator& X, const Operator& Y );  

  // scalar-tensor commutator with 3b
  bool Mscheme_Test_comm331st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm223st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm231st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm232st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm232st_amc( const Operator& X, const Operator& Y );
  bool Test_comm232st_amc( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm133st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm132st( const Operator& X, const Operator& Y );
  /// Raw leftover m-scheme Wick (physical m). X is 2b tensor ω, Y is 3b tensor W=[ω,H]_3.
  /// Same strings as Mscheme_Test_comm231tts / 232tts / 132tts; no m-average, no CG projection.
  double Mscheme_comm231tts_wick(const Operator &X, const Operator &Y, int i, int mi, int j, int mj);
  double Mscheme_comm132tts_wick(const Operator &X, const Operator &Y, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  double Mscheme_comm232tts_wick(const Operator &X, const Operator &Y, int i, int mi, int j, int mj, int k, int mk, int l, int ml);
  /// Pure m-scheme Wick for [X,Y]_3 (=comm223tts kernel); scalar 3b m-component, no J-scheme.
  double Mscheme_comm223tts_wick(const Operator &X, const Operator &Y,
                                 int i, int mi, int j, int mj, int k, int mk,
                                 int l, int ml, int m, int mm, int n, int mn);
  bool Mscheme_Test_comm231tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm110tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm220tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm111tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm121tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm122tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm221tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_pp_hhtts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_phtts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm132tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm232tts( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm223tts( const Operator& X, const Operator& Y );
  /// Project m-scheme gold for comm232tts_bare Eq. eq onto J:
  /// S=Σ CG CG Z(m), reduced RME = S/Ĵ. PASS if GetTBME_J = −S/Ĵ on all stored TBMEs.
  bool RME_Test_comm232tts_bare( const Operator& X, const Operator& Y, int eq );
  bool Mscheme_Test_comm332_ppph_hhhpst(const Operator& X, const Operator& Y); 
  bool Mscheme_Test_comm332_pphhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm233_pp_hhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm233_phst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm333_ppp_hhhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm333_pph_hhpst(const Operator &X, const Operator &Y); 


  bool Test_comm330ss( const Operator& X, const Operator& Y );
  bool Test_comm331ss( const Operator& X, const Operator& Y );
  bool Test_comm231ss( const Operator& X, const Operator& Y );

  bool Test_comm132ss( const Operator& X, const Operator& Y );
  bool Test_comm232ss( const Operator& X, const Operator& Y );
  bool Test_comm332_ppph_hhhpss( const Operator& X, const Operator& Y ); 
  bool Test_comm332_pphhss( const Operator& X, const Operator& Y ); 

  bool Test_comm223ss( const Operator& X, const Operator& Y );
  bool Test_comm133ss( const Operator& X, const Operator& Y );

  bool Test_comm233_pp_hhss( const Operator& X, const Operator& Y );    
  bool Test_comm233_phss( const Operator& X, const Operator& Y );      
  bool Test_comm333_ppp_hhhss( const Operator& X, const Operator& Y );  
  bool Test_comm333_pph_hhpss( const Operator& X, const Operator& Y );  


  bool Test_comm211sd(        const Operator& X, const Operator& Y   );
  bool Test_comm231sd(        const Operator& X, const Operator& Y   );
  bool Test_comm431sd(        const Operator& X, const Operator& Y   );
  bool Test_comm233sd(        const Operator& X, const Operator& Y   );
  bool Test_comm413sd(        const Operator& X, const Operator& Y   );
  bool Test_comm433_pp_hh_sd( const Operator& X, const Operator& Y   );
  bool Test_comm433sd_ph(     const Operator& X, const Operator& Y   );

  bool TestRPAEffectiveCharge( const Operator& H, const Operator& OpIn, size_t k, size_t l);

//  bool TestFactorizedDoubleCommutators(ModelSpace& ms);
//  bool TestFactorizedDoubleCommutators();
  bool TestFactorizedDoubleCommutators(Operator& eta, Operator& H);

  bool TestPerturbativeTriples();

  bool SanityCheck();

 private:
  bool Mscheme_Test_comm232st_core(const Operator &X, const Operator &Y,
                                  commutator_func comm, const char *tag);

};


#endif
